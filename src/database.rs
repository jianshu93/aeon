use anndists::dist::DistHamming;
use log::debug;
use rayon::prelude::*;
use rust_diskann::{DiskANN, DiskAnnParams};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::metadata::{
    PrefixParams, genomes_path, idmap_path, index_path, load_params, save_params,
};
use crate::sketch::sketch_from_params;
use crate::utils::{read_list_file, write_genome_list, write_idmap_tsv};

pub(crate) fn build_database(
    prefix: &str,
    reference_list: &str,
    params: PrefixParams,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let genomes = read_list_file(reference_list)?;
    if genomes.is_empty() {
        return Err("reference list is empty".into());
    }

    eprintln!(
        "Building database: prefix={}, genomes={}, seq_type={}, k={}, sketch_size={}, R={}, build_beam={}",
        prefix,
        genomes.len(),
        params.seq_type,
        params.kmer_size,
        params.sketch_size,
        params.max_degree,
        params.build_beam_width
    );
    let mut unique = HashSet::with_capacity(genomes.len());
    for genome in &genomes {
        if !unique.insert(genome.as_str()) {
            return Err(format!("reference list contains duplicate path: {genome}").into());
        }
        if !Path::new(genome).is_file() {
            return Err(
                format!("reference genome does not exist or is not a file: {genome}").into(),
            );
        }
    }

    let vectors = sketch_from_params(&genomes, &params);
    let diskann_params = DiskAnnParams {
        max_degree: params.max_degree,
        build_beam_width: params.build_beam_width,
        alpha: params.alpha,
        extra_seeds: params.extra_seeds,
    };
    let index = DiskANN::<u16, DistHamming>::build_index_with_params(
        &vectors,
        DistHamming,
        &index_path(prefix),
        diskann_params,
    )?;
    sanity_check_index_mapping(&index, &vectors, 3, params.build_beam_width);
    write_genome_list(&genomes_path(prefix), &genomes)?;
    write_idmap_tsv(&idmap_path(prefix), &genomes)?;
    save_params(prefix, &params)?;
    eprintln!(
        "Built database: prefix={}, genomes={}, elapsed={:.3}s",
        prefix,
        genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

pub(crate) fn update_database(
    prefix: &str,
    insert_paths: Vec<String>,
    delete_names: Vec<String>,
    beam_width: Option<usize>,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    if insert_paths.is_empty() && delete_names.is_empty() {
        return Err("update contains neither insertions nor deletions".into());
    }
    let params = load_params(prefix)?;
    let old_genomes = read_list_file(&genomes_path(prefix))?;
    let delete_ids = resolve_delete_ids(&old_genomes, &delete_names)?;
    let deleted_names = delete_names
        .iter()
        .map(String::as_str)
        .collect::<HashSet<_>>();
    let retained_names = old_genomes
        .iter()
        .filter(|name| !deleted_names.contains(name.as_str()))
        .map(String::as_str)
        .collect::<HashSet<_>>();
    let mut inserted = HashSet::with_capacity(insert_paths.len());
    for path in &insert_paths {
        if !Path::new(path).is_file() {
            return Err(format!("insert genome does not exist or is not a file: {path}").into());
        }
        if !inserted.insert(path.as_str()) {
            return Err(format!("insert list contains duplicate path: {path}").into());
        }
        if retained_names.contains(path.as_str()) {
            return Err(format!("genome is already present and not being deleted: {path}").into());
        }
    }

    let insertion_count = insert_paths.len();
    let insert_vectors = sketch_from_params(&insert_paths, &params);
    let source = DiskANN::<u16, DistHamming>::open_index_with(&index_path(prefix), DistHamming)?;
    if source.num_vectors != old_genomes.len() {
        return Err(format!(
            "index has {} vectors but genome map has {} entries",
            source.num_vectors,
            old_genomes.len()
        )
        .into());
    }
    let final_count = old_genomes.len() - delete_ids.len() + insert_paths.len();
    if final_count == 0 {
        return Err("an update cannot produce an empty database".into());
    }
    let capacity = old_genomes.len().max(final_count);
    let work_path = format!("{prefix}.aeon-work-{}", std::process::id());
    let temporary_index = format!("{}.aeon-next-{}", index_path(prefix), std::process::id());
    eprintln!(
        "Updating database: prefix={}, current={}, delete={}, insert={}, final={}",
        prefix,
        old_genomes.len(),
        delete_ids.len(),
        insertion_count,
        final_count
    );
    let mut update = source.begin_updates(capacity, params.alpha, &work_path)?;
    let mut slot_names = old_genomes.into_iter().map(Some).collect::<Vec<_>>();
    slot_names.resize(capacity, None);

    if !delete_ids.is_empty() {
        let stats = update.delete_batch(&delete_ids)?;
        for id in &delete_ids {
            slot_names[*id as usize] = None;
        }
        let recovered: usize = stats.iter().map(|stat| stat.recovered_in_neighbors).sum();
        let candidates: usize = stats.iter().map(|stat| stat.repair_candidates).sum();
        debug!(
            "Deleted {} genomes; MERIT averages: {:.2} recovered in-neighbors, {:.2} repair candidates",
            stats.len(),
            recovered as f64 / stats.len() as f64,
            candidates as f64 / stats.len() as f64
        );
        eprintln!("Deleted {} genomes", stats.len());
    }
    if !insert_vectors.is_empty() {
        let beam = beam_width.unwrap_or(params.build_beam_width);
        let inserted_ids = update.insert_batch(insert_vectors, beam)?;
        for (id, name) in inserted_ids.into_iter().zip(insert_paths) {
            slot_names[id as usize] = Some(name);
        }
        eprintln!("Inserted {insertion_count} genomes with beam {beam}");
    }

    let (committed, old_to_new) = update.commit_updates_to_static(&temporary_index)?;
    if committed.num_vectors != final_count {
        return Err(format!(
            "committed index has {} vectors; expected {final_count}",
            committed.num_vectors
        )
        .into());
    }
    let new_genomes = rebuild_genome_order(&slot_names, &old_to_new, final_count)?;
    drop(committed);
    install_updated_files(prefix, &temporary_index, &new_genomes)?;
    eprintln!(
        "Updated database: prefix={}, genomes={}, elapsed={:.3}s",
        prefix,
        new_genomes.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

pub(crate) fn search_database(
    prefix: &str,
    query_list: &str,
    k: usize,
    beam: usize,
    output_path: Option<&str>,
) -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let params = load_params(prefix)?;
    let references = read_list_file(&genomes_path(prefix))?;
    let queries = read_list_file(query_list)?;
    if queries.is_empty() {
        return Err("query list is empty".into());
    }
    eprintln!(
        "Searching database: prefix={}, queries={}, references={}, k={}, beam={}",
        prefix,
        queries.len(),
        references.len(),
        k,
        beam
    );
    let vectors = sketch_from_params(&queries, &params);
    let index = DiskANN::<u16, DistHamming>::open_index_with(&index_path(prefix), DistHamming)?;
    if index.num_vectors != references.len() {
        return Err("index and genome mapping have different lengths".into());
    }
    let hits = vectors
        .par_iter()
        .enumerate()
        .map(|(query, vector)| (query, index.search_with_dists(vector, k, beam)))
        .collect::<Vec<_>>();
    let output: Box<dyn Write> = match output_path {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    write_search_results(output, &queries, &hits, &references)?;
    eprintln!(
        "Searched {} queries in {:.3}s",
        queries.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn resolve_delete_ids(genomes: &[String], names: &[String]) -> Result<Vec<u32>, Box<dyn Error>> {
    let mut by_name = HashMap::with_capacity(genomes.len());
    for (id, name) in genomes.iter().enumerate() {
        if by_name.insert(name.as_str(), id as u32).is_some() {
            return Err(format!("database contains duplicate genome name: {name}").into());
        }
    }
    let mut seen = HashSet::with_capacity(names.len());
    names
        .iter()
        .map(|name| {
            if !seen.insert(name.as_str()) {
                return Err(format!("delete list contains duplicate name: {name}").into());
            }
            by_name
                .get(name.as_str())
                .copied()
                .ok_or_else(|| format!("genome name is not present in the database: {name}").into())
        })
        .collect()
}

fn rebuild_genome_order(
    slot_names: &[Option<String>],
    old_to_new: &[Option<u32>],
    expected: usize,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut genomes = vec![None; expected];
    for (old_id, new_id) in old_to_new.iter().enumerate() {
        if let Some(new_id) = new_id {
            let name = slot_names
                .get(old_id)
                .and_then(Option::as_ref)
                .ok_or_else(|| format!("missing genome name for live slot {old_id}"))?;
            genomes[*new_id as usize] = Some(name.clone());
        }
    }
    genomes
        .into_iter()
        .enumerate()
        .map(|(id, name)| name.ok_or_else(|| format!("missing genome name for new ID {id}").into()))
        .collect()
}

fn install_updated_files(
    prefix: &str,
    temporary_index: &str,
    genomes: &[String],
) -> Result<(), Box<dyn Error>> {
    let suffix = format!("aeon-next-{}", std::process::id());
    let temporary_genomes = format!("{}.{}", genomes_path(prefix), suffix);
    let temporary_idmap = format!("{}.{}", idmap_path(prefix), suffix);
    write_genome_list(&temporary_genomes, genomes)?;
    write_idmap_tsv(&temporary_idmap, genomes)?;
    fs::rename(&temporary_genomes, genomes_path(prefix))?;
    fs::rename(&temporary_idmap, idmap_path(prefix))?;
    fs::rename(temporary_index, index_path(prefix))?;
    Ok(())
}

fn write_search_results(
    mut output: Box<dyn Write>,
    query_paths: &[String],
    hits: &[(usize, Vec<(u32, f32)>)],
    references: &[String],
) -> io::Result<()> {
    writeln!(output, "Query\tHit\tVectorID\tRawJaccard\tHammingDist")?;
    for (query, neighbors) in hits {
        for (id, distance) in neighbors {
            let raw_jaccard = (1.0 - *distance).clamp(0.0, 1.0);
            writeln!(
                output,
                "{}\t{}\t{}\t{:.6}\t{:.6}",
                query_paths[*query], references[*id as usize], id, raw_jaccard, distance
            )?;
        }
    }
    Ok(())
}

fn sanity_check_index_mapping(
    index: &DiskANN<u16, DistHamming>,
    vectors: &[Vec<u16>],
    num_checks: usize,
    build_beam_width: usize,
) {
    let checks = vectors.len().min(num_checks.max(1));
    let beam = build_beam_width.max(64);
    for (id, vector) in vectors.iter().take(checks).enumerate() {
        match index.search_with_dists(vector, 1, beam).first().copied() {
            Some((found, distance)) => debug!(
                "index mapping check reference_id={id} nearest_id={found} distance={distance:.6}"
            ),
            None => debug!("index mapping check reference_id={id} returned no neighbors"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delete_names_require_exact_matches() {
        let genomes = vec!["/db/a.fna.gz".into(), "/db/b.fna.gz".into()];
        assert_eq!(
            resolve_delete_ids(&genomes, &["/db/b.fna.gz".into()]).unwrap(),
            vec![1]
        );
        assert!(resolve_delete_ids(&genomes, &["b.fna.gz".into()]).is_err());
    }

    #[test]
    fn rebuilds_names_from_compacted_ids() {
        let slots = vec![Some("a".into()), None, Some("c".into()), Some("d".into())];
        let mapping = vec![Some(0), None, Some(1), Some(2)];
        assert_eq!(
            rebuild_genome_order(&slots, &mapping, 3).unwrap(),
            vec!["a", "c", "d"]
        );
    }
}
