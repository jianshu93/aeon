use rayon::ThreadPoolBuilder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::sync::Once;

static INIT_RAYON: Once = Once::new();

pub(crate) fn init_rayon_global(threads: usize) {
    INIT_RAYON.call_once(|| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .expect("failed to initialize the global Rayon thread pool");
    });
    log::debug!("Rayon thread pool ready threads={threads}");
}

pub(crate) fn requested_threads(matches: &clap::ArgMatches) -> usize {
    matches
        .get_one::<usize>("threads")
        .copied()
        .unwrap_or_else(num_cpus::get)
}

pub(crate) fn read_list_file(filepath: &str) -> io::Result<Vec<String>> {
    let reader = BufReader::new(File::open(filepath)?);
    reader
        .lines()
        .map(|line| line.map(|value| value.trim().to_owned()))
        .filter(|line| line.as_ref().map_or(true, |value| !value.is_empty()))
        .collect()
}

pub(crate) fn write_genome_list(filepath: &str, genomes: &[String]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(filepath)?);
    for genome in genomes {
        writeln!(writer, "{genome}")?;
    }
    Ok(())
}

pub(crate) fn write_idmap_tsv(filepath: &str, genomes: &[String]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(filepath)?);
    for (id, genome) in genomes.iter().enumerate() {
        writeln!(writer, "{id}\t{genome}")?;
    }
    Ok(())
}
