use needletail::{Sequence, parse_fastx_file};
use num_traits::{NumCast, PrimInt, ToPrimitive};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;
use xxhash_rust::xxh3::xxh3_64_with_seed;

// DNA / base k-mer machinery
use kmerutils::base::{
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit, KmerBuilder,
    alphabet::Alphabet2b,
    kmergenerator::{KmerGenerationPattern, KmerGenerator},
    sequence::Sequence as SequenceStruct,
};
use kmerutils::sketcharg::{DataType, SeqSketcherParams, SketchAlgo};
use kmerutils::sketching::setsketchert::{OptDensHashSketch, RevOptDensHashSketch, SeqSketcherT};

// Amino-acid k-mers and AA sketchers (AA-specific generator + pattern)
use kmerutils::aautils::kmeraa::{
    Alphabet as AAAlphabet, KmerAA32bit, KmerAA64bit,
    KmerGenerationPattern as AAKmerGenerationPattern, KmerGenerator as AAKmerGenerator, SequenceAA,
};

use kmerutils::aautils::setsketchert as aasketch;
use kmerutils::aautils::setsketchert::SeqSketcherAAT;

use crate::metadata::PrefixParams;

/// Converts ASCII-encoded bases (from Needletail) into our `SequenceStruct`.
fn ascii_to_seq(bases: &[u8]) -> Result<SequenceStruct, ()> {
    let alphabet = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alphabet);
    Ok(seq)
}

/// Converts ASCII-encoded amino-acid sequence into `SequenceAA`.
/// - uppercases all residues
/// - drops any residue not in the 20-AA alphabet (ACDEFGHIKLMNPQRSTVWY),
///   including '*', 'X', gaps, etc.
fn ascii_to_seq_aa(bases: &[u8]) -> Result<SequenceAA, ()> {
    // Normalize to uppercase to match Alphabet's bases
    let mut upper = Vec::with_capacity(bases.len());
    for &b in bases {
        upper.push(b.to_ascii_uppercase());
    }

    let alphabet = AAAlphabet::new();
    Ok(SequenceAA::new_filtered(&upper, &alphabet))
}

/// Read lines (one path per line) from a text file.
/// Strict, order-preserving parallel sketching:
/// returns Vec<Vec<u16>> aligned to `file_paths` order.
///
/// IMPORTANT: we REQUIRE Sketcher::Sig = u16 at compile time (for DNA).
fn sketch_files_ordered_u16<Kmer, Sketcher, F>(
    file_paths: &[String],
    sketcher: &Sketcher,
    kmer_hash_fn: F,
) -> Vec<Vec<u16>>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + std::fmt::Debug,
    KmerGenerator<Kmer>: KmerGenerationPattern<Kmer>,
    Sketcher: SeqSketcherT<Kmer, Sig = u16> + Sync,
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    let completed = AtomicUsize::new(0);
    let total = file_paths.len();
    // We keep a local Vec<(original_index, signature)> and then sort by index
    // to guarantee that final vectors follow file_paths order exactly.
    let mut indexed: Vec<(usize, Vec<u16>)> = file_paths
        .par_iter()
        .enumerate()
        .map(|(i, path)| {
            let mut sequences: Vec<SequenceStruct> = Vec::new();
            let mut reader =
                parse_fastx_file(path).unwrap_or_else(|e| panic!("Invalid FASTA/Q {path}: {e}"));

            while let Some(record) = reader.next() {
                let rec = record.unwrap_or_else(|e| panic!("Error reading record in {path}: {e}"));
                let seq_norm = rec.normalize(false).into_owned();
                let seq = ascii_to_seq(&seq_norm).unwrap();
                sequences.push(seq);
            }

            let sequences_ref: Vec<&SequenceStruct> = sequences.iter().collect();

            // returns Vec<Vec<u16>>, and for "seqs" interface inner vec has size 1
            let signature = sketcher.sketch_compressedkmer_seqs(&sequences_ref, kmer_hash_fn);
            let sig_u16 = signature
                .first()
                .cloned()
                .unwrap_or_else(|| panic!("sketcher returned empty signature for {path}"));

            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Sketched {done}/{total} genomes");
            }

            (i, sig_u16)
        })
        .collect();

    indexed.sort_by_key(|(i, _)| *i);
    indexed.into_iter().map(|(_, v)| v).collect()
}

/// Strict, order-preserving parallel AA sketching (Sig = u16).
/// Returns Vec<Vec<u16>> aligned to `file_paths` order.
fn sketch_files_ordered_aa_u16<Kmer, Sketcher, F>(
    file_paths: &[String],
    sketcher: &Sketcher,
    kmer_hash_fn: F,
) -> Vec<Vec<u16>>
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Send + Sync,
    <Kmer as CompressedKmerT>::Val: num::PrimInt + Send + Sync + std::fmt::Debug,
    // *** AA-specific generator + pattern bound (this fixes E0277) ***
    AAKmerGenerator<Kmer>: AAKmerGenerationPattern<Kmer>,
    Sketcher: SeqSketcherAAT<Kmer, Sig = u16> + Sync,
    F: Fn(&Kmer) -> <Kmer as CompressedKmerT>::Val + Send + Sync + Copy,
{
    let completed = AtomicUsize::new(0);
    let total = file_paths.len();
    let mut indexed: Vec<(usize, Vec<u16>)> = file_paths
        .par_iter()
        .enumerate()
        .map(|(i, path)| {
            let mut sequences: Vec<SequenceAA> = Vec::new();
            let mut reader =
                parse_fastx_file(path).unwrap_or_else(|e| panic!("Invalid FASTA/Q {path}: {e}"));

            let mut record_num = 0_u64;
            while let Some(record) = reader.next() {
                let rec = record.unwrap_or_else(|e| panic!("Error reading record in {path}: {e}"));
                // For AA we do NOT normalize as DNA; use raw sequence bytes.
                let seq_bytes = rec.seq();
                let sequence_id = String::from_utf8_lossy(rec.id());
                if seq_bytes.is_empty() {
                    eprintln!(
                        "ERROR: sequence of null length, file: {path:?}, record num: {record_num}, sequence id: {sequence_id}"
                    );
                    record_num += 1;
                    continue;
                }

                let seq = ascii_to_seq_aa(&seq_bytes)
                    .unwrap_or_else(|_| panic!("AA parse error in {path}"));
                if seq.is_empty() {
                    eprintln!(
                        "ERROR: null encoded sequence, file: {path:?}, record num: {record_num}, sequence id: {sequence_id}"
                    );
                    record_num += 1;
                    continue;
                }
                sequences.push(seq);
                record_num += 1;
            }

            let sequences_ref: Vec<&SequenceAA> = sequences.iter().collect();

            // AA sketcher interface: sketch_compressedkmeraa_seqs
            let signature = sketcher.sketch_compressedkmeraa_seqs(&sequences_ref, kmer_hash_fn);
            let sig_u16 = signature
                .first()
                .cloned()
                .unwrap_or_else(|| panic!("sketcher returned empty signature for {path}"));

            let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if done % 1000 == 0 {
                eprintln!("Sketched {done}/{total} genomes");
            }

            (i, sig_u16)
        })
        .collect();

    indexed.sort_by_key(|(i, _)| *i);
    indexed.into_iter().map(|(_, v)| v).collect()
}

/// Canonicalize (min(kmer, revcomp)) -> pack -> XXH3 -> return Kmer::Val (u32/u64).
/// Safe for k<=32 (never shifts by 64).
pub fn make_xxh3_canonical_kmer_hash_fn<Kmer>(
    seed: u64,
) -> impl Fn(&Kmer) -> Kmer::Val + Copy + Send + Sync
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer> + Copy,
    Kmer::Val: PrimInt + ToPrimitive + NumCast,
{
    move |kmer: &Kmer| -> Kmer::Val {
        // canonicalize
        let rc = kmer.reverse_complement();
        let canonical = if rc < *kmer { rc } else { *kmer };

        // DNA packing: 2 bits per base
        let k: usize = canonical.get_nb_base() as usize;
        let bits: usize = 2usize * k;

        // safe mask (avoid 1<<64)
        let mask_u64: u64 = if bits >= 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };

        // convert packed bits to u64 safely (works for u32/u64)
        let packed_u64: u64 = canonical
            .get_compressed_value()
            .to_u64()
            .expect("Kmer::Val must fit in u64")
            & mask_u64;

        // hash packed bytes
        let h64: u64 = xxh3_64_with_seed(&packed_u64.to_le_bytes(), seed);

        // output type:
        // - if Val is u32-ish => use low 32 bits
        // - if Val is u64-ish => keep full 64
        if std::mem::size_of::<Kmer::Val>() <= 4 {
            let low32 = h64 as u32 as u64; // intentional truncation
            NumCast::from(low32).expect("casting hash->Val failed")
        } else {
            NumCast::from(h64).expect("casting hash->Val failed")
        }
    }
}

/// Amino-acid packed k-mer hashing (no reverse-complement).
/// AA alphabet needs 5 bits (2^5 = 32) for 20 residues.
pub fn make_xxh3_aa_kmer_hash_fn<Kmer>(
    seed: u64,
) -> impl Fn(&Kmer) -> Kmer::Val + Copy + Send + Sync
where
    Kmer: CompressedKmerT + KmerBuilder<Kmer>,
    Kmer::Val: PrimInt + NumCast + ToPrimitive,
{
    const AA_BITS_PER_RESIDUE: usize = 5;

    move |kmer: &Kmer| -> Kmer::Val {
        let k: usize = kmer.get_nb_base() as usize;
        let bits: usize = AA_BITS_PER_RESIDUE * k;

        let mask_u64: u64 = if bits >= 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };

        let packed_u64: u64 = kmer
            .get_compressed_value()
            .to_u64()
            .expect("Kmer::Val must be convertible to u64")
            & mask_u64;

        let h64: u64 = xxh3_64_with_seed(&packed_u64.to_le_bytes(), seed);

        let out_u64 = if std::mem::size_of::<Kmer::Val>() <= 4 {
            (h64 as u32) as u64
        } else {
            h64
        };

        NumCast::from(out_u64).expect("cast to Kmer::Val must succeed")
    }
}

/// Dispatch sketching by k-mer size.
/// Fixed Sig=u16 for the sketcher output (DNA).
fn sketch_with_kmer_dispatch_u16(
    paths: &[String],
    kmer_size: usize,
    sketch_size: usize,
    densification: usize, // 0 optdens, 1 revoptdens
    hash_seed: u64,
) -> Vec<Vec<u16>> {
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS, // label; actual is controlled by sketcher type
        DataType::DNA,
    );

    if kmer_size <= 14 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size == 16 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer16b32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer16b32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer16b32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer16b32bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer16b32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size <= 32 {
        let kmer_hash_fn = make_xxh3_canonical_kmer_hash_fn::<Kmer64bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = OptDensHashSketch::<Kmer64bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer64bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            1 => {
                let sketcher = RevOptDensHashSketch::<Kmer64bit, f32>::new(&sketch_args);
                sketch_files_ordered_u16::<Kmer64bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else {
        panic!("kmer_size cannot exceed 32 and must not be 15!");
    }
}

/// Amino-acid k-mer dispatcher (Sig = u16).
/// AA k must be <= 12 (KmerAA32bit for k<=6, KmerAA64bit for 7..=12).
fn sketch_aa_with_kmer_dispatch_u16(
    paths: &[String],
    kmer_size: usize,
    sketch_size: usize,
    densification: usize, // 0 optdens, 1 revoptdens
    hash_seed: u64,
) -> Vec<Vec<u16>> {
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS, // label; actual controlled by sketcher type
        DataType::AA,
    );

    if kmer_size <= 6 {
        let kmer_hash_fn = make_xxh3_aa_kmer_hash_fn::<KmerAA32bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = aasketch::OptDensHashSketch::<KmerAA32bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            1 => {
                let sketcher =
                    aasketch::RevOptDensHashSketch::<KmerAA32bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA32bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else if kmer_size <= 12 {
        let kmer_hash_fn = make_xxh3_aa_kmer_hash_fn::<KmerAA64bit>(hash_seed);

        match densification {
            0 => {
                let sketcher = aasketch::OptDensHashSketch::<KmerAA64bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA64bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            1 => {
                let sketcher =
                    aasketch::RevOptDensHashSketch::<KmerAA64bit, f32>::new(&sketch_args);
                sketch_files_ordered_aa_u16::<KmerAA64bit, _, _>(paths, &sketcher, kmer_hash_fn)
            }
            _ => panic!("densification must be 0 or 1"),
        }
    } else {
        panic!("kmer_size for amino acids must be <= 12");
    }
}

pub(crate) fn sketch_from_params(paths: &[String], params: &PrefixParams) -> Vec<Vec<u16>> {
    let started = Instant::now();
    log::debug!(
        "sketching files={} seq_type={} kmer_size={} sketch_size={} densification={}",
        paths.len(),
        params.seq_type,
        params.kmer_size,
        params.sketch_size,
        params.densification
    );
    let sketches = match params.seq_type.as_str() {
        "dna" => sketch_with_kmer_dispatch_u16(
            paths,
            params.kmer_size,
            params.sketch_size,
            params.densification,
            params.hash_seed,
        ),
        "aa" => sketch_aa_with_kmer_dispatch_u16(
            paths,
            params.kmer_size,
            params.sketch_size,
            params.densification,
            params.hash_seed,
        ),
        other => panic!("unknown sequence type in parameters: {other}"),
    };
    eprintln!(
        "Sketched {}/{} genomes in {:.3}s",
        paths.len(),
        paths.len(),
        started.elapsed().as_secs_f64()
    );
    sketches
}

// Database operations, CLI parsing, and process orchestration live in their
// own modules. This module only owns sequence parsing and sketch generation.

#[cfg(test)]
mod tests {
    use super::ascii_to_seq_aa;

    #[test]
    fn aa_conversion_exposes_empty_and_fully_filtered_records() {
        assert!(ascii_to_seq_aa(b"").unwrap().is_empty());
        assert!(ascii_to_seq_aa(b"***XXX---").unwrap().is_empty());
        assert!(!ascii_to_seq_aa(b"MTEQIELIK").unwrap().is_empty());
    }
}
