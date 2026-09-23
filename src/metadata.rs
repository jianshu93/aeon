use serde::{Deserialize, Serialize};
use std::{fs, io};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct PrefixParams {
    pub(crate) kmer_size: usize,
    pub(crate) sketch_size: usize,
    pub(crate) densification: usize,
    pub(crate) seq_type: String,
    pub(crate) max_degree: usize,
    pub(crate) build_beam_width: usize,
    pub(crate) alpha: f32,
    pub(crate) extra_seeds: usize,
    pub(crate) sketch_elem_size: u8,
    pub(crate) sketch_elem_type: String,
    pub(crate) distance: String,
    pub(crate) hash_seed: u64,
}

pub(crate) fn params_path(prefix: &str) -> String {
    format!("{prefix}.params.json")
}

pub(crate) fn index_path(prefix: &str) -> String {
    format!("{prefix}.diskann")
}

pub(crate) fn genomes_path(prefix: &str) -> String {
    format!("{prefix}.genomes.txt")
}

pub(crate) fn idmap_path(prefix: &str) -> String {
    format!("{prefix}.idmap.tsv")
}

pub(crate) fn save_params(prefix: &str, params: &PrefixParams) -> io::Result<()> {
    let json = serde_json::to_string_pretty(params).expect("PrefixParams must serialize");
    fs::write(params_path(prefix), json)
}

pub(crate) fn load_params(prefix: &str) -> io::Result<PrefixParams> {
    let json = fs::read_to_string(params_path(prefix))?;
    serde_json::from_str(&json).map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))
}
