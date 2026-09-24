mod database;
mod metadata;
mod sketch;
mod utils;

use std::error::Error;

use crate::database::{build_database, search_database, update_database};
use crate::metadata::PrefixParams;
use crate::utils::{init_rayon_global, read_list_file, requested_threads};

use clap::{Arg, ArgAction, Command};

fn prefix_arg() -> Arg {
    Arg::new("prefix")
        .long("prefix")
        .short('p')
        .required(true)
        .action(ArgAction::Set)
        .value_name("PREFIX")
        .help("Database prefix shared by the .diskann, .genomes.txt, .idmap.tsv, and .params.json files")
}

fn threads_arg() -> Arg {
    Arg::new("threads")
        .long("threads")
        .short('t')
        .value_name("N")
        .value_parser(clap::value_parser!(usize))
        .help("Worker threads [default: all logical CPUs]")
}

fn update_beam_arg() -> Arg {
    Arg::new("beam_width")
        .long("beam-width")
        .value_name("N")
        .value_parser(clap::value_parser!(usize))
        .help("Candidate beam for inserted genomes [default: original build beam]")
        .long_help("Candidate beam used to connect inserted genomes to the existing Vamana graph. Larger values generally improve graph quality but increase update time. When omitted, Aeon reuses the build beam stored in PREFIX.params.json.")
}

fn build_cli() -> Command {
    Command::new("aeon")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Efficient and robust genome similarity search with evolving disk-based proximity graphs")
        .long_about("Build and search mmap-backed genome sketch indexes. Insertions and MERIT deletions run in a temporary dynamic workspace and are committed back to the ordinary static DiskANN format.")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("todisk")
                .about("Build a new static diskANN index using the Vamana graph construction algorithm")
                .long_about("Sketch every sequence file in the reference list, construct a Vamana graph, and write the static mmap-searchable index plus its genome-name mapping and parameter metadata.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("reference_list")
                        .long("reference-list")
                        .short('r')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing one reference FASTA/FASTQ path per line")
                        .long_help("Text file containing one reference FASTA/FASTQ path per line. Needletail-supported compressed inputs are accepted. Paths are retained as the genome names used by search and exact deletion."),
                )
                .arg(
                    Arg::new("kmer_size")
                        .long("kmer-size")
                        .short('k')
                        .default_value("16")
                        .value_name("K")
                        .value_parser(clap::value_parser!(usize))
                        .help("k-mer length; DNA supports <=32 except 15, amino acid supports <=12"),
                )
                .arg(
                    Arg::new("sketch_size")
                        .long("sketch-size")
                        .short('s')
                        .default_value("8192")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of u16 values in each MinHash sketch"),
                )
                .arg(
                    Arg::new("densification")
                        .long("densification")
                        .short('d')
                        .default_value("0")
                        .value_name("MODE")
                        .value_parser(["0", "1"])
                        .help("Densification algorithm: 0=OPTDENS, 1=RevOPTDENS"),
                )
                .arg(
                    Arg::new("seq_type")
                        .long("seq-type")
                        .default_value("dna")
                        .value_name("TYPE")
                        .value_parser(["dna", "aa"])
                        .help("Input alphabet: dna or aa (amino acid)"),
                )
                .arg(
                    Arg::new("max_degree")
                        .long("max-degree")
                        .default_value("64")
                        .value_name("R")
                        .value_parser(clap::value_parser!(usize))
                        .help("Maximum Vamana out-degree; larger values use more space and usually improve recall"),
                )
                .arg(
                    Arg::new("build_beam_width")
                        .long("build-beam-width")
                        .default_value("256")
                        .value_name("L")
                        .value_parser(clap::value_parser!(usize))
                        .help("Construction beam; larger values improve graph quality at higher build cost"),
                )
                .arg(
                    Arg::new("alpha")
                        .long("alpha")
                        .default_value("1.2")
                        .value_name("FLOAT")
                        .value_parser(clap::value_parser!(f32))
                        .help("Vamana robust-pruning alpha parameter"),
                )
                .arg(
                    Arg::new("extra_seeds")
                        .long("extra-seeds")
                        .default_value("1")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Additional random construction-search entry points per node"),
                )
                .arg(
                    Arg::new("hash_seed")
                        .long("hash-seed")
                        .default_value("1337")
                        .value_name("INTEGER")
                        .value_parser(clap::value_parser!(u64))
                        .help("XXH3 seed stored in metadata and reused for every future sketch"),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("insert")
                .about("Insert sketch vectors into an existing diskANN index after sketching new genomes")
                .long_about("Sketch new genome files with the original database parameters, connect them to a temporary dynamic Vamana graph, and commit a compact static index.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("genome_list")
                        .long("genome-list")
                        .short('i')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing new genome FASTA/FASTQ file paths to sketch and then insert, one per line"),
                )
                .arg(update_beam_arg())
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("delete")
                .about("Delete sketch vectors in the diskANN index with MERIT graph repair and version invalidation")
                .long_about("Delete sketch vectors with the exact stored genome names, repair affected graph neighborhoods with MERIT, compact vector IDs, and commit an ordinary static index. Sequence files are not read during deletion.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("name_list")
                        .long("name-list")
                        .short('d')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing exact names from PREFIX.genomes.txt")
                        .long_help("Text file containing one stored genome name per line. Matching is exact, including directory components; basename-only matching is intentionally not performed."),
                )
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("update")
                .about("Delete sketch vectors in vamana graph and then insert new sketch vectors in one transaction")
                .long_about("Apply MERIT deletion and Vamana insertion in one temporary update session, then perform one compact static commit. A deleted path may be reinserted in the same command.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("delete_list")
                        .long("delete-list")
                        .short('d')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing exact stored genome names to delete"),
                )
                .arg(
                    Arg::new("insert_list")
                        .long("insert-list")
                        .short('i')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing genome FASTA/FASTQ paths to sketch and then insert"),
                )
                .arg(update_beam_arg())
                .arg(threads_arg()),
        )
        .subcommand(
            Command::new("search")
                .about("Search query genomes against a static diskANN index")
                .long_about("Sketch query genomes using parameters recovered from database metadata, search the mmap-backed static DiskANN graph in parallel, and report ranked neighbors.")
                .arg(prefix_arg())
                .arg(
                    Arg::new("query_list")
                        .long("query-list")
                        .short('q')
                        .required(true)
                        .value_name("FILE")
                        .help("Text file containing one query FASTA/FASTQ path per line"),
                )
                .arg(
                    Arg::new("k")
                        .long("k")
                        .default_value("10")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Number of nearest neighbors reported per query"),
                )
                .arg(
                    Arg::new("beam_width")
                        .long("beam-width")
                        .default_value("512")
                        .value_name("N")
                        .value_parser(clap::value_parser!(usize))
                        .help("Search beam; larger values usually improve recall but reduce throughput"),
                )
                .arg(
                    Arg::new("output")
                        .long("output")
                        .short('o')
                        .value_name("FILE")
                        .help("Write TSV results to FILE instead of standard output"),
                )
                .arg(threads_arg()),
        )
}

fn main() -> Result<(), Box<dyn Error>> {
    let _ = env_logger::Builder::from_default_env().try_init();
    let matches = build_cli().get_matches();
    let (command, args) = matches.subcommand().expect("subcommand is required");
    init_rayon_global(requested_threads(args));
    log::info!("starting command={command}");

    match command {
        "todisk" => build_database(
            args.get_one::<String>("prefix").unwrap(),
            args.get_one::<String>("reference_list").unwrap(),
            PrefixParams {
                kmer_size: *args.get_one("kmer_size").unwrap(),
                sketch_size: *args.get_one("sketch_size").unwrap(),
                densification: args
                    .get_one::<String>("densification")
                    .unwrap()
                    .parse()
                    .unwrap(),
                seq_type: args.get_one::<String>("seq_type").unwrap().clone(),
                max_degree: *args.get_one("max_degree").unwrap(),
                build_beam_width: *args.get_one("build_beam_width").unwrap(),
                alpha: *args.get_one("alpha").unwrap(),
                extra_seeds: *args.get_one("extra_seeds").unwrap(),
                sketch_elem_size: 2,
                sketch_elem_type: "u16".into(),
                distance: "anndists::dist::DistHamming".into(),
                hash_seed: *args.get_one("hash_seed").unwrap(),
            },
        )?,
        "insert" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            read_list_file(args.get_one::<String>("genome_list").unwrap())?,
            Vec::new(),
            args.get_one::<usize>("beam_width").copied(),
        )?,
        "delete" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            Vec::new(),
            read_list_file(args.get_one::<String>("name_list").unwrap())?,
            None,
        )?,
        "update" => update_database(
            args.get_one::<String>("prefix").unwrap(),
            read_list_file(args.get_one::<String>("insert_list").unwrap())?,
            read_list_file(args.get_one::<String>("delete_list").unwrap())?,
            args.get_one::<usize>("beam_width").copied(),
        )?,
        "search" => search_database(
            args.get_one::<String>("prefix").unwrap(),
            args.get_one::<String>("query_list").unwrap(),
            *args.get_one::<usize>("k").unwrap(),
            *args.get_one::<usize>("beam_width").unwrap(),
            args.get_one::<String>("output").map(String::as_str),
        )?,
        _ => unreachable!(),
    }
    Ok(())
}
