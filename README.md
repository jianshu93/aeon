



# Aeon: efficient and robust genome similarity search with evolving disk-based proximity graphs

***"Built once, updated forever. Nodes and edges evolve like cells, while Aeon remains eternal."***


Aeon is a dynamic genome similarity search platform built on mmap-backed DiskANN index, MERIT deletion repair, and b-bit one-permutation MinHash sketches. Search indexes remain in the ordinary static DiskANN format: updates use a temporary dynamic workspace and commit a compact static index.

## Install

```bash
git clone https://github.com/jianshu93/aeon
cd aeon
cargo build --release
./target/release/aeon --help
```

Input lists contain one genome path or stored genome name per line. FASTA/FASTQ and compressed inputs supported by Needletail are accepted. DNA and amino-acid databases are supported.

Run `aeon COMMAND --help` for the full description of any command and its options. `--threads` defaults to all logical CPUs for every command.

## Commands

### Build

```bash
aeon todisk \
  --prefix streptomyces \
  --reference-list genomes.txt \
  --kmer-size 16 \
  --sketch-size 12000 \
  --max-degree 64 \
  --build-beam-width 256
```

This creates:

- `streptomyces.diskann`: static mmap-searchable DiskANN index.
- `streptomyces.genomes.txt`: genome names in current vector-ID order.
- `streptomyces.idmap.tsv`: vector ID to genome name mapping.
- `streptomyces.params.json`: sketch and graph parameters.

Important build options:

| Option | Default | Meaning |
|---|---:|---|
| `--seq-type <dna\|aa>` | `dna` | Select the DNA or amino-acid sketcher. |
| `--kmer-size <K>` | `16` | k-mer length. DNA supports up to 32 except 15; amino acid supports up to 12. |
| `--sketch-size <N>` | `8192` | Number of `u16` values per b-bit MinHash sketch. Larger sketches improve resolution but increase index size and distance cost. |
| `--densification <0\|1>` | `0` | `0` selects OPTDENS; `1` selects RevOPTDENS. |
| `--max-degree <R>` | `64` | Maximum Vamana out-degree. Larger values usually improve recall while increasing storage and traversal work. |
| `--build-beam-width <L>` | `256` | Construction beam. Larger values generally improve graph quality at higher build cost. |
| `--alpha <FLOAT>` | `1.2` | Vamana robust-pruning alpha. |
| `--extra-seeds <N>` | `1` | Additional random construction-search entry points per node. |
| `--hash-seed <INTEGER>` | `1337` | XXH3 seed saved in metadata and reused by insert and search. |

### Insert

The insertion list must contain paths to genome files. Aeon sketches them with the database's original parameters.

```bash
aeon insert --prefix streptomyces --genome-list insert.txt
```

`--beam-width` controls candidate exploration while connecting inserted genomes. It defaults to the original build beam stored in `PREFIX.params.json`; increasing it may improve graph quality at the cost of update time.

### Delete

Deletion does not read or sketch genomes. Every line must exactly match a name stored in `PREFIX.genomes.txt`.

```bash
aeon delete --prefix streptomyces --name-list delete.txt
```

Aeon uses MERIT versioned-edge invalidation and local repair, then removes dynamic version state while committing the updated ordinary static index.

`--name-list` must contain the complete stored strings from `PREFIX.genomes.txt`; matching is exact and basename-only matching is not performed.

### Combined Update

Delete and insert in one update session and one static commit:

```bash
aeon update \
  --prefix streptomyces \
  --delete-list delete.txt \
  --insert-list insert.txt
```

A path may be deleted and reinserted in the same command. Names that remain in the database cannot be inserted again.

`--beam-width` has the same meaning as for `insert`. Deletions use the rust-diskann MERIT defaults (`repair beam = 2R`, `k_r = 2`).

### Search

```bash
aeon search \
  --prefix streptomyces \
  --query-list queries.txt \
  --k 10 \
  --beam-width 512 \
  --output neighbors.tsv
```

Search always opens the static index. Results contain query path, matched genome name, current vector ID, raw sketch Jaccard estimate, and normalized Hamming distance (or hash collision probability).

Search options:

| Option | Default | Meaning |
|---|---:|---|
| `--k <N>` | `10` | Nearest neighbors reported for each query. |
| `--beam-width <N>` | `512` | Graph-search beam. Larger values generally improve recall but reduce throughput. |
| `--output <FILE>` | stdout | Write tab-separated results to a file. |

## Update Semantics

- Updates are written to temporary files before replacing database files.
- Vector IDs may change after deletion because live records are compacted.
- Aeon rebuilds `.genomes.txt` and `.idmap.tsv` from rust-diskann's old-to-new ID mapping after every update.
- `insert`, `delete`, and `update` finish with the same static format produced by `todisk`.

## Benchmark

Apple Silicon M4 Max benchmark on 2,799 RefSeq Streptomycetaceae genomes (6.4 GB input), using 8 threads, DNA k=16, OPTDENS sketch size 8,192, graph degree 64, and construction beam/ef 256:

| Operation | aeon | GSearch | Comparison |
|---|---:|---:|---:|
| Build wall time | 50.38 s | 42.22 s | aeon 19% slower |
| Build peak RSS | 675 MiB | 991 MiB | aeon 32% lower |
| Database size, including mappings | 48.23 MB | 47.43 MB | aeon 1.7% larger |
| Search 100 genomes | 1.63 s | 2.40 s | aeon 1.47x faster |
| Search peak RSS | 394 MiB | 571 MiB | aeon 31% lower |
| Top-1 zero-distance matches | 100/100 | 100/100 | equal |

A combined 1% replacement (28 MERIT deletions followed by 28 insertions and one static commit) took 4.53 seconds with 222 MiB peak RSS. This was 11.1x faster than rebuilding the aeon database, and all 28 replacement queries returned a zero-distance top hit.

## References

- Subramanya et al., 2019 *DiskANN: Fast Accurate Billion-point Nearest Neighbor Search on a Single Node*. NeurIPS 2019.
- Wu et.al., 2026. *MERIT: Efficient In-Place Deletion for Dynamic Graph-Based Approximate Nearest Neighbor Indexes*. arXiv preprint arXiv:2607.29173.
- Li et.al., 2012. *One permutation hashing*. Advances in Neural Information Processing Systems 25.
- Li et al., 2010. *b-Bit Minwise Hashing*. WWW 2010.
