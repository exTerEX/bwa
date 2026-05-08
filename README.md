# BWA-Style FM-Index Aligner (C++20)

A modern C++20 implementation of core Burrows-Wheeler Aligner building blocks:

- SA-IS suffix array construction
- Burrows-Wheeler Transform (BWT)
- FM-index with rank/select acceleration
- Exact and bounded inexact matching (mismatches + simple indels)
- Seed-and-extend inexact mode for faster candidate generation
- Compressed suffix array sampling with LF-based locate recovery

This project is structured similarly to the sibling repositories `bruijn` and `water`:

- CMake-based C++20 build
- `include/`, `src/`, `tests/`, `data/`, `scripts/`
- MIT license and runnable example test script

## Features

- Suffix Array (SA-IS)
  - Induced sorting over integer alphabets
  - Sentinel-terminated text (`$`) for BWT/FM-index construction

- BWT + FM-index
  - Backward search over compressed transform
  - Per-symbol succinct bit-vectors
  - Fast `rank1` and `select1` via superblocks + word-level popcount
  - Optional compressed SA representation (`--sa-sample-rate`)
    - Stores sampled SA rows/values only
    - Recovers unsampled suffix positions via LF-mapping walks

- Inexact matching
  - Two strategies:
    - bounded DFS over FM-index intervals (`backtracking`)
    - seed-and-extend candidate refinement (`seed`)
  - `auto` mode selects seed-and-extend on longer reads and falls back to backtracking
  - Supports:
    - substitutions (mismatches)
    - deletions in query
    - insertions in query
  - Tunable `--max-edits`, `--search-mode`, and `--seed-candidates`

## Project Layout

- `include/bwa/dna.hpp`
  - DNA symbol encode/decode and normalization
- `include/bwa/sais.hpp`
  - SA-IS suffix array implementation
- `include/bwa/bit_vector.hpp`
  - Rank/select bit-vector structure
- `include/bwa/fm_index.hpp`
  - BWT/FM-index and exact/inexact search
- `include/bwa/fasta.hpp`, `src/fasta.cpp`
  - FASTA/FASTQ parser
- `include/bwa/index_io.hpp`, `src/index_io.cpp`
  - Binary index serialization
- `src/main.cpp`
  - CLI entry point
- `tests/test_fm_index.cpp`
  - Unit-style correctness tests

## Build

With CMake:

```bash
cmake -S . -B build
cmake --build build -j
```

Or via Make wrapper:

```bash
make
```

## CLI

### Build Index

```bash
./build/bwa index --reference data/example_reference.fasta --output build/example.bwaidx
```

Build a sampled index (smaller memory footprint at the cost of slower `locate`):

```bash
./build/bwa index --reference data/example_reference.fasta --output build/example_sampled.bwaidx --sa-sample-rate 8
```

### Align Reads (Exact)

```bash
./build/bwa align --index build/example.bwaidx --reads data/example_reads.fasta --max-edits 0
```

### Align Reads (Inexact)

```bash
./build/bwa align --index build/example.bwaidx --reads data/example_reads.fasta --max-edits 1 --max-hits 10
```

Force seed-and-extend mode:

```bash
./build/bwa align --index build/example.bwaidx --reads data/example_reads.fasta --max-edits 2 --search-mode seed --seed-candidates 8192
```

## Output Format

The align command prints one metadata line and one line per read:

```text
command=align reads=3 reference_bases=210 max_edits=1 max_hits=10 search_mode=auto seed_candidates=4096
read=read_exact len=32 hits=3 best_edits=0 matches=48:0,47:1,49:1
```

For exact mode (`--max-edits 0`), output uses `positions=`; for inexact mode it uses `matches=position:edits`.

## Tests

Run CTest suite:

```bash
ctest --test-dir build --output-on-failure
```

Run example integration script:

```bash
./scripts/test_example.sh
```

## Benchmarking

Compare inexact strategies (`backtracking` vs `seed`) and full-SA vs sampled-SA locate latency:

```bash
./scripts/benchmark_compare.sh \
  --reference data/example_reference.fasta \
  --reads data/example_reads.fasta \
  --repeat 20 \
  --max-edits 1 \
  --sa-sample-rate 8
```

Write CSV output for downstream plotting or analysis:

```bash
./scripts/benchmark_compare.sh \
  --reference data/example_reference.fasta \
  --reads data/example_reads.fasta \
  --repeat 20 \
  --max-edits 1 \
  --sa-sample-rate 8 \
  --csv-out build/benchmark.csv
```

Append to an existing CSV file:

```bash
./scripts/benchmark_compare.sh \
  --reference data/example_reference.fasta \
  --reads data/example_reads.fasta \
  --csv-out build/benchmark.csv \
  --csv-append
```

The script prints two sections:

- Inexact search comparison: full/sampled index crossed with backtracking/seed mode
- Exact locate latency comparison: full SA vs sampled SA

Typical output fields:

- `total_ms`: total wall-clock time across all repeats
- `avg_ms`: average wall-clock time per run
- `total_hits`: total reported hits across all reads (from the final run)

CSV columns:

- `section,label,index_kind,search_mode,max_edits,repeat,max_hits,sa_sample_rate,seed_candidates,total_ms,avg_ms,total_hits,reference,reads`

Use your real datasets by passing your own `--reference` and `--reads` paths.

## Notes and Limitations

- This implementation is focused on algorithmic clarity and correctness in C++20.
- Inexact matching uses bounded backtracking and is best suited for small edit budgets.
- Inexact matching now supports both backtracking and seed-and-extend; `auto` mode chooses based on query length.
- SA sampling reduces index memory usage and index file size but increases locate cost because positions are recovered by LF walks.
- The current suffix array path uses 32-bit SA-IS internals, suitable for medium/large references but not multi-billion base-scale production indexing yet.
