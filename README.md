# mopdepth

[![build](https://github.com/kojix2/mopdepth/actions/workflows/build.yml/badge.svg)](https://github.com/kojix2/mopdepth/actions/workflows/build.yml)
[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Fkojix2%2Fmopdepth%2Flines)](https://tokei.kojix2.net/github/kojix2/mopdepth)
![Static Badge](https://img.shields.io/badge/PURE-VIBE_CODING-magenta)

A fast BAM/CRAM depth calculation tool written in Crystal, inspired by [mosdepth](https://github.com/brentp/mosdepth).

**This is an experiment to see if well-known tools can be ported to Crystal using “vibe coding”.**

## Features

- Fast depth calculation for BAM/CRAM files
- Multiple processing modes (fast mode, fragment mode, CIGAR-based)
- Per-base and region-based depth analysis
- BED file support for custom regions
- Plain, gzip, and BGZF BED input
- Explicit FASTA references for CRAM input
- Optional D4 per-base output
- Window-based analysis
- Comprehensive filtering options (MAPQ, fragment length, flags)

## Installation

### Prerequisites

- Crystal
- hts-lib (for BAM/CRAM support)
- Make

### Build from source

```bash
git clone https://github.com/kojix2/mopdepth
cd mopdepth
make
```

The default target always builds an optimized release binary at `bin/mopdepth`.
Use `make debug` only when a non-optimized development binary is needed.

To build with D4 support, install Rust/Cargo and run:

```bash
make d4
```

This downloads the pinned D4 source, builds `libd4binding.a`, and links that
archive into the optimized mopdepth binary. The D4 binding is therefore not a
runtime shared-library dependency. Other platform libraries, including HTSlib,
follow the normal Crystal and system linker configuration. Run `make test-d4`
to execute the D4 write/read-back integration tests.

## Usage

```bash
./mopdepth [options] <prefix> <BAM-or-CRAM>
```

### Basic example

```bash
./mopdepth output sample.bam
```

### Options

- `-t, --threads THREADS`: BAM decompression threads
- `-c, --chrom CHROM`: Restrict to chromosome
- `-b, --by BY`: BED file or numeric window size
- `-f, --fasta FASTA`: FASTA reference for CRAM input (defaults to `REF_PATH`)
- `-n, --no-per-base`: Skip per-base output
- `--d4`: Write per-base depth as D4 (requires a `make d4` build)
- `-Q, --mapq MAPQ`: MAPQ threshold
- `-l, --min-frag-len MIN`: Minimum fragment length
- `-u, --max-frag-len MAX`: Maximum fragment length
- `-x, --fast-mode`: Fast mode (read start/end positions only)
- `-a, --fragment-mode`: Count full fragment (proper pairs only)
- `-m, --use-median`: Use median for region stats instead of mean
- `-q, --quantize QUANTIZE`: Write quantized output (for example, `0:1:4:`)
- `-T, --thresholds THRESHOLDS`: Comma-separated thresholds for region coverage
- `-F, --flag FLAG`: Exclude reads with FLAG bits set
- `-i, --include-flag FLAG`: Include only reads with FLAG bits set
- `-R, --read-groups GROUPS`: Comma-separated read group IDs
- `-M, --mos`: Use mosdepth-compatible filenames (mosdepth.*); default is mopdepth.*
- `-v, --version`: Show version
- `-h, --help`: Show help message

### Processing modes

- **Default mode**: CIGAR-based depth calculation (most accurate)
- **Fast mode** (`-x`): Uses read start/end positions (faster but less accurate)
- **Fragment mode** (`-a`): Counts full fragments for paired-end reads

**Note**: Fast mode and fragment mode cannot be used together.

Window sizes must be positive. BED intervals are validated against the alignment
header: coordinates must be 0-based half-open, non-empty, non-negative, and
contained in the named reference. References absent from the alignment header are
reported and ignored. Overlapping BED intervals remain independent observations.

`--chrom` also accepts a 1-based inclusive range such as `chr1:100-200`. Unlike
mosdepth 0.3.x, mopdepth applies that range to calculation and output. In fragment
mode it still examines the complete chromosome so fragments spanning the selected
range are not missed.

### Output files

- Summary: `<prefix>.(mopdepth|mosdepth).summary.txt`
- Per-base: `<prefix>.per-base.bed.gz` (unless `-n`)
- D4 per-base: `<prefix>.per-base.d4` (`--d4`, replacing per-base BED.gz)
- Global dist: `<prefix>.(mopdepth|mosdepth).global.dist.txt`
- Regions: `<prefix>.regions.bed.gz` (when `--by`)
- Region dist: `<prefix>.(mopdepth|mosdepth).region.dist.txt` (when `--by`)
- Quantized: `<prefix>.quantized.bed.gz` (when `--quantize`)
- Thresholds: `<prefix>.thresholds.bed.gz` (when `--thresholds` and `--by`)

By default, files are named with the `mopdepth.*` label. Use `-M/--mos` to switch to `mosdepth.*`.

### mosdepth compatibility

mopdepth follows mosdepth's coverage, filtering, NoData, threshold, quantization,
and window-distribution conventions where those conventions are well-defined.
Compressed output and CSI files are compared by decoded content and query results,
not by their binary bytes.

Two known mosdepth edge-case errors are intentionally not reproduced: contained
read pairs use the true intersection when removing mate overlap, and quantization
examines the final real base of every reference. Region means are accumulated as
integer depth sums before formatting to avoid floating-point accumulation drift.

D4-enabled builds support mosdepth-style `--d4` per-base output and create an
embedded secondary frame index after successfully closing the file. `-M/--mos`
changes compatible output names; it does not enable bug-for-bug compatibility.

### Summary file format

The summary file contains the following columns:

- `chrom`: Chromosome name
- `length`: Chromosome length
- `bases`: Total depth (sum of all depths)
- `mean`: Mean depth
- `min`: Minimum depth
- `max`: Maximum depth

## Examples

### Basic depth calculation

```bash
./mopdepth output sample.bam
```

### With BED regions

```bash
./mopdepth -b regions.bed output sample.bam
```

### Window-based analysis (1kb windows)

```bash
./mopdepth -b 1000 output sample.bam
```

### Fast mode with MAPQ filtering

```bash
./mopdepth -x -Q 20 output sample.bam
```

### Fragment mode for paired-end data

```bash
./mopdepth -a -l 100 -u 1000 output sample.bam
```

## License

MIT License
