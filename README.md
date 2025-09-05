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
- Window-based analysis
- Comprehensive filtering options (MAPQ, fragment length, flags)

## Installation

### Prerequisites

- Crystal
- hts-lib (for BAM/CRAM support)

### Build from source

```bash
git clone https://github.com/kojix2/mopdepth
cd mopdepth
shards install
shares build --release
```

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
- `-n, --no-per-base`: Skip per-base output
- `-Q, --mapq MAPQ`: MAPQ threshold
- `-l, --min-frag-len MIN`: Minimum fragment length
- `-u, --max-frag-len MAX`: Maximum fragment length
- `-x, --fast-mode`: Fast mode (read start/end positions only)
- `-a, --fragment-mode`: Count full fragment (proper pairs only)
- `-m, --use-median`: Use median for region stats instead of mean
- `-M, --mos`: Use mosdepth-compatible filenames (mosdepth.*); default is depth.*
- `-v, --version`: Show version
- `-h, --help`: Show help message

### Processing modes

- **Default mode**: CIGAR-based depth calculation (most accurate)
- **Fast mode** (`-x`): Uses read start/end positions (faster but less accurate)
- **Fragment mode** (`-a`): Counts full fragments for paired-end reads

**Note**: Fast mode and fragment mode cannot be used together.

### Output files

- Summary: `<prefix>.(mopdepth|mosdepth).summary.txt`
- Per-base: `<prefix>.per-base.bed` (unless `-n`)
- Global dist: `<prefix>.(mopdepth|mosdepth).global.dist.txt`
- Regions: `<prefix>.regions.bed` (when `--by`)
- Region dist: `<prefix>.(mopdepth|mosdepth).region.dist.txt` (when `--by`)

By default, files are named with the `mopdepth.*` label. Use `-M/--mos` to switch to `mosdepth.*`.

### Summary file format

The summary file contains the following columns:

- `chrom`: Chromosome name
- `length`: Chromosome length
- `sum_depth`: Total depth (sum of all depths)
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
