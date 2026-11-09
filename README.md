# mopdepth

[![build](https://github.com/kojix2/mopdepth/actions/workflows/build.yml/badge.svg)](https://github.com/kojix2/mopdepth/actions/workflows/build.yml)
[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Fkojix2%2Fmopdepth%2Flines)](https://tokei.kojix2.net/github/kojix2/mopdepth)
![Static Badge](https://img.shields.io/badge/PURE-VIBE_CODING-magenta)

An LLM-assisted port of [mosdepth](https://github.com/brentp/mosdepth) to Crystal
for BAM/CRAM depth calculation. LLMs were used to translate and adapt the original
implementation. It provides per-base depth, region statistics, coverage
distributions, quantization, and threshold counts.

## Installation

Building from source requires Crystal, Shards, HTSlib, and Make.

```sh
git clone https://github.com/kojix2/mopdepth
cd mopdepth
make
bin/mopdepth --help
```

`make` builds an optimized release binary at `bin/mopdepth`. Use `make debug`
for a development build. D4 output is included through
[d4.cr](https://github.com/bio-cr/d4.cr), with no Rust or native D4 dependency.

## Usage

Input must be an indexed BAM or CRAM file. For CRAM, supply the reference FASTA
with `-f` or `REF_PATH`.

```sh
bin/mopdepth [options] <prefix> <BAM-or-CRAM>

# Per-base depth and summary statistics
bin/mopdepth output sample.bam

# BED regions or 1 kb windows, without per-base output
bin/mopdepth -n -b regions.bed regions sample.bam
bin/mopdepth -n -b 1000 windows sample.bam

# CRAM input or D4 output
bin/mopdepth -f reference.fa cram sample.cram
bin/mopdepth --d4 depth sample.bam
```

Use `bin/mopdepth --help` for the full option list.

| Option | Purpose |
| --- | --- |
| `-t, --threads N` | BAM decompression threads |
| `-c, --chrom REGION` | Select a chromosome or a 1-based inclusive range, such as `chr1:100-200` |
| `-b, --by BED_OR_SIZE` | BED regions or a positive window size |
| `-f, --fasta FASTA` | CRAM reference; defaults to `REF_PATH` |
| `-n, --no-per-base` | Skip per-base output |
| `--d4` | Replace per-base BED.gz with D4 |
| `-Q, --mapq N` | Minimum mapping quality |
| `-l, --min-frag-len N` | Minimum fragment length |
| `-u, --max-frag-len N` | Maximum fragment length |
| `-x, --fast-mode` | Calculate depth from read start/end positions |
| `-a, --fragment-mode` | Count full fragments from proper pairs |
| `-m, --use-median` | Use median instead of mean for region statistics |
| `-q, --quantize BINS` | Quantized output, for example `0:1:4:` |
| `-T, --thresholds LIST` | Region coverage counts at comma-separated thresholds; requires `--by` |
| `-F, --flag MASK` | Exclude records with any specified FLAG bit; default `1796` |
| `-i, --include-flag MASK` | Include records with any specified FLAG bit |
| `-R, --read-groups LIST` | Select comma-separated read group IDs |
| `--exclude-tag TAG:VALUE` | Exclude exact string AUX tag matches; comma-separated and repeatable |
| `-M, --mos` | Use `mosdepth` instead of `mopdepth` in summary and distribution filenames |

The default mode uses CIGAR operations to calculate depth. Fast mode and fragment
mode are mutually exclusive. Fragment mode examines the complete chromosome even
with `--chrom`, so fragments spanning the selected range are included.

### Regions

BED input may be plain text, gzip, or BGZF. Coordinates must be zero-based,
half-open, non-empty, and within the named reference. Unknown references are
reported and ignored. Overlapping intervals are counted independently.

With `-n -b regions.bed`, chromosomes absent from the BED are skipped unless
`-q` or `-T` is also requested. Skipped chromosomes do not contribute to summaries,
distributions, or `total`. For processed chromosomes, `total` and the global
distribution cover the entire chromosome or selected `--chrom` range;
`total_region` and the region distribution cover the BED intervals.

### Tag filtering

To keep PCR/library duplicates while excluding optical duplicates:

```sh
bin/mopdepth -F 772 --exclude-tag DT:SQ output sample.bam
```

The default FLAG mask excludes all duplicate-flagged records, so changing `-F`
is necessary for this use. Tag conditions are combined with OR and apply to
individual records in every mode. Matching is exact and case-sensitive; missing
or non-string tags do not match. Tag names must match `[A-Za-z][A-Za-z0-9]`.
Values may contain colons, spaces, or be empty (`XX:`); quote shell arguments
containing spaces. Commas separate conditions.

## Output

`<label>` is `mopdepth` by default and `mosdepth` with `-M`.

| File | Written when |
| --- | --- |
| `<prefix>.<label>.summary.txt` | Always |
| `<prefix>.<label>.global.dist.txt` | Always |
| `<prefix>.per-base.bed.gz` | Per-base output is enabled without `--d4` |
| `<prefix>.per-base.d4` | `--d4` replaces per-base BED.gz; `-n` disables both |
| `<prefix>.regions.bed.gz` | `--by` is set |
| `<prefix>.<label>.region.dist.txt` | `--by` is set |
| `<prefix>.quantized.bed.gz` | `--quantize` is set |
| `<prefix>.thresholds.bed.gz` | `--thresholds` and `--by` are set |

BED outputs are BGZF-compressed and receive CSI indexes after successful
completion. D4 output receives an embedded secondary frame index.

The summary columns are `chrom`, `length`, `bases`, `mean`, `min`, and `max`.
`length` is the number of bases summarized; `bases` is the sum of their depths.

Histogram-based statistics have fixed limits:

| Statistic | Limit |
| --- | --- |
| Region median (`-m`) | Depths above 65,535 share bin 65,535 |
| Numeric-window region distribution | Rounded means or medians above 511 share bin 511 |
| Global and BED-region distributions | Depths above 400,000 share bin 399,990 |

These limits leave per-base depth, threshold counts, quantization, and summary
mean and maximum depth unchanged.

## Library API

```crystal
require "mopdepth"

options = Mopdepth::Options.new(minimum_mapping_quality: 20)
control = Mopdepth::Control.new { |progress| STDERR.puts progress }

Mopdepth::Reader.open("sample.bam", threads: 4) do |reader|
  reader.each_coverage(options: options, control: control) do |coverage|
    puts "#{coverage.region.chromosome}: #{coverage.values.size} bases"
  end
end
```

`Reader` provides `coverage`, `coverage_into`, and `each_coverage`. Coordinates
are zero-based and half-open. `coverage_into` reuses a caller-owned
`Array(Int32)`; results contain one depth per base. Use one reader per concurrent
worker. The block form closes the reader automatically.

`Options` supports `Mode::Standard`, `Mode::Fast`, and `Mode::Fragment`. Tag
exclusions use `exclude_tags: [Mopdepth::TagFilter.new("DT", "SQ")]`; FLAG
filtering is set separately through `exclude_flags`. Pass `control:` for progress
callbacks and cooperative cancellation; another fiber or thread can call
`control.cancel`. For CRAM, pass `reference: "reference.fa"` to `Reader.open`.

## mosdepth compatibility

Coverage, filtering, NoData, thresholds, quantization, and window distributions
follow mosdepth conventions. Compatibility checks compare decoded output and
indexed query results.

Mopdepth uses the true intersection when removing overlap from contained read
pairs, includes the final reference base in quantization, and applies `--chrom`
ranges to calculation and output. These intentionally differ from mosdepth
0.3.x edge cases. Region means use integer depth sums before formatting.

## Development

```sh
make test
make test-d4
make lint
make format-check
```

## License

[MIT](LICENSE)
