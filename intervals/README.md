# Interval Lists for Batch Effect Analysis

## Motivation

The batch-effect pipeline measures per-sample variant statistics across genomic
interval classes. A ~10,000x spread in total base-pair coverage between
interval sets (ACMG59 at ~2 Mbp vs HighConf_Genome at ~2.7 Gbp) creates two
problems:

1. **Hail partition pruning fails** — `hl.filter_intervals` on a set covering
   ~90% of the genome cannot skip any partitions, so the "optimization" loads
   the entire MatrixTable anyway.
2. **Class imbalance in batch statistics** — Cohen's d from millions of
   variants in one class vs hundreds in another makes cross-interval
   comparisons misleading.

Downsampling the largest interval set to match the second-largest restores
effective partition pruning and brings coverage into the same order of
magnitude.

## Interval Inventory

| Interval | File | Count | Total bp | % GRCh38 | Source |
|----------|------|------:|----------|----------|--------|
| ACMG59 | `acmg59_allofus_19dec2019.GRC38.wGenes.NEW.bed` | 1,221 | ~2.3 Mbp | 0.07% | AoU clinical gene list |
| GC_gt_85 | `GRCh38_gc85_slop50.bed` | 1,998 | ~0.4 Mbp | 0.01% | GC content >85% ± 50bp slop |
| GC_lt_25 | `GRCh38_gclt25_merged.bed` | 487,211 | ~258 Mbp | 8.35% | GC content <25%, merged |
| Low_Mappability | `GRCh38_lowmappabilityall.bed` | 673,815 | ~284 Mbp | 9.20% | Low mappability regions |
| HighConf_Genome | `giab_highconf_wgs_calling_regions_hg38_intersection.bed` | 1,047,093 | ~2.7 Gbp | 87.52% | GIAB high-confidence calling regions |

## Downsampling Strategy

`downsample_bed.py sample` uses **chromosome-stratified, bp-proportional
sampling**:

1. Compute each chromosome's share of total bp in the input file.
2. Allocate the target bp budget proportionally (minimum 1 interval per
   chromosome).
3. Within each chromosome, randomly shuffle intervals and greedily accumulate
   until the budget is reached.

**Why base-pairs, not interval count?** Interval widths span 4 orders of
magnitude. Sampling by count would over-represent narrow intervals and
under-represent wide ones, distorting the coverage profile.

**Reproducibility:** All runs use `--seed 42` by default. The manifest JSON
records the exact seed, command, and per-chromosome breakdown.

## Usage

### Analyze existing interval files

```bash
# Auto-discover all .bed files in this directory
python downsample_bed.py probe --interval-dir intervals/

# Explicit NAME=path pairs
python downsample_bed.py probe \
    --beds ACMG59=intervals/acmg59.bed Low_Map=intervals/lowmap.bed
```

Outputs:
- Summary table to stdout
- `intervals/probe_report.md` — Markdown report with embedded figures

### Downsample a BED file

```bash
# Match another file's total bp coverage
python downsample_bed.py sample \
    --input intervals/giab_highconf_wgs_calling_regions_hg38_intersection.bed \
    --match-file intervals/GRCh38_lowmappabilityall.bed \
    --output intervals/giab_highconf_downsampled.bed

# Explicit target
python downsample_bed.py sample \
    --input big.bed --target-bp 250000000 --output out.bed

# Skip the report
python downsample_bed.py sample \
    --input big.bed --match-file ref.bed --output out.bed --no-report
```

Outputs alongside `--output`:
- `*.bed` — downsampled BED file
- `*.bed.manifest.json` — machine-readable provenance
- `*.bed.report.md` — before/after Markdown report with embedded figures

### Convert reports to PDF

```bash
pandoc intervals/probe_report.md -o probe_report.pdf
pandoc intervals/giab_highconf_downsampled.bed.report.md -o sample_report.pdf
```

## Dependencies

- Python 3.10+
- pandas
- numpy
- matplotlib
- seaborn
