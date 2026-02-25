# wgs-batch-effect-analysis

Measuring and characterizing sequencing center batch effects in All of Us whole-genome sequencing data. Computes per-sample variant statistics across genomic interval classes and performs pairwise comparisons between sequencing centers, stratified by genetic ancestry, using Hail on Spark.

## Overview

Large-scale WGS consortia aggregate data from multiple sequencing centers, each with potentially different instruments, library prep protocols, and bioinformatics pipelines. These technical differences can introduce systematic batch effects that confound downstream analyses, particularly for rare variant studies.

This pipeline implements a **within-ancestry stratified analysis** to separate technical batch effects from population structure:

1. Stratify samples by genetic ancestry to avoid confounding
2. Compute per-sample variant statistics (SNPs, indels, Ti/Tv ratios) across genomic regions
3. Perform pairwise statistical comparisons between sequencing centers (Cohen's d, Welch's t-test)
4. Focus on genomic regions known to be technically challenging (low mappability, GC extremes, etc.)

## Pipeline Steps

| Step | Description |
|------|-------------|
| 1 | Load sample metadata (ancestry, batch labels) and subsample |
| 2 | Load variant data (VCF or pre-built MatrixTable), filter to target intervals |
| 3 | Annotate samples with batch and ancestry labels |
| 4 | Compute per-sample variant counts by interval (explode + group_by aggregation) |
| 5 | Extract results to pandas; compute derived metrics and pairwise comparisons |
| 6 | Save results and timing metrics |

## Run Modes

| Mode | Samples/Group | Chromosomes | Use Case |
|------|---------------|-------------|----------|
| `smoke` | 100 | chr17 | Quick validation (~10 min) |
| `dev` | 1,000 | chr7, chr13, chr17 | Development (~2-3 hrs) |
| `medium` | 5,000 | All | Statistical validation |
| `full` | 10,000 | All | Production analysis |

## Genomic Intervals

| Interval | Description | Count |
|----------|-------------|------:|
| ACMG59 | Clinically actionable genes | 1,221 |
| Low_Mappability | Regions prone to alignment artifacts | 673,815 |
| GC_gt_85 | High GC content (>85%) | 1,998 |
| GC_lt_25 | Low GC content (<25%) | 487,211 |
| HighConf_Genome | GIAB high-confidence calling regions | 1,047,093 |

## Usage

### In a notebook (Terra/Dataproc)

```python
from batch_e import PipelineConfig, run_pipeline

cfg = PipelineConfig(
    mode='dev',
    ancestries=['eur', 'afr'],
    batches=['bi_S4', 'bcm_S4', 'uw_S4'],
)

results = run_pipeline(cfg)

# Pairwise comparisons with effect sizes
comparisons = results['comparisons']

# Filter to appreciable effects
significant = comparisons[comparisons['cohens_d'].abs() > 0.5]
```

### From the command line

```bash
python batch_e.py --mode dev \
    --batches bi_S4 bcm_S4 uw_S4 \
    --ancestries eur afr \
    --data-source vcf \
    --output-dir gs://my-bucket/results/
```

### Generate a report

```bash
python batch_e_reporter.py gs://my-bucket/results/ -o report.html
```

## Metrics

**Per-sample counts** (by interval and zygosity):
- SNP transitions / transversions
- Insertions / deletions

**Derived metrics:**
- Ti/Tv ratio
- Del/Ins ratio
- Total SNP and indel counts

**Comparison statistics:**
- Cohen's d (pooled standard deviation)
- Welch's t-test p-values

## Requirements

- Python 3.10+
- [Hail](https://hail.is/) (0.2.x)
- numpy, pandas, scipy
- matplotlib, seaborn (for reporting)
- gcsfs (for GCS access)

Designed to run on Google Cloud Dataproc via [Terra](https://terra.bio/).

## Project Structure

```
.
├── batch_e.py              # Main pipeline
├── batch_e_reporter.py     # HTML report generator
└── intervals/
    ├── *.bed               # Genomic interval files
    ├── scripts/            # Interval download and processing utilities
    └── downsample_report/  # Interval downsampling analysis
```
