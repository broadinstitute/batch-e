# Batch Effect Analysis Pipeline

A general-purpose pipeline for measuring group-level batch effects in whole-genome sequencing (WGS) variant callsets. Compares a user-supplied grouping variable (e.g., sequencing center, instrument, lab) across genomic interval classes, stratified by genetic ancestry, using Hail on Spark.

## Overview

Large-scale WGS consortia aggregate data from multiple sources, each with potentially different instruments, library prep protocols, and bioinformatics pipelines. These technical differences can introduce systematic batch effects that confound downstream analyses, particularly for rare variant studies.

This pipeline implements a **within-ancestry stratified analysis** to separate technical batch effects from population structure:

1. Stratify samples by genetic ancestry to avoid confounding
2. Compute per-sample variant statistics (SNPs, indels, Ti/Tv ratios) across genomic regions
3. Perform pairwise statistical comparisons between groups (Cohen's d, Welch's t-test)
4. Focus on genomic regions known to be technically challenging (low mappability, GC extremes, etc.)

## Pipeline Steps

| Step | Description |
|------|-------------|
| 1 | Load ancestry + comparison metadata TSVs, merge, subsample |
| 2 | Load variant data (VCF or pre-built MatrixTable), filter to target intervals |
| 3 | Annotate samples with ancestry and comparison group labels |
| 4 | Compute per-sample variant counts by interval (explode + group_by aggregation) |
| 5 | Extract results to pandas; compute derived metrics and pairwise comparisons |
| 6 | Save results and timing metrics |

## Usage

### In a notebook (Terra/Dataproc)

```python
from batch_e import PipelineConfig, run_pipeline

cfg = PipelineConfig(
    input_path='gs://bucket/data.mt/',
    ancestry_tsv='gs://bucket/ancestry.tsv',
    ancestry_col='ancestry_pred_other',
    comparison_tsv='gs://bucket/site_labels.tsv',
    comparison_col='site_id',
    comparison_name='sequencing_center',
    ancestries=['eur', 'afr'],
    interval_dict={
        'ACMG59': 'gs://bucket/interval_lists/acmg59.bed',
        'Low_Mappability': 'gs://bucket/interval_lists/lowmap.bed.gz',
    },
    samples_per_group=5000,
)

results = run_pipeline(cfg)

# Pairwise comparisons with effect sizes
comparisons = results['comparisons']

# Filter to appreciable effects
significant = comparisons[comparisons['cohens_d'].abs() > 0.5]
```

### From the command line

```bash
python batch_e.py \
    --input-path gs://bucket/data.mt/ \
    --ancestry-tsv gs://bucket/ancestry.tsv \
    --ancestry-col ancestry_pred_other \
    --comparison-tsv gs://bucket/site_labels.tsv \
    --comparison-col site_id \
    --comparison-name sequencing_center \
    --interval ACMG59=gs://bucket/acmg59.bed \
    --interval Low_Mappability=gs://bucket/lowmap.bed.gz \
    --samples-per-group 5000 \
    --output-dir gs://bucket/results/
```

### Generate a report

```bash
python batch_e_reporter.py gs://bucket/results/ -o report.html
```

### Via WDL on Terra (hailrunner)

The pipeline can be run as a WDL workflow on Terra, using [hailrunner](https://github.com/broadinstitute/hailrunner) to manage ephemeral Dataproc clusters. The workflow runs analysis on Dataproc and generates an HTML report in a lightweight follow-up task.

Import the workflow from [Dockstore](https://dockstore.org/) or directly from `wdl/batch_e.wdl`.

Intervals default to the 5 standard interval sets (ACMG59, Low_Mappability, GC_gt_85, GC_lt_25, HighConf_Genome) served from GitHub. HTTPS URLs are automatically staged to GCS before the Dataproc analysis runs. The data source is auto-detected from which path is provided (`vcf_path` or `mt_path`).

Example input JSON:

```json
{
    "batch_e.input_path": "gs://bucket/data.mt/",
    "batch_e.ancestry_tsv": "gs://bucket/ancestry.tsv",
    "batch_e.comparison_tsv": "gs://bucket/site_labels.tsv",
    "batch_e.comparison_col": "site_id",
    "batch_e.output_dir": "gs://staging-bucket/batch_effect_results/run_001",
    "batch_e.staging_bucket": "gs://staging-bucket",
    "batch_e.workers": 32
}
```

Override intervals if needed:

```json
"batch_e.intervals": ["MyRegion=gs://bucket/my_regions.bed.gz"]
```

## Configuration

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `ancestry_tsv` | Path to TSV with sample ancestry predictions |
| `comparison_tsv` | Path to TSV with comparison group labels |
| `comparison_col` | Column name for the grouping variable |
| `interval_dict` | Dict of interval name → BED file path (non-empty) |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `comparison_name` | `comparison_col` | Human-readable label for the comparison |
| `ancestry_col` | `ancestry_pred_other` | Column name for ancestry in ancestry TSV |
| `sample_id_col` | `research_id` | Sample ID column shared across all TSVs |
| `ancestries` | all | Ancestries to include |
| `comparison_values` | all | Comparison group values to include |
| `samples_per_group` | None (all) | Max samples per (group, ancestry) |
| `input_path` | — | VCF glob or MatrixTable path (auto-detected from `.mt` suffix) |
| `filter_to_pass` | True | Keep only PASS variants |
| `filter_FT_pass` | True | Set GT to missing unless FT contains PASS |
| `pruning_subsample_n` | 50,000 | Subsample large interval lists for partition pruning |

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

Designed to run on Google Cloud Dataproc via [Terra](https://terra.bio/) or via [hailrunner](https://github.com/broadinstitute/hailrunner).

## Project Structure

```
.
├── batch_e/
│   ├── batch_e.py              # Main pipeline
│   └── batch_e_reporter.py     # HTML report generator
├── wdl/
│   └── batch_e.wdl             # WDL workflow (hailrunner + reporter)
├── intervals/
│   ├── *.bed.gz                # Genomic interval files
│   ├── scripts/                # Interval download and processing utilities
│   └── downsample_report/      # Interval downsampling analysis
└── .dockstore.yml              # Dockstore workflow registration
```

## License

This project is licensed under the [MIT License](LICENSE).
