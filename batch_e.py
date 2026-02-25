#!/usr/bin/env python3
"""
Batch Effect Analysis Pipeline for AoU WGS Data

A modular pipeline for measuring sequencing center batch effects in VCF data.
Supports multiple run modes (smoke, dev, medium, full) and both within-ancestry
and pooled analysis strategies.

Usage in notebook:
    # Paste this code into a cell, then:
    cfg = PipelineConfig(
        mode='dev',
        ancestries=['eur', 'afr'],
        batches=['bcm_v8_delta', 'bi_v8_delta'],
    )
    results = run_pipeline(cfg)

Usage from command line:
    python batch_effect_pipeline.py --mode dev --output-dir gs://bucket/results/
"""

from __future__ import annotations

import os
import json
import logging
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from enum import Enum

import numpy as np
import pandas as pd
import gcsfs

# Hail import is deferred to avoid startup cost when not needed
hl = None

def _init_hail():
    """Initialize Hail with proper Spark configuration."""
    import hail as hl
    
    # Configure Spark for large broadcast variables
    spark_conf = {
        'spark.rpc.message.maxSize': '1024',  # 1GB instead of 512MB
        'spark.driver.maxResultSize': '8g',
    }
    
    hl.init(
        default_reference='GRCh38', 
        idempotent=True,
        spark_conf=spark_conf
    )
    return hl

# =============================================================================
# Configuration
# =============================================================================

class RunMode(Enum):
    """Pipeline run modes with different scale/speed tradeoffs."""
    SMOKE = "smoke"      # ~100 samples, 1 interval - catch syntax errors
    DEV = "dev"          # ~1000 samples, 2 intervals - catch logic errors
    MEDIUM = "medium"    # ~5000 samples, all intervals - catch scaling issues
    FULL = "full"        # All samples, all intervals - production


@dataclass
class PipelineConfig:
    """Configuration for batch effect analysis pipeline."""
    
    # Run identification
    run_name: Optional[str] = None
    mode: str = "dev"  # smoke, dev, medium, full - controls variant filtering only
    
    # Sample subsetting
    samples_per_group: Optional[int] = None  # None = all samples (no subsampling)
        
    # Data source
    data_source: str = "vcf"          # "vcf" (primary) or "mt" (fast path using pre-built MT)
    vcf_path: Optional[str] = None    # glob pattern, e.g., "gs://.../*.vcf.bgz"
    hail_mt_path: Optional[str] = None  # pre-built MT path (used when data_source="mt")

    # VCF preprocessing
    split_multi: bool = True           # split multi-allelic sites (required for VCF input)
    acaf_filter: bool = False          # apply ACAF-style filtering (AF>threshold or AC>threshold)
    acaf_min_af: float = 0.01
    acaf_min_ac: int = 100

    # Caching (VCF import → MT for reuse)
    cache_mt: bool = True              # write imported VCF as MT for reuse
    cache_mt_path: Optional[str] = None  # auto-generated if None
    force_reimport: bool = False       # ignore cache, re-import from VCF

    # Other data paths (defaults from environment)
    workspace_bucket: Optional[str] = None
    metadata_path: Optional[str] = None
    
    # Batch and stratification
    batch_col: str = "batch"
    ancestry_col: str = "ancestry_pred_other"
    sample_id_col: str = "research_id"
    
    # Which batches and ancestries to include
    batches: Optional[List[str]] = None  # None = all
    ancestries: Optional[List[str]] = None  # None = all
    
    # Analysis parameters
    analysis_strategy: str = "within_ancestry"  # "within_ancestry" or "pooled"
    min_samples_per_group: int = 50  # Skip groups with fewer samples
    
    # Interval lists for stratified analysis
    interval_dict: Optional[Dict[str, str]] = None  # Auto-populated if None
    pruning_subsample_n: Optional[int] = 50_000  # Subsample large interval lists for partition pruning (None = no subsampling). Annotation always uses full intervals.

    # Sampling parameters (for non-full modes)
    samples_per_group: Optional[int] = None  # Auto-set based on mode
    random_seed: int = 42
    
    # Variant filtering
    filter_to_pass: bool = True  # Keep only PASS variants
    filter_FT_pass: bool = True  # Set GT to missing unless FT contains PASS
    
    # Output control
    output_dir: Optional[str] = None
    export_per_sample: bool = False  # Large output, off by default
    export_sampled_per_sample: bool = True
    sampled_per_sample_n: int = 2000
    
    # Compute parameters
    hail_min_partitions: int = 1000
    
    def __post_init__(self):
        """Set defaults based on mode and environment."""
        # Resolve environment variables
        if self.workspace_bucket is None:
            self.workspace_bucket = os.environ.get('WORKSPACE_BUCKET', '')
        
        if self.metadata_path is None:
            self.metadata_path = f"{self.workspace_bucket}/batchE_tables/sample_metadata.tsv"
        
        if self.vcf_path is None:
            self.vcf_path = os.environ.get('WGS_ACAF_THRESHOLD_VCF_PATH', '')

        if self.hail_mt_path is None:
            self.hail_mt_path = os.environ.get('WGS_ACAF_THRESHOLD_SPLIT_HAIL_PATH', '')

        if self.cache_mt_path is None:
            self.cache_mt_path = f"{self.workspace_bucket}/batch_effect_cache/imported.mt"

        # Set default samples_per_group based on mode if not explicitly provided
        if self.samples_per_group is None:
            mode_defaults = {
                "smoke": 100,
                "dev": 1000,
                "medium": 5000,
                "full": 10000,  # Even "full" mode defaults to 10k per group
            }
            self.samples_per_group = mode_defaults.get(self.mode, 5000)
    
        # Generate run name
        if self.run_name is None:
            timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
            self.run_name = f"batche_{self.mode}_{timestamp}"
        
        # Set output directory
        if self.output_dir is None:
            self.output_dir = f"{self.workspace_bucket}/batch_effect_results/{self.run_name}"
        
        # Default interval dict
        if self.interval_dict is None:
            self.interval_dict = self._default_interval_dict()
        
        # Adjust intervals for smaller modes
        if self.mode == "smoke":
            # Just use one small interval
            self.interval_dict = {"ACMG59": self.interval_dict.get("ACMG59", list(self.interval_dict.values())[0])}
        elif self.mode == "dev":
            # Use two intervals
            keys = list(self.interval_dict.keys())[:2]
            self.interval_dict = {k: self.interval_dict[k] for k in keys}
    
    def _default_interval_dict(self) -> Dict[str, str]:
        """Default genomic interval lists for analysis."""
        bucket = os.environ.get('WORKSPACE_BUCKET', '')
        return {
            "ACMG59": f"{bucket}/interval_lists/acmg59_allofus_19dec2019.GRC38.wGenes.NEW.bed",
            "Low_Mappability": f"{bucket}/interval_lists/GRCh38_lowmappabilityall.bed.gz",
            "GC_gt_85": f"{bucket}/interval_lists/GRCh38_gc85_slop50.bed.gz",
            "GC_lt_25": f"{bucket}/interval_lists/GRCh38_gclt25_merged.bed",
            "HighConf_Genome": f"{bucket}/interval_lists/giab_highconf_wgs_calling_regions_hg38_intersection.bed",
        }


# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(run_name: str, level: int = logging.INFO) -> logging.Logger:
    """Configure logging for the pipeline."""
    logger = logging.getLogger(f"batche.{run_name}")
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(level)
    return logger

def log_cluster_info(logger: logging.Logger) -> Dict[str, Any]:
    """Log Spark cluster configuration and size."""
    try:
        from pyspark import SparkContext
        sc = SparkContext.getOrCreate()
        
        conf = sc.getConf()
        
        # Get executor info
        executor_memory = conf.get("spark.executor.memory", "unknown")
        executor_cores = conf.get("spark.executor.cores", "unknown")
        driver_memory = conf.get("spark.driver.memory", "unknown")
        
        # Get number of executors (this is the actual count of active executors)
        # _jsc.sc().getExecutorMemoryStatus() returns a map of executor -> memory
        executor_status = sc._jsc.sc().getExecutorMemoryStatus()
        num_executors = executor_status.size() - 1  # subtract 1 for driver
        
        # Default parallelism
        parallelism = sc.defaultParallelism
        
        cluster_info = {
            'num_executors': num_executors,
            'executor_memory': executor_memory,
            'executor_cores': executor_cores,
            'driver_memory': driver_memory,
            'default_parallelism': parallelism,
            'app_name': sc.appName,
            'master': sc.master,
        }
        
        logger.info("=" * 60)
        logger.info("SPARK CLUSTER INFO")
        logger.info("=" * 60)
        logger.info(f"  Executors (workers): {num_executors}")
        logger.info(f"  Executor memory:     {executor_memory}")
        logger.info(f"  Executor cores:      {executor_cores}")
        logger.info(f"  Driver memory:       {driver_memory}")
        logger.info(f"  Default parallelism: {parallelism}")
        logger.info(f"  Total cores:         ~{num_executors * int(executor_cores) if executor_cores != 'unknown' else 'unknown'}")
        logger.info("=" * 60)
        
        return cluster_info
        
    except Exception as e:
        logger.warning(f"Could not get cluster info: {e}")
        return {'error': str(e)}


def collect_spark_metrics(logger: logging.Logger) -> Dict[str, Any]:
    """Collect Spark metrics after computation completes.

    Uses the Java SparkContext API directly since the PySpark StatusTracker
    has limited methods in Spark 3.5.x.
    """
    metrics = {}
    try:
        from pyspark import SparkContext
        sc = SparkContext.getOrCreate()

        # Executor count (subtract 1 for driver)
        try:
            status = sc._jsc.sc().getExecutorMemoryStatus()
            metrics['num_executors_at_completion'] = status.size() - 1
        except Exception:
            pass

        # Default parallelism (reflects actual cluster capacity)
        try:
            metrics['default_parallelism'] = sc.defaultParallelism
        except Exception:
            pass

        # Application ID for cross-referencing with Spark UI / YARN logs
        try:
            metrics['application_id'] = sc.applicationId
        except Exception:
            pass

        logger.info("=" * 60)
        logger.info("SPARK METRICS (post-computation)")
        logger.info("=" * 60)
        for k, v in metrics.items():
            logger.info(f"  {k}: {v}")
        logger.info("=" * 60)

    except Exception as e:
        logger.warning(f"Could not collect Spark metrics: {e}")
        metrics['error'] = str(e)

    return metrics


# =============================================================================
# Data Loading
# =============================================================================

def get_gcs_filesystem() -> gcsfs.GCSFileSystem:
    """Get GCS filesystem with requester pays enabled."""
    project = os.environ.get('GOOGLE_PROJECT', '')
    return gcsfs.GCSFileSystem(project=project, requester_pays=True)


def load_metadata(cfg: PipelineConfig, logger: logging.Logger) -> pd.DataFrame:
    """Load and filter sample metadata."""
    logger.info(f"Loading metadata from {cfg.metadata_path}")
    
    fs = get_gcs_filesystem()
    path = cfg.metadata_path.replace('gs://', '')
    
    with fs.open(path, 'r') as f:
        df = pd.read_csv(f, sep='\t')
    
    df[cfg.sample_id_col] = df[cfg.sample_id_col].astype(str)
    logger.info(f"Loaded {len(df)} samples")
    
    # Filter to requested batches
    if cfg.batches is not None:
        df = df[df[cfg.batch_col].isin(cfg.batches)]
        logger.info(f"Filtered to batches {cfg.batches}: {len(df)} samples")
    
    # Filter to requested ancestries
    if cfg.ancestries is not None:
        df = df[df[cfg.ancestry_col].isin(cfg.ancestries)]
        logger.info(f"Filtered to ancestries {cfg.ancestries}: {len(df)} samples")
    
    return df


def subsample_metadata(
    df: pd.DataFrame, 
    cfg: PipelineConfig, 
    logger: logging.Logger
) -> pd.DataFrame:
    """Subsample metadata for non-full runs."""
    if cfg.samples_per_group is None:
        logger.info("No subsampling: using all samples")
        return df
    
    logger.info(f"Subsampling to {cfg.samples_per_group} per (batch, ancestry) group")
    
    rng = np.random.RandomState(cfg.random_seed)
    groups = []
    
    for (batch, ancestry), group_df in df.groupby([cfg.batch_col, cfg.ancestry_col]):
        n = len(group_df)
        k = min(cfg.samples_per_group, n)
        if k < cfg.min_samples_per_group:
            logger.warning(f"Skipping ({batch}, {ancestry}): only {n} samples")
            continue
        idx = rng.choice(n, size=k, replace=False)
        groups.append(group_df.iloc[idx])
    
    if not groups:
        raise ValueError("No groups with sufficient samples after filtering")
    
    subsampled = pd.concat(groups, ignore_index=True)
    logger.info(f"Subsampled to {len(subsampled)} total samples")
    
    # Log group sizes
    group_counts = subsampled.groupby([cfg.batch_col, cfg.ancestry_col]).size()
    logger.info(f"Group sizes:\n{group_counts.to_string()}")
    
    return subsampled

# =============================================================================
# Hail Operations
# =============================================================================
def import_from_vcf(cfg: PipelineConfig, logger: logging.Logger):
    """Import VCFs, split multi-allelics, optionally apply ACAF filter.

    Args:
        cfg: Pipeline configuration with vcf_path, split_multi, acaf_filter settings.
        logger: Logger instance.

    Returns:
        Hail MatrixTable with biallelic variants ready for downstream processing.
    """
    hl = _init_hail()

    logger.info(f"Importing VCFs from {cfg.vcf_path}")
    mt = hl.import_vcf(
        cfg.vcf_path,
        force_bgz=True,
        reference_genome='GRCh38',
        array_elements_required=False,
    )
    logger.info(f"Imported VCF: {mt.count_cols()} samples, {mt.n_partitions()} partitions")

    if cfg.split_multi:
        mt = hl.split_multi_hts(mt)
        logger.info("Split multi-allelic sites")

    if cfg.acaf_filter:
        # INFO/AC and INFO/AF field names may vary — verify against VCF header
        mt = mt.filter_rows(
            (hl.max(mt.info.AF) > cfg.acaf_min_af) |
            (hl.max(mt.info.AC) > cfg.acaf_min_ac)
        )
        logger.info(f"Applied ACAF filter: AF>{cfg.acaf_min_af} or AC>{cfg.acaf_min_ac}")

    return mt


def load_and_filter_mt(
    cfg: PipelineConfig,
    sample_ids: List[str],
    logger: logging.Logger
) -> Tuple[Any, Dict[str, Any]]:
    """Load variant data (from VCF or pre-built MT), filter, and return interval tables.

    When data_source="vcf": imports from VCF, checks for cached MT first.
    When data_source="mt": reads pre-built MT directly (legacy fast path).

    Returns:
        Tuple of (filtered MatrixTable, dict of interval Hail Tables for annotation).
    """
    hl = _init_hail()

    if cfg.data_source == "vcf":
        # --- VCF path: check cache first, then import ---
        if cfg.cache_mt and not cfg.force_reimport:
            try:
                mt = hl.read_matrix_table(cfg.cache_mt_path)
                logger.info(f"Loaded cached MT from {cfg.cache_mt_path}")
            except Exception:
                logger.info("No cached MT found, importing from VCF")
                mt = import_from_vcf(cfg, logger)
                if cfg.cache_mt:
                    logger.info(f"Writing cache to {cfg.cache_mt_path}")
                    mt.write(cfg.cache_mt_path, overwrite=True)
                    mt = hl.read_matrix_table(cfg.cache_mt_path)
                    logger.info(f"Cached imported MT to {cfg.cache_mt_path}")
        else:
            mt = import_from_vcf(cfg, logger)
            if cfg.cache_mt:
                logger.info(f"Writing cache to {cfg.cache_mt_path} (force_reimport={cfg.force_reimport})")
                mt.write(cfg.cache_mt_path, overwrite=True)
                mt = hl.read_matrix_table(cfg.cache_mt_path)
                logger.info(f"Cached imported MT to {cfg.cache_mt_path}")
        logger.info(f"MT ready: {mt.n_partitions()} partitions (row/col counts deferred)")
    elif cfg.data_source == "mt":
        logger.info(f"Loading pre-built MatrixTable from {cfg.hail_mt_path}")
        mt = hl.read_matrix_table(cfg.hail_mt_path)
        logger.info(f"MT loaded: {mt.n_partitions()} partitions (row/col counts deferred)")
    else:
        raise ValueError(f"Unknown data_source: {cfg.data_source!r} (expected 'vcf' or 'mt')")

    # --- Drop unused entry fields early (only GT and FT are needed) ---
    mt = mt.select_entries(mt.GT, mt.FT)
    logger.info("Selected entry fields: GT, FT (dropped unused fields)")

    # --- Filter columns (samples) BEFORE rows for efficiency ---
    sample_set = hl.literal(set(sample_ids))
    mt = mt.filter_cols(sample_set.contains(mt.s))
    logger.info(f"Filtering to {len(sample_ids)} requested samples (lazy)")

    # --- Load intervals once (used for both partition pruning and annotation) ---
    interval_tables = {}
    all_intervals_for_pruning = []

    if cfg.interval_dict:
        for name, path in cfg.interval_dict.items():
            logger.info(f"Loading interval list: {name}")
            is_gzipped = path.endswith('.gz')

            try:
                if '.interval_list' in path:
                    interval_table = hl.import_locus_intervals(
                        path, reference_genome='GRCh38',
                        force_bgz=is_gzipped, filter=r'^(#|@)'
                    )
                else:
                    interval_table = hl.import_bed(
                        path, reference_genome='GRCh38',
                        force_bgz=is_gzipped,
                        skip_invalid_intervals=True,
                        filter=r'^#'
                    )

                # Keep the table for annotation later
                interval_tables[name] = interval_table

                # Collect intervals for partition pruning
                intervals = interval_table.interval.collect()
                logger.info(f"  {name}: {len(intervals)} intervals")

                # Subsample HighConf_Genome if configured
                if cfg.pruning_subsample_n is not None and len(intervals) > cfg.pruning_subsample_n:
                    rng = np.random.RandomState(cfg.random_seed)
                    idx = rng.choice(len(intervals), size=cfg.pruning_subsample_n, replace=False)
                    intervals = [intervals[i] for i in sorted(idx)]
                    logger.info(f"  Subsampled {name} to {len(intervals)} intervals for pruning")

                all_intervals_for_pruning.append(intervals)

            except Exception as e:
                logger.warning(f"  Failed to load {name}: {e}")
                continue

        # Flatten all intervals and apply partition pruning
        flat_intervals = [iv for sublist in all_intervals_for_pruning for iv in sublist]
        logger.info(f"Total: {len(flat_intervals)} intervals for partition pruning")

        if flat_intervals:
            partitions_before = mt.n_partitions()
            mt = hl.filter_intervals(mt, flat_intervals, keep=True)
            partitions_after = mt.n_partitions()
            logger.info(f"Partition pruning: {partitions_before} → {partitions_after} partitions ({partitions_before - partitions_after} pruned)")

    # For smoke/dev mode, additionally filter to specific chromosomes
    if cfg.mode == 'smoke':
        mt = mt.filter_rows(mt.locus.contig == 'chr17')
        logger.info("SMOKE MODE: Filtered to chr17 only")
    elif cfg.mode == 'dev':
        dev_chroms = {'chr7', 'chr17', 'chr13'}
        mt = mt.filter_rows(hl.literal(dev_chroms).contains(mt.locus.contig))
        logger.info(f"DEV MODE: Filtered to {dev_chroms}")

    # Apply variant filters
    if cfg.filter_to_pass:
        mt = mt.filter_rows(hl.is_missing(mt.filters) | (mt.filters.length() == 0))
        logger.info("Filtered to PASS variants")

    # Apply genotype filter (FT field)
    if cfg.filter_FT_pass:
        mt = mt.annotate_entries(
            GT=hl.or_missing(
                hl.is_missing(mt.FT) | mt.FT.contains("PASS"),
                mt.GT
            )
        )
        logger.info("Set GT to missing where FT is not PASS")

    return mt, interval_tables


# NOTE: load_interval_tables() has been merged into load_and_filter_mt()
# to avoid loading interval BED files twice. Interval tables are now
# returned as the second element of load_and_filter_mt()'s return tuple.

def annotate_mt_with_metadata(
    mt,
    metadata: pd.DataFrame,
    cfg: PipelineConfig,
    logger: logging.Logger
):
    """Annotate MatrixTable columns with batch and ancestry labels."""
    hl = _init_hail()
    
    logger.info("Annotating MT with sample metadata")
    
    # Create Hail table from metadata
    meta_ht = hl.Table.from_pandas(metadata[[
        cfg.sample_id_col, cfg.batch_col, cfg.ancestry_col
    ]])
    meta_ht = meta_ht.key_by(cfg.sample_id_col)
    
    # Annotate MT columns
    mt = mt.annotate_cols(
        batch=meta_ht[mt.s][cfg.batch_col],
        ancestry=meta_ht[mt.s][cfg.ancestry_col]
    )
    
    return mt



def compute_variant_stats(
    mt,
    interval_tables: Dict[str, Any],
    cfg: PipelineConfig,
    logger: logging.Logger
):
    """Compute per-sample variant statistics within each interval.

    Uses hl.agg.explode + hl.agg.group_by to produce compact JVM bytecode.
    This avoids the 64KB method size limit that occurs when unrolling
    N intervals × 8 aggregators into separate hl.agg.filter blocks.

    Args:
        mt: MatrixTable with GT entry field and locus row key.
        interval_tables: Dict mapping interval name to Hail Table (keyed by locus interval).

    Returns:
        MT with `interval_stats` column: dict<str, struct<8 counters>> per sample.
    """
    hl = _init_hail()

    logger.info(f"Computing variant stats for {len(interval_tables)} intervals (explode+group_by)")

    # For each variant, build array of matching interval names
    interval_flags = hl.array([
        hl.if_else(hl.is_defined(table[mt.locus]), name, hl.missing(hl.tstr))
        for name, table in interval_tables.items()
    ]).filter(hl.is_defined)

    # Pre-compute shared conditions (referenced once in bytecode template)
    is_ti = hl.is_transition(mt.alleles[0], mt.alleles[1])
    is_tv = hl.is_transversion(mt.alleles[0], mt.alleles[1])
    is_ins = hl.is_insertion(mt.alleles[0], mt.alleles[1])
    is_del = hl.is_deletion(mt.alleles[0], mt.alleles[1])
    is_het = mt.GT.is_het()
    is_hom = mt.GT.is_hom_var()

    # Single aggregation — group_by dispatches at runtime, not compile time
    mt = mt.annotate_cols(
        interval_stats=hl.agg.explode(
            lambda iv: hl.agg.group_by(iv, hl.struct(
                snp_ti_het=hl.agg.count_where(is_ti & is_het),
                snp_tv_het=hl.agg.count_where(is_tv & is_het),
                snp_ti_hom=hl.agg.count_where(is_ti & is_hom),
                snp_tv_hom=hl.agg.count_where(is_tv & is_hom),
                indel_ins_het=hl.agg.count_where(is_ins & is_het),
                indel_del_het=hl.agg.count_where(is_del & is_het),
                indel_ins_hom=hl.agg.count_where(is_ins & is_hom),
                indel_del_hom=hl.agg.count_where(is_del & is_hom),
            )),
            interval_flags
        )
    )

    return mt

def extract_sample_stats(mt, cfg: PipelineConfig, logger: logging.Logger) -> pd.DataFrame:
    """Extract per-sample statistics as a pandas DataFrame.

    Unpacks the `interval_stats` dict<str, struct<8 counters>> produced by
    compute_variant_stats into flat columns (e.g., ACMG59_snp_ti_het).
    """
    hl = _init_hail()

    logger.info("Extracting sample statistics to DataFrame")

    # Get column table and unkey to avoid key conflicts
    cols_ht = mt.cols()
    cols_ht = cols_ht.key_by()  # Remove key to allow selecting 's'

    # Flatten the dict<str, struct> into flat columns per interval
    interval_names = list(cfg.interval_dict.keys())

    # Default struct for intervals with no variants for a sample
    zero_stats = hl.struct(
        snp_ti_het=hl.int64(0), snp_tv_het=hl.int64(0),
        snp_ti_hom=hl.int64(0), snp_tv_hom=hl.int64(0),
        indel_ins_het=hl.int64(0), indel_del_het=hl.int64(0),
        indel_ins_hom=hl.int64(0), indel_del_hom=hl.int64(0),
    )

    # Build selection
    select_exprs = {
        's': cols_ht.s,
        'batch': cols_ht.batch,
        'ancestry': cols_ht.ancestry,
    }

    stat_fields = [
        'snp_ti_het', 'snp_tv_het', 'snp_ti_hom', 'snp_tv_hom',
        'indel_ins_het', 'indel_del_het', 'indel_ins_hom', 'indel_del_hom',
    ]

    for interval in interval_names:
        stats = hl.coalesce(cols_ht.interval_stats.get(interval), zero_stats)
        for f in stat_fields:
            select_exprs[f"{interval}_{f}"] = stats[f]

    result_ht = cols_ht.select(**select_exprs)

    # Convert to pandas
    df = result_ht.to_pandas()

    logger.info(f"Extracted stats for {len(df)} samples")

    return df

# =============================================================================
# Summary Statistics
# =============================================================================

def compute_derived_metrics(df: pd.DataFrame, cfg: PipelineConfig) -> pd.DataFrame:
    """Compute derived metrics (Ti/Tv ratio, del/ins ratio, etc.)."""
    df = df.copy()
    
    for interval in cfg.interval_dict.keys():
        prefix = interval
        
        # SNP counts
        df[f"{prefix}_snp_het"] = df[f"{prefix}_snp_ti_het"] + df[f"{prefix}_snp_tv_het"]
        df[f"{prefix}_snp_hom"] = df[f"{prefix}_snp_ti_hom"] + df[f"{prefix}_snp_tv_hom"]
        df[f"{prefix}_snp_total"] = df[f"{prefix}_snp_het"] + df[f"{prefix}_snp_hom"]
        
        # Ti/Tv ratios
        tv_het = df[f"{prefix}_snp_tv_het"].replace(0, np.nan)
        tv_hom = df[f"{prefix}_snp_tv_hom"].replace(0, np.nan)
        tv_total = (df[f"{prefix}_snp_tv_het"] + df[f"{prefix}_snp_tv_hom"]).replace(0, np.nan)
        
        df[f"{prefix}_titv_het"] = df[f"{prefix}_snp_ti_het"] / tv_het
        df[f"{prefix}_titv_hom"] = df[f"{prefix}_snp_ti_hom"] / tv_hom
        df[f"{prefix}_titv_total"] = (df[f"{prefix}_snp_ti_het"] + df[f"{prefix}_snp_ti_hom"]) / tv_total
        
        # Indel counts
        df[f"{prefix}_indel_het"] = df[f"{prefix}_indel_ins_het"] + df[f"{prefix}_indel_del_het"]
        df[f"{prefix}_indel_hom"] = df[f"{prefix}_indel_ins_hom"] + df[f"{prefix}_indel_del_hom"]
        df[f"{prefix}_indel_total"] = df[f"{prefix}_indel_het"] + df[f"{prefix}_indel_hom"]
        
        # Del/Ins ratios
        ins_het = df[f"{prefix}_indel_ins_het"].replace(0, np.nan)
        ins_hom = df[f"{prefix}_indel_ins_hom"].replace(0, np.nan)
        ins_total = (df[f"{prefix}_indel_ins_het"] + df[f"{prefix}_indel_ins_hom"]).replace(0, np.nan)
        
        df[f"{prefix}_delins_het"] = df[f"{prefix}_indel_del_het"] / ins_het
        df[f"{prefix}_delins_hom"] = df[f"{prefix}_indel_del_hom"] / ins_hom
        df[f"{prefix}_delins_total"] = (df[f"{prefix}_indel_del_het"] + df[f"{prefix}_indel_del_hom"]) / ins_total
    
    return df


def compute_group_summaries(
    df: pd.DataFrame,
    cfg: PipelineConfig,
    logger: logging.Logger
) -> pd.DataFrame:
    """Compute mean/std for each metric by batch (and optionally ancestry)."""
    
    # Identify metric columns (everything except s, batch, ancestry)
    id_cols = ['s', 'batch', 'ancestry']
    metric_cols = [c for c in df.columns if c not in id_cols]
    
    if cfg.analysis_strategy == "within_ancestry":
        group_cols = ['batch', 'ancestry']
    else:
        group_cols = ['batch']
    
    logger.info(f"Computing summaries grouped by {group_cols}")
    
    # Compute aggregations
    agg_funcs = ['mean', 'std', 'count']
    summary = df.groupby(group_cols)[metric_cols].agg(agg_funcs)
    
    # Flatten column names
    summary.columns = ['_'.join(col).strip() for col in summary.columns]
    summary = summary.reset_index()
    
    return summary


def compute_pairwise_comparisons(
    df: pd.DataFrame,
    cfg: PipelineConfig,
    logger: logging.Logger
) -> pd.DataFrame:
    """Compute pairwise batch comparisons (effect sizes and p-values)."""
    from scipy import stats as scipy_stats
    
    batches = df['batch'].unique()
    
    if cfg.analysis_strategy == "within_ancestry":
        ancestries = df['ancestry'].unique()
    else:
        ancestries = [None]
    
    # Focus on key derived metrics
    key_metrics = [c for c in df.columns if c.endswith(('_snp_total', '_indel_total', '_titv_total', '_delins_total'))]
    
    logger.info(f"Computing comparisons for metrics: {key_metrics}")
    
    results = []
    
    for ancestry in ancestries:
        if ancestry is not None:
            sub_df = df[df['ancestry'] == ancestry]
        else:
            sub_df = df
        
        for i, batch_x in enumerate(batches):
            for batch_y in list(batches)[i+1:]:
                x_data = sub_df[sub_df['batch'] == batch_x]
                y_data = sub_df[sub_df['batch'] == batch_y]
                
                if len(x_data) < cfg.min_samples_per_group or len(y_data) < cfg.min_samples_per_group:
                    continue
                
                for metric in key_metrics:
                    # Convert to numpy float64 arrays explicitly
                    x_vals = pd.to_numeric(x_data[metric], errors='coerce').astype('float64').dropna().values
                    y_vals = pd.to_numeric(y_data[metric], errors='coerce').astype('float64').dropna().values
                    
                    if len(x_vals) < 10 or len(y_vals) < 10:
                        logger.warning(f"Skipping {metric} for {batch_x} vs {batch_y}: insufficient data")
                        continue
                    
                    mean_x, std_x = x_vals.mean(), x_vals.std()
                    mean_y, std_y = y_vals.mean(), y_vals.std()
                    n_x, n_y = len(x_vals), len(y_vals)

                    # Cohen's d with pooled standard deviation
                    pooled_std = np.sqrt(
                        ((n_x - 1) * std_x**2 + (n_y - 1) * std_y**2) / (n_x + n_y - 2)
                    )
                    cohens_d = (mean_x - mean_y) / pooled_std if pooled_std > 0 else np.nan
                    
                    # Welch's t-test
                    try:
                        t_stat, p_val = scipy_stats.ttest_ind(x_vals, y_vals, equal_var=False)
                    except Exception as e:
                        logger.warning(f"t-test failed for {metric}: {e}")
                        p_val = np.nan
                    
                    results.append({
                        'ancestry': ancestry,
                        'batch_x': batch_x,
                        'batch_y': batch_y,
                        'metric': metric,
                        'mean_x': mean_x,
                        'std_x': std_x,
                        'n_x': len(x_vals),
                        'mean_y': mean_y,
                        'std_y': std_y,
                        'n_y': len(y_vals),
                        'cohens_d': cohens_d,
                        'p_value': p_val,
                        'neg_log10_p': -np.log10(p_val) if p_val and p_val > 0 else np.nan,
                    })
    
    comparison_df = pd.DataFrame(results)
    logger.info(f"Computed {len(comparison_df)} pairwise comparisons")
    
    return comparison_df



# =============================================================================
# Output
# =============================================================================

def save_results(
    sample_stats: pd.DataFrame,
    group_summaries: pd.DataFrame,
    comparisons: pd.DataFrame,
    cfg: PipelineConfig,
    logger: logging.Logger
) -> Dict[str, str]:
    """Save all results to GCS."""
    fs = get_gcs_filesystem()
    output_dir = cfg.output_dir.replace('gs://', '')
    
    outputs = {}
    
    # Save sample stats (optionally)
    if cfg.export_per_sample or cfg.mode in ('smoke', 'dev'):
        path = f"{output_dir}/sample_stats.tsv"
        logger.info(f"Saving sample stats to gs://{path}")
        with fs.open(path, 'w') as f:
            sample_stats.to_csv(f, sep='\t', index=False)
        outputs['sample_stats'] = f"gs://{path}"
    
    # Save group summaries
    path = f"{output_dir}/group_summaries.tsv"
    logger.info(f"Saving group summaries to gs://{path}")
    with fs.open(path, 'w') as f:
        group_summaries.to_csv(f, sep='\t', index=False)
    outputs['group_summaries'] = f"gs://{path}"
    
    # Save comparisons
    path = f"{output_dir}/pairwise_comparisons.tsv"
    logger.info(f"Saving pairwise comparisons to gs://{path}")
    with fs.open(path, 'w') as f:
        comparisons.to_csv(f, sep='\t', index=False)
    outputs['comparisons'] = f"gs://{path}"
    
    # Save config
    path = f"{output_dir}/config.json"
    logger.info(f"Saving config to gs://{path}")
    config_dict = asdict(cfg)
    with fs.open(path, 'w') as f:
        json.dump(config_dict, f, indent=2, default=str)
    outputs['config'] = f"gs://{path}"
    
    return outputs

# =============================================================================
# Main Pipeline
# =============================================================================

def run_pipeline(cfg: PipelineConfig) -> Dict[str, Any]:
    """
    Run the complete batch effect analysis pipeline.

    Args:
        cfg: Pipeline configuration

    Returns:
        Dictionary with output paths, summary statistics, and timing data
    """
    logger = setup_logging(cfg.run_name)
    step_timings = {}
    pipeline_start = time.time()

    logger.info("=" * 60)
    logger.info(f"Starting batch effect pipeline: {cfg.run_name}")
    logger.info(f"Mode: {cfg.mode}")
    logger.info(f"Data source: {cfg.data_source}")
    if cfg.data_source == "vcf":
        logger.info(f"VCF path: {cfg.vcf_path}")
        logger.info(f"Cache MT: {cfg.cache_mt} (path: {cfg.cache_mt_path})")
        logger.info(f"Force reimport: {cfg.force_reimport}")
    else:
        logger.info(f"MT path: {cfg.hail_mt_path}")
    logger.info(f"Pruning subsample threshold: {cfg.pruning_subsample_n}")
    logger.info("=" * 60)

    # Log cluster info after Hail init
    _init_hail()
    cluster_info = log_cluster_info(logger)

    # Step 1: Load and filter metadata
    logger.info("\n[Step 1/6] Loading metadata...")
    t0 = time.time()
    metadata = load_metadata(cfg, logger)
    metadata = subsample_metadata(metadata, cfg, logger)
    sample_ids = metadata[cfg.sample_id_col].tolist()
    step_timings['step1_load_metadata'] = round(time.time() - t0, 1)
    logger.info(f"  Step 1 took {step_timings['step1_load_metadata']}s")

    # Step 2: Load variant data (VCF or pre-built MT) and intervals
    logger.info(f"\n[Step 2/6] Loading variant data (source: {cfg.data_source})...")
    t0 = time.time()
    mt, interval_tables = load_and_filter_mt(cfg, sample_ids, logger)
    step_timings['step2_load_and_filter_mt'] = round(time.time() - t0, 1)
    logger.info(f"  Step 2 took {step_timings['step2_load_and_filter_mt']}s")

    # Step 3: Annotate with metadata
    logger.info("\n[Step 3/6] Annotating with batch/ancestry labels...")
    t0 = time.time()
    mt = annotate_mt_with_metadata(mt, metadata, cfg, logger)
    step_timings['step3_annotate_metadata'] = round(time.time() - t0, 1)
    logger.info(f"  Step 3 took {step_timings['step3_annotate_metadata']}s")

    # Step 4: Compute variant stats (explode+group_by — compact bytecode)
    logger.info("\n[Step 4/6] Computing variant statistics...")
    t0 = time.time()
    mt = compute_variant_stats(mt, interval_tables, cfg, logger)
    step_timings['step4_compute_stats'] = round(time.time() - t0, 1)
    logger.info(f"  Step 4 took {step_timings['step4_compute_stats']}s (graph construction only)")

    # Step 5: Extract and process results (this triggers the actual Spark computation)
    logger.info("\n[Step 5/6] Extracting results (triggers Spark execution)...")
    t0 = time.time()
    sample_stats = extract_sample_stats(mt, cfg, logger)
    step_timings['step5a_extract_to_pandas'] = round(time.time() - t0, 1)
    logger.info(f"  Step 5a (to_pandas) took {step_timings['step5a_extract_to_pandas']}s — THIS IS THE MAIN COMPUTE COST")

    t0 = time.time()
    sample_stats = compute_derived_metrics(sample_stats, cfg)
    group_summaries = compute_group_summaries(sample_stats, cfg, logger)
    comparisons = compute_pairwise_comparisons(sample_stats, cfg, logger)
    step_timings['step5b_pandas_comparisons'] = round(time.time() - t0, 1)
    logger.info(f"  Step 5b (pandas) took {step_timings['step5b_pandas_comparisons']}s")

    # Collect Spark metrics after the main computation
    spark_metrics = collect_spark_metrics(logger)

    # Step 6: Save outputs
    logger.info("\n[Step 6/6] Saving results...")
    t0 = time.time()
    outputs = save_results(sample_stats, group_summaries, comparisons, cfg, logger)
    step_timings['step6_save_results'] = round(time.time() - t0, 1)
    logger.info(f"  Step 6 took {step_timings['step6_save_results']}s")

    step_timings['total_pipeline'] = round(time.time() - pipeline_start, 1)

    # Save timing and metrics to output
    fs = get_gcs_filesystem()
    output_dir = cfg.output_dir.replace('gs://', '')
    timing_path = f"{output_dir}/timing_metrics.json"
    timing_data = {
        'step_timings': step_timings,
        'spark_metrics': spark_metrics,
        'cluster_info': cluster_info,
        'config_summary': {
            'mode': cfg.mode,
            'data_source': cfg.data_source,
            'samples_per_group': cfg.samples_per_group,
            'n_intervals': len(cfg.interval_dict) if cfg.interval_dict else 0,
            'interval_names': list(cfg.interval_dict.keys()) if cfg.interval_dict else [],
            'pruning_subsample_n': cfg.pruning_subsample_n,
            'n_samples_requested': len(sample_ids),
            'n_samples_extracted': len(sample_stats),
        },
    }
    with fs.open(timing_path, 'w') as f:
        json.dump(timing_data, f, indent=2, default=str)
    outputs['timing_metrics'] = f"gs://{timing_path}"
    logger.info(f"Saved timing metrics to gs://{timing_path}")

    # Print timing summary
    logger.info("\n" + "=" * 60)
    logger.info("TIMING SUMMARY")
    logger.info("=" * 60)
    for step, elapsed in step_timings.items():
        logger.info(f"  {step}: {elapsed}s")
    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"Output directory: {cfg.output_dir}")
    logger.info("=" * 60)

    return {
        'outputs': outputs,
        'sample_stats': sample_stats,
        'group_summaries': group_summaries,
        'comparisons': comparisons,
        'config': cfg,
        'step_timings': step_timings,
        'spark_metrics': spark_metrics,
    }


# =============================================================================
# CLI Entry Point
# =============================================================================

def cli():
    import argparse
    
    parser = argparse.ArgumentParser(description='Batch Effect Analysis Pipeline')
    parser.add_argument('--mode', choices=['smoke', 'dev', 'medium', 'full'], default='dev')
    parser.add_argument('--batches', nargs='+', help='Batches to include')
    parser.add_argument('--ancestries', nargs='+', help='Ancestries to include')
    parser.add_argument('--strategy', choices=['within_ancestry', 'pooled'], default='within_ancestry')
    parser.add_argument('--output-dir', help='Output directory')
    parser.add_argument('--data-source', choices=['vcf', 'mt'], default='vcf',
                        help='Data source: "vcf" (primary) or "mt" (pre-built MatrixTable)')
    parser.add_argument('--vcf-path', help='VCF path/glob (default: $WGS_ACAF_THRESHOLD_VCF_PATH)')
    parser.add_argument('--mt-path', help='Pre-built MT path (for --data-source mt)')
    parser.add_argument('--no-cache', action='store_true', help='Disable MT caching of VCF import')
    parser.add_argument('--force-reimport', action='store_true',
                        help='Ignore cached MT, re-import from VCF')

    args = parser.parse_args()

    cfg = PipelineConfig(
        mode=args.mode,
        batches=args.batches,
        ancestries=args.ancestries,
        analysis_strategy=args.strategy,
        output_dir=args.output_dir,
        data_source=args.data_source,
        vcf_path=args.vcf_path,
        hail_mt_path=args.mt_path,
        cache_mt=not args.no_cache,
        force_reimport=args.force_reimport,
    )
    
    results = run_pipeline(cfg)
    print(f"\nOutputs saved to: {cfg.output_dir}")