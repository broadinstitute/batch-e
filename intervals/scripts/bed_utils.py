"""Shared BED-file loading utilities for interval analysis scripts."""

import pandas as pd

GRCH38_TOTAL_BP = 3_088_269_832

KARYOTYPE_ORDER = {f"chr{i}": i for i in range(1, 23)}
KARYOTYPE_ORDER["chrX"] = 23
KARYOTYPE_ORDER["chrY"] = 24

# Maps BED filenames (stem) to human-readable interval names.
DEFAULT_NAME_MAP = {
    "acmg59_allofus_19dec2019.GRC38.wGenes.NEW": "ACMG59",
    "GRCh38_lowmappabilityall": "Low_Mappability",
    "GRCh38_gc85_slop50": "GC_gt_85",
    "GRCh38_gclt25_merged": "GC_lt_25",
    "giab_highconf_wgs_calling_regions_hg38_intersection": "HighConf_Genome",
}


def load_bed(path: str) -> pd.DataFrame:
    """Parse a BED file into a DataFrame with [chrom, start, end, width, chrom_order].

    Skips lines starting with '#'. Handles 3+ column BED files (extra cols dropped).
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": int, "end": int},
    )
    df["width"] = df["end"] - df["start"]
    df["chrom_order"] = df["chrom"].map(KARYOTYPE_ORDER)
    # Drop non-standard contigs (alt, random, Un, etc.)
    df = df.dropna(subset=["chrom_order"]).copy()
    df["chrom_order"] = df["chrom_order"].astype(int)
    return df


def format_bp(bp: int) -> str:
    """Human-readable base-pair formatting: '2,386 Mbp', '545 Kbp', '312 bp'."""
    if bp >= 1_000_000_000:
        return f"{bp / 1_000_000_000:,.2f} Gbp"
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:,.1f} Mbp"
    if bp >= 1_000:
        return f"{bp / 1_000:,.1f} Kbp"
    return f"{bp:,} bp"
