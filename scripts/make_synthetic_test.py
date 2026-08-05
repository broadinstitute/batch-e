#!/usr/bin/env python3
"""
Generate a tiny synthetic dataset for end-to-end WDL smoke tests.

Produces under --gcs-prefix:
    synthetic.mt/                  Hail MatrixTable (300 samples, ~5K variants)
    ancestry.tsv                   research_id, ancestry_pred_other
    comparison.tsv                 research_id, sample_source
    intervals/*.bed.gz             5 tiny BEDs matching the MT's chr20 window

The MT covers chr20:1,000,000-2,000,000. Sample matrix is
2 ancestries (EUR, AFR) x 2 sources (blood, saliva) x 75 samples = 300.
A per-group AF shift is baked into GT sampling so downstream comparisons
produce non-trivial Cohen's d.

Run via hailrunner (same pattern as batch_e.py):

    hailrunner run \\
        --staging-bucket gs://fc-secure-.../ \\
        --script https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/scripts/make_synthetic_test.py \\
        --workers 2 --worker-type n1-standard-4 --driver-type n1-standard-4 \\
        -- --gcs-prefix gs://fc-secure-.../batch_e_synth_v1

Runs equally well on a laptop with hail installed:

    python scripts/make_synthetic_test.py \\
        --gcs-prefix gs://fc-secure-.../batch_e_synth_v1
"""

from __future__ import annotations

import argparse
import gzip
import io

import hail as hl
import numpy as np
import pandas as pd


CONTIG = "chr20"
WINDOW_START = 1_000_000
WINDOW_END = 2_000_000

N_VARIANTS = 5_000
N_PARTITIONS = 8

ANCESTRIES = ["EUR", "AFR"]
SOURCES = ["blood", "saliva"]
PER_CELL = 75  # 2 x 2 x 75 = 300 samples

SEED = 42


def build_sample_table() -> pd.DataFrame:
    """Deterministic 300-row sample metadata with a per-group AF shift."""
    rows = []
    idx = 0
    for ancestry in ANCESTRIES:
        for source in SOURCES:
            # Per-cell allele-frequency shift. Blood carriers get a small
            # positive shift, saliva a negative one; EUR/AFR differ slightly
            # so the ancestry stratification also produces visible structure.
            shift = 0.0
            shift += 0.010 if source == "blood" else -0.010
            shift += 0.005 if ancestry == "EUR" else -0.005
            for _ in range(PER_CELL):
                idx += 1
                rows.append({
                    "research_id": f"SYNTH_{idx:06d}",
                    "ancestry_pred_other": ancestry,
                    "sample_source": source,
                    "af_shift": shift,
                })
    df = pd.DataFrame(rows)
    # Shuffle so col ordering doesn't line up with ancestry/source blocks.
    df = df.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    return df


def build_intervals() -> dict[str, list[tuple[str, int, int]]]:
    """Five pretend interval classes carving up the synthetic window."""
    rng = np.random.RandomState(SEED + 1)

    def random_slabs(start: int, end: int, n: int, width: int) -> list[tuple[str, int, int]]:
        step = (end - start) // n
        out = []
        for i in range(n):
            s = start + i * step + rng.randint(0, max(1, step - width))
            e = min(end, s + width)
            out.append((CONTIG, s, e))
        return out

    return {
        "ACMG59":          random_slabs(1_050_000, 1_150_000, 20, 500),
        "Low_Mappability": random_slabs(1_200_000, 1_400_000, 30, 1000),
        "GC_gt_85":        random_slabs(1_450_000, 1_550_000, 15, 300),
        "GC_lt_25":        random_slabs(1_600_000, 1_700_000, 15, 800),
        "HighConf_Genome": [(CONTIG, WINDOW_START, WINDOW_END)],
    }


def write_tsv(df: pd.DataFrame, path: str) -> None:
    with hl.hadoop_open(path, "w") as f:
        df.to_csv(f, sep="\t", index=False)


def write_bed_gz(records: list[tuple[str, int, int]], path: str) -> None:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        for chrom, start, end in records:
            gz.write(f"{chrom}\t{start}\t{end}\n".encode())
    with hl.hadoop_open(path, "wb") as f:
        f.write(buf.getvalue())


def build_matrix_table(sample_df: pd.DataFrame, mt_path: str) -> None:
    hl.reset_global_randomness_seed(SEED)

    mt = hl.utils.range_matrix_table(
        n_rows=N_VARIANTS, n_cols=len(sample_df), n_partitions=N_PARTITIONS
    )

    # --- Row fields: locus, alleles, filters ---
    pos = WINDOW_START + mt.row_idx * ((WINDOW_END - WINDOW_START) // N_VARIANTS)
    allele_class = mt.row_idx % 5
    ref = (hl.switch(allele_class)
           .when(0, "A").when(1, "A").when(2, "C").when(3, "A").when(4, "AT")
           .or_missing())
    alt = (hl.switch(allele_class)
           .when(0, "G").when(1, "C").when(2, "T").when(3, "AT").when(4, "A")
           .or_missing())
    # 5 % of rows carry a non-PASS filter tag so filter_to_pass has real work.
    filters = hl.if_else(
        mt.row_idx % 20 == 0,
        hl.set(["RF"]),
        hl.empty_set(hl.tstr),
    )
    mt = mt.annotate_rows(
        locus=hl.locus(CONTIG, pos),
        alleles=hl.array([ref, alt]),
        filters=filters,
    )
    mt = mt.key_rows_by("locus", "alleles")

    # --- Col fields: s (sample id), ancestry, sample_source, af_shift ---
    ids = hl.literal(sample_df["research_id"].tolist())
    ancestries = hl.literal(sample_df["ancestry_pred_other"].tolist())
    sources = hl.literal(sample_df["sample_source"].tolist())
    shifts = hl.literal(sample_df["af_shift"].tolist())
    mt = mt.annotate_cols(
        s=ids[mt.col_idx],
        ancestry=ancestries[mt.col_idx],
        sample_source=sources[mt.col_idx],
        af_shift=shifts[mt.col_idx],
    )
    mt = mt.key_cols_by("s")

    # --- Entry fields: GT, FT ---
    # Base AF 5 % plus the per-cell shift; independent Bernoulli per haplotype
    # yields a diploid genotype with realistic het/hom mix.
    p = hl.max(0.001, hl.min(0.999, 0.05 + mt.af_shift))
    a1 = hl.int(hl.rand_bool(p))
    a2 = hl.int(hl.rand_bool(p))
    mt = mt.annotate_entries(
        GT=hl.call(a1, a2),
        FT=hl.if_else(hl.rand_bool(0.95), hl.set(["PASS"]), hl.set(["LowGQ"])),
    )

    mt = mt.drop("af_shift")

    mt.write(mt_path, overwrite=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--gcs-prefix",
        required=True,
        help="Destination prefix, e.g. gs://bucket/batch_e_synth_v1",
    )
    args = ap.parse_args()

    hl.init(default_reference="GRCh38", idempotent=True)

    prefix = args.gcs_prefix.rstrip("/")
    mt_path = f"{prefix}/synthetic.mt"
    ancestry_path = f"{prefix}/ancestry.tsv"
    comparison_path = f"{prefix}/comparison.tsv"
    intervals_prefix = f"{prefix}/intervals"

    print(f"Building sample table ({len(ANCESTRIES) * len(SOURCES) * PER_CELL} samples)")
    sample_df = build_sample_table()

    print(f"Uploading ancestry TSV -> {ancestry_path}")
    write_tsv(sample_df[["research_id", "ancestry_pred_other"]], ancestry_path)

    print(f"Uploading comparison TSV -> {comparison_path}")
    write_tsv(sample_df[["research_id", "sample_source"]], comparison_path)

    print("Building interval BEDs")
    intervals = build_intervals()
    for name, records in intervals.items():
        gcs_path = f"{intervals_prefix}/{name}.bed.gz"
        print(f"  {name}: {len(records)} regions -> {gcs_path}")
        write_bed_gz(records, gcs_path)

    print(f"Building MatrixTable ({N_VARIANTS} variants) -> {mt_path}")
    build_matrix_table(sample_df, mt_path)

    print("\nDone.")
    print(f"input_path:     {mt_path}")
    print(f"ancestry_tsv:   {ancestry_path}")
    print(f"comparison_tsv: {comparison_path}")
    print("intervals:")
    for name in intervals:
        print(f"  {name}={intervals_prefix}/{name}.bed.gz")


if __name__ == "__main__":
    main()
