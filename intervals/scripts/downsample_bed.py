#!/usr/bin/env python3
"""Analyze and downsample BED interval files.

Two subcommands:

    # Analyze all BED files in a directory
    python downsample_bed.py probe --interval-dir intervals/

    # Downsample one BED file to match another's coverage
    python downsample_bed.py sample \
        --input big.bed --match-file ref.bed --output out.bed

Both produce self-contained Markdown reports with base64-embedded figures
suitable for pandoc PDF conversion.
"""

import argparse
import base64
import io
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from bed_utils import (
    DEFAULT_NAME_MAP,
    GRCH38_TOTAL_BP,
    KARYOTYPE_ORDER,
    format_bp,
    load_bed,
)

DEFAULT_DPI = 150


# ---------------------------------------------------------------------------
# Figure utilities
# ---------------------------------------------------------------------------

def fig_to_base64(fig: plt.Figure, dpi: int = DEFAULT_DPI) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def fig_to_markdown(fig: plt.Figure, alt: str = "", dpi: int = DEFAULT_DPI) -> str:
    """Return a Markdown image tag with an inline base64-encoded PNG."""
    return f"![{alt}](data:image/png;base64,{fig_to_base64(fig, dpi)})"


# ---------------------------------------------------------------------------
# BED discovery / loading
# ---------------------------------------------------------------------------

def discover_beds(interval_dir: str) -> dict[str, str]:
    """Auto-discover BED files and map to human-readable names.

    Skips files whose stem ends with '_downsampled' (our own outputs).
    """
    beds = {}
    for f in sorted(Path(interval_dir).iterdir()):
        if not f.name.endswith(".bed"):
            continue
        if f.stem.endswith("_downsampled"):
            continue
        name = DEFAULT_NAME_MAP.get(f.stem, f.stem)
        beds[name] = str(f)
    return beds


def parse_bed_specs(specs: list[str]) -> dict[str, str]:
    """Parse NAME=path pairs from CLI."""
    beds = {}
    for spec in specs:
        if "=" not in spec:
            print(f"Error: --beds expects NAME=path, got: {spec}", file=sys.stderr)
            sys.exit(1)
        name, path = spec.split("=", 1)
        beds[name] = path
    return beds


def load_all(beds: dict[str, str]) -> dict[str, pd.DataFrame]:
    """Load all BED files, return {name: DataFrame}."""
    loaded = {}
    for name, path in beds.items():
        df = load_bed(path)
        print(f"  Loaded {name}: {len(df):,} intervals from {path}")
        loaded[name] = df
    return loaded


# ---------------------------------------------------------------------------
# Summary helpers
# ---------------------------------------------------------------------------

def summary_rows(data: dict[str, pd.DataFrame]) -> list[dict]:
    """Build summary rows for a set of loaded BED DataFrames."""
    rows = []
    for name, df in data.items():
        total_bp = int(df["width"].sum())
        rows.append({
            "Interval": name,
            "Count": len(df),
            "Total bp": total_bp,
            "Total bp (fmt)": format_bp(total_bp),
            "Mean width": int(df["width"].mean()),
            "Median width": int(df["width"].median()),
            "% GRCh38": 100 * total_bp / GRCH38_TOTAL_BP,
        })
    return rows


def print_summary_table(data: dict[str, pd.DataFrame]) -> None:
    """Print a formatted summary table to stdout."""
    rows = []
    for name, df in data.items():
        total_bp = df["width"].sum()
        rows.append({
            "Interval": name,
            "Count": f"{len(df):>10,}",
            "Total bp": f"{format_bp(total_bp):>12}",
            "Mean width": f"{df['width'].mean():>10,.0f}",
            "Median width": f"{df['width'].median():>10,.0f}",
            "% GRCh38": f"{100 * total_bp / GRCH38_TOTAL_BP:>8.2f}%",
        })
    summary = pd.DataFrame(rows)
    print("\n" + summary.to_string(index=False) + "\n")


def summary_markdown_table(data: dict[str, pd.DataFrame]) -> str:
    """Return a Markdown table of summary stats."""
    rows = summary_rows(data)
    lines = [
        "| Interval | Count | Total bp | Mean width | Median width | % GRCh38 |",
        "|----------|------:|----------|------------|--------------|----------|",
    ]
    for r in rows:
        lines.append(
            f"| {r['Interval']} | {r['Count']:,} | {r['Total bp (fmt)']} "
            f"| {r['Mean width']:,} | {r['Median width']:,} | {r['% GRCh38']:.2f}% |"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plot functions (return Figure, caller decides what to do with it)
# ---------------------------------------------------------------------------

def plot_coverage_bars(data: dict[str, pd.DataFrame]) -> plt.Figure:
    """Horizontal bar chart of total bp per BED file (log scale)."""
    names = list(data.keys())
    totals = [data[n]["width"].sum() for n in names]
    order = sorted(range(len(names)), key=lambda i: totals[i])
    names = [names[i] for i in order]
    totals = [totals[i] for i in order]

    fig, ax = plt.subplots(figsize=(10, max(3, len(names) * 0.7)))
    colors = sns.color_palette("viridis", len(names))
    bars = ax.barh(names, totals, color=colors)

    ax.set_xscale("log")
    ax.set_xlabel("Total base pairs (log scale)")
    ax.set_title("Genomic Coverage by Interval Set")

    for bar, bp in zip(bars, totals):
        pct = 100 * bp / GRCH38_TOTAL_BP
        label = f"  {format_bp(bp)} ({pct:.2f}%)"
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                label, va="center", fontsize=9)

    fig.tight_layout()
    return fig


def plot_width_distributions(data: dict[str, pd.DataFrame]) -> plt.Figure:
    """Violin + box plots of interval width distributions (log scale)."""
    records = []
    for name, df in data.items():
        records.append(pd.DataFrame({"Interval": name, "width": df["width"]}))
    combined = pd.concat(records, ignore_index=True)

    medians = combined.groupby("Interval")["width"].median().sort_values()
    order = medians.index.tolist()

    fig, ax = plt.subplots(figsize=(10, max(3, len(order) * 0.8)))
    sns.violinplot(data=combined, y="Interval", x="width", order=order,
                   inner="box", density_norm="width", cut=0, ax=ax, orient="h")
    ax.set_xscale("log")
    ax.set_xlabel("Interval width (bp, log scale)")
    ax.set_title("Interval Width Distributions")

    for i, name in enumerate(order):
        med = medians[name]
        ax.text(med, i, f"  med={med:,.0f}", va="center", fontsize=8, alpha=0.7)

    fig.tight_layout()
    return fig


def plot_chromosome_coverage(data: dict[str, pd.DataFrame]) -> plt.Figure:
    """Grouped bar chart of per-chromosome bp coverage."""
    chroms = sorted(KARYOTYPE_ORDER.keys(), key=lambda c: KARYOTYPE_ORDER[c])
    names = list(data.keys())

    per_chrom = {}
    for name, df in data.items():
        bp_by_chr = df.groupby("chrom")["width"].sum()
        per_chrom[name] = [bp_by_chr.get(c, 0) for c in chroms]

    x = np.arange(len(chroms))
    width = 0.8 / len(names)
    colors = sns.color_palette("tab10", len(names))

    fig, ax = plt.subplots(figsize=(16, 6))
    for i, name in enumerate(names):
        offset = (i - len(names) / 2 + 0.5) * width
        ax.bar(x + offset, per_chrom[name], width, label=name, color=colors[i])

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in chroms], fontsize=9)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Total base pairs")
    ax.set_yscale("log")
    ax.set_title("Per-Chromosome Coverage by Interval Set")
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_probe_report(data: dict[str, pd.DataFrame], interval_dir: str,
                          dpi: int = DEFAULT_DPI) -> str:
    """Generate a self-contained Markdown report for the probe subcommand."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# BED Interval Analysis",
        "",
        f"**Generated:** {now}  |  **Directory:** `{interval_dir}`",
        "",
        "## Summary",
        "",
        summary_markdown_table(data),
        "",
        "## Genomic Coverage",
        "",
        fig_to_markdown(plot_coverage_bars(data), alt="Coverage bar chart", dpi=dpi),
        "",
        "## Interval Width Distributions",
        "",
        fig_to_markdown(plot_width_distributions(data), alt="Width distributions", dpi=dpi),
        "",
        "## Per-Chromosome Coverage",
        "",
        fig_to_markdown(plot_chromosome_coverage(data), alt="Chromosome coverage", dpi=dpi),
        "",
    ]

    # Recommendation: identify outlier by total bp
    rows = summary_rows(data)
    if len(rows) >= 2:
        rows_sorted = sorted(rows, key=lambda r: r["Total bp"], reverse=True)
        biggest = rows_sorted[0]
        second = rows_sorted[1]
        ratio = biggest["Total bp"] / second["Total bp"] if second["Total bp"] > 0 else float("inf")
        if ratio > 5:
            lines.extend([
                "## Recommendation",
                "",
                f"**`{biggest['Interval']}`** is {ratio:.0f}x larger than the next-biggest interval "
                f"set (`{second['Interval']}`). This extreme spread can defeat Hail partition pruning "
                "and create class imbalance in batch-effect statistics.",
                "",
                "Consider downsampling to match the second-largest set:",
                "",
                "```bash",
                "python downsample_bed.py sample \\",
                f"    --input <{biggest['Interval']}_path> \\",
                f"    --match-file <{second['Interval']}_path> \\",
                f"    --output <{biggest['Interval']}_downsampled.bed>",
                "```",
                "",
            ])

    return "\n".join(lines)


def generate_sample_report(
    input_name: str,
    input_df: pd.DataFrame,
    output_df: pd.DataFrame,
    target_bp: int,
    target_desc: str,
    seed: int,
    stratified: bool,
    context_data: dict[str, pd.DataFrame],
    output_path: str,
    dpi: int = DEFAULT_DPI,
) -> str:
    """Generate a self-contained Markdown report for the sample subcommand."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    input_bp = int(input_df["width"].sum())
    output_bp = int(output_df["width"].sum())
    out_name = Path(output_path).stem

    lines = [
        "# BED Downsampling Report",
        "",
        f"**Generated:** {now}",
        "",
        f"| | Value |",
        f"|---|---|",
        f"| **Input** | `{input_name}` — {len(input_df):,} intervals, {format_bp(input_bp)} |",
        f"| **Target** | {format_bp(target_bp)} ({target_desc}) |",
        f"| **Output** | `{Path(output_path).name}` — {len(output_df):,} intervals, {format_bp(output_bp)} |",
        f"| **Reduction** | {100 * (1 - output_bp / input_bp):.1f}% of original bp |",
        "",
        "## Before / After",
        "",
        "| Metric | Before | After |",
        "|--------|--------|-------|",
        f"| Intervals | {len(input_df):,} | {len(output_df):,} |",
        f"| Total bp | {format_bp(input_bp)} | {format_bp(output_bp)} |",
        f"| Mean width | {int(input_df['width'].mean()):,} | {int(output_df['width'].mean()):,} |",
        f"| Median width | {int(input_df['width'].median()):,} | {int(output_df['width'].median()):,} |",
        f"| % GRCh38 | {100 * input_bp / GRCH38_TOTAL_BP:.2f}% | {100 * output_bp / GRCH38_TOTAL_BP:.2f}% |",
        "",
    ]

    # Build combined data for plots: context + before/after
    plot_data = dict(context_data)
    plot_data[f"{input_name} (before)"] = input_df
    plot_data[out_name] = output_df

    lines.extend([
        "## All Intervals — Coverage",
        "",
        fig_to_markdown(plot_coverage_bars(plot_data), alt="Coverage comparison", dpi=dpi),
        "",
        "## Width Distributions",
        "",
        fig_to_markdown(plot_width_distributions(plot_data), alt="Width distributions", dpi=dpi),
        "",
        "## Per-Chromosome Coverage",
        "",
        fig_to_markdown(plot_chromosome_coverage(plot_data), alt="Chromosome coverage", dpi=dpi),
        "",
        "## Parameters",
        "",
        f"| Parameter | Value |",
        f"|-----------|-------|",
        f"| Seed | {seed} |",
        f"| Stratified | {stratified} |",
        f"| Target source | {target_desc} |",
        "",
        "## Provenance",
        "",
        f"See `{Path(output_path).name}.manifest.json` for full machine-readable provenance.",
        "",
    ])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Downsampling core (unchanged logic from original)
# ---------------------------------------------------------------------------

def compute_target_bp(args) -> tuple[int, str]:
    """Resolve the target bp from whichever CLI mode was used.

    Returns (target_bp, description_string).
    """
    if args.target_bp is not None:
        return args.target_bp, f"explicit --target-bp {args.target_bp}"

    if args.match_file is not None:
        ref = load_bed(args.match_file)
        bp = int(ref["width"].sum())
        return bp, f"--match-file {Path(args.match_file).name} ({format_bp(bp)})"

    # --match-files mode
    totals = {}
    for path in args.match_files:
        ref = load_bed(path)
        totals[path] = int(ref["width"].sum())

    stat = args.match_stat
    if stat == "max":
        bp = max(totals.values())
    elif stat == "min":
        bp = min(totals.values())
    elif stat == "median":
        bp = int(np.median(list(totals.values())))
    elif stat == "mean":
        bp = int(np.mean(list(totals.values())))
    else:
        raise ValueError(f"Unknown --match-stat: {stat}")

    files_str = ", ".join(f"{Path(p).name}={format_bp(v)}" for p, v in totals.items())
    return bp, f"--match-files {stat}({files_str}) = {format_bp(bp)}"


def downsample_stratified(df: pd.DataFrame, target_bp: int,
                          rng: np.random.Generator) -> pd.DataFrame:
    """Downsample with proportional allocation per chromosome."""
    total_bp = df["width"].sum()
    if target_bp >= total_bp:
        return df

    chrom_groups = df.groupby("chrom")
    chrom_bp = chrom_groups["width"].sum()

    chrom_targets = {}
    for chrom, bp in chrom_bp.items():
        chrom_targets[chrom] = max(1, int(target_bp * bp / total_bp))

    kept = []
    for chrom, group_df in chrom_groups:
        budget = chrom_targets[chrom]
        shuffled = group_df.sample(frac=1, random_state=rng.integers(2**31)).reset_index(drop=True)
        cumulative = shuffled["width"].cumsum()
        mask = cumulative <= budget
        if not mask.any():
            kept.append(shuffled.nsmallest(1, "width"))
        else:
            kept.append(shuffled[mask])

    return pd.concat(kept, ignore_index=True)


def downsample_uniform(df: pd.DataFrame, target_bp: int,
                       rng: np.random.Generator) -> pd.DataFrame:
    """Downsample uniformly without chromosome stratification."""
    total_bp = df["width"].sum()
    if target_bp >= total_bp:
        return df

    shuffled = df.sample(frac=1, random_state=rng.integers(2**31)).reset_index(drop=True)
    cumulative = shuffled["width"].cumsum()
    mask = cumulative <= target_bp
    if not mask.any():
        return shuffled.head(1)
    return shuffled[mask]


def sort_karyotype(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by karyotype order, then by start position."""
    return df.sort_values(["chrom_order", "start"]).reset_index(drop=True)


def write_bed(df: pd.DataFrame, path: str) -> None:
    """Write a 3-column BED file."""
    df[["chrom", "start", "end"]].to_csv(path, sep="\t", header=False, index=False)


def write_manifest(
    output_path: str,
    input_path: str,
    input_df: pd.DataFrame,
    output_df: pd.DataFrame,
    target_bp: int,
    target_desc: str,
    seed: int,
    stratified: bool,
) -> str:
    """Write a JSON manifest with full provenance. Returns the manifest path."""
    input_chrom = input_df.groupby("chrom")["width"].agg(["count", "sum"]).rename(
        columns={"count": "intervals", "sum": "total_bp"}
    )
    output_chrom = output_df.groupby("chrom")["width"].agg(["count", "sum"]).rename(
        columns={"count": "intervals", "sum": "total_bp"}
    )

    chrom_breakdown = {}
    for chrom in sorted(set(input_chrom.index) | set(output_chrom.index),
                        key=lambda c: KARYOTYPE_ORDER.get(c, 99)):
        entry = {}
        if chrom in input_chrom.index:
            entry["input_intervals"] = int(input_chrom.loc[chrom, "intervals"])
            entry["input_bp"] = int(input_chrom.loc[chrom, "total_bp"])
        if chrom in output_chrom.index:
            entry["output_intervals"] = int(output_chrom.loc[chrom, "intervals"])
            entry["output_bp"] = int(output_chrom.loc[chrom, "total_bp"])
        chrom_breakdown[chrom] = entry

    manifest = {
        "created": datetime.now(timezone.utc).isoformat(),
        "command": f"python {' '.join(sys.argv)}",
        "input": {
            "path": str(Path(input_path).resolve()),
            "intervals": len(input_df),
            "total_bp": int(input_df["width"].sum()),
            "total_bp_fmt": format_bp(int(input_df["width"].sum())),
        },
        "output": {
            "path": str(Path(output_path).resolve()),
            "intervals": len(output_df),
            "total_bp": int(output_df["width"].sum()),
            "total_bp_fmt": format_bp(int(output_df["width"].sum())),
        },
        "target": {
            "bp": target_bp,
            "bp_fmt": format_bp(target_bp),
            "source": target_desc,
        },
        "params": {
            "seed": seed,
            "stratified": stratified,
        },
        "per_chromosome": chrom_breakdown,
    }

    manifest_path = output_path + ".manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Manifest: {manifest_path}")
    return manifest_path


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def cmd_probe(args) -> None:
    """Handle the 'probe' subcommand: analyze BED files and generate report."""
    dpi = args.dpi

    # Resolve BED files
    if args.interval_dir:
        beds = discover_beds(args.interval_dir)
        interval_dir = args.interval_dir
    else:
        beds = parse_bed_specs(args.beds)
        interval_dir = "."

    if not beds:
        print("No BED files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {len(beds)} BED files...")
    data = load_all(beds)

    print_summary_table(data)

    # Report path
    if args.report:
        report_path = args.report
    else:
        report_path = os.path.join(interval_dir, "probe_report.md")

    os.makedirs(os.path.dirname(os.path.abspath(report_path)), exist_ok=True)

    print(f"Generating report: {report_path}")
    report = generate_probe_report(data, interval_dir, dpi=dpi)
    with open(report_path, "w") as f:
        f.write(report)

    print(f"  Report: {report_path}")
    print("Done.")


def cmd_sample(args) -> None:
    """Handle the 'sample' subcommand: downsample + document."""
    dpi = args.dpi

    # Load input
    print(f"Loading input: {args.input}")
    input_df = load_bed(args.input)
    input_bp = int(input_df["width"].sum())
    input_name = Path(args.input).stem
    print(f"  {len(input_df):,} intervals, {format_bp(input_bp)}")

    # Resolve target
    target_bp, target_desc = compute_target_bp(args)
    print(f"Target: {format_bp(target_bp)} ({target_desc})")

    if target_bp >= input_bp:
        print(f"  WARNING: target ({format_bp(target_bp)}) >= input ({format_bp(input_bp)})")
        print("  Returning full file unchanged.")
        output_df = input_df
    else:
        rng = np.random.default_rng(args.seed)
        if args.stratify:
            print("  Downsampling (stratified by chromosome)...")
            output_df = downsample_stratified(input_df, target_bp, rng)
        else:
            print("  Downsampling (uniform)...")
            output_df = downsample_uniform(input_df, target_bp, rng)

    output_df = sort_karyotype(output_df)

    parent = Path(args.output).parent
    if parent != Path("."):
        os.makedirs(str(parent), exist_ok=True)

    write_bed(output_df, args.output)
    achieved_bp = int(output_df["width"].sum())
    print(f"  Output: {len(output_df):,} intervals, {format_bp(achieved_bp)}")
    pct_reduction = 100 * (1 - achieved_bp / input_bp)
    print(f"  Reduction: {pct_reduction:.1f}% of original bp coverage")

    # Manifest
    write_manifest(args.output, args.input, input_df, output_df,
                   target_bp, target_desc, args.seed, args.stratify)

    # Report (unless --no-report)
    if not args.no_report:
        # Discover sibling BED files for comparison context
        if args.interval_dir:
            sibling_beds = discover_beds(args.interval_dir)
        elif args.beds:
            sibling_beds = parse_bed_specs(args.beds)
        else:
            # Auto-discover from input's parent directory
            sibling_beds = discover_beds(str(Path(args.input).parent))

        # Remove the input file itself from siblings (we show it as before/after)
        sibling_beds = {
            name: path for name, path in sibling_beds.items()
            if os.path.abspath(path) != os.path.abspath(args.input)
        }

        print(f"Loading {len(sibling_beds)} sibling BED files for context...")
        context_data = load_all(sibling_beds) if sibling_beds else {}

        report_path = args.output + ".report.md"
        print(f"Generating report: {report_path}")
        report = generate_sample_report(
            input_name=input_name,
            input_df=input_df,
            output_df=output_df,
            target_bp=target_bp,
            target_desc=target_desc,
            seed=args.seed,
            stratified=args.stratify,
            context_data=context_data,
            output_path=args.output,
            dpi=dpi,
        )
        with open(report_path, "w") as f:
            f.write(report)
        print(f"  Report: {report_path}")

    print("Done.")


# ---------------------------------------------------------------------------
# CLI parser
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze and downsample BED interval files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python downsample_bed.py probe --interval-dir intervals/\n"
            "  python downsample_bed.py sample --input big.bed --match-file ref.bed --output out.bed\n"
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- probe ---
    probe_parser = subparsers.add_parser(
        "probe", help="Analyze BED files and generate a summary report",
    )
    probe_input = probe_parser.add_mutually_exclusive_group(required=True)
    probe_input.add_argument("--interval-dir", help="Directory to auto-discover BED files")
    probe_input.add_argument("--beds", nargs="+", metavar="NAME=PATH",
                             help="Explicit NAME=path pairs")
    probe_parser.add_argument("--report", help="Output report path (default: <interval-dir>/probe_report.md)")
    probe_parser.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Figure DPI (default: 150)")

    # --- sample ---
    sample_parser = subparsers.add_parser(
        "sample", help="Downsample a BED file and generate a before/after report",
    )
    sample_parser.add_argument("--input", required=True, help="Input BED file")
    sample_parser.add_argument("--output", required=True, help="Output BED file path")
    sample_parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

    # Stratification
    strat_group = sample_parser.add_mutually_exclusive_group()
    strat_group.add_argument("--stratify-by-chrom", dest="stratify", action="store_true",
                             default=True, help="Proportional sampling per chromosome (default)")
    strat_group.add_argument("--no-stratify", dest="stratify", action="store_false",
                             help="Uniform random sampling ignoring chromosomes")

    # Target specification
    target_group = sample_parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument("--target-bp", type=int, help="Explicit target total bp")
    target_group.add_argument("--match-file", help="Match one reference BED file's coverage")
    target_group.add_argument("--match-files", nargs="+",
                              help="Match a stat across multiple reference BED files")
    sample_parser.add_argument("--match-stat", default="max",
                               choices=["max", "min", "median", "mean"],
                               help="Statistic for --match-files (default: max)")

    # Report control
    sample_parser.add_argument("--no-report", action="store_true",
                               help="Skip generating the Markdown report")
    sample_parser.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                               help="Figure DPI (default: 150)")

    # Context BED files for report
    context_group = sample_parser.add_mutually_exclusive_group()
    context_group.add_argument("--interval-dir",
                               help="Directory for sibling BED files (default: input's parent)")
    context_group.add_argument("--beds", nargs="+", metavar="NAME=PATH",
                               help="Explicit sibling BED files for report context")

    args = parser.parse_args()

    if args.command == "probe":
        cmd_probe(args)
    elif args.command == "sample":
        cmd_sample(args)


if __name__ == "__main__":
    main()
