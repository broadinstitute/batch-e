#!/usr/bin/env python3
"""
batch_e_reporter.py — Self-contained HTML report generator for batch effect results.

Reads pairwise_comparisons.tsv, group_summaries.tsv, and optionally sample_stats.tsv,
config.json, and timing_metrics.json from a results directory (local or GCS). Generates
a single HTML file with embedded matplotlib/seaborn charts as base64 PNGs.

Usage:
    python batch_e_reporter.py <input_path> [-o output.html] [--title "..."] \
        [--no-sample-stats] [--effect-threshold 0.5]
"""

import base64
import io
import json
import logging
import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

APPRECIABLE_THRESHOLD = 0.5  # |d| >= 0.5 = "appreciable batch effect"

KNOWN_METRIC_SUFFIXES = [
    "_snp_ti_het", "_snp_tv_het", "_snp_ti_hom", "_snp_tv_hom",
    "_snp_het", "_snp_hom", "_snp_total",
    "_titv_het", "_titv_hom", "_titv_total",
    "_indel_ins_het", "_indel_del_het", "_indel_ins_hom", "_indel_del_hom",
    "_indel_het", "_indel_hom", "_indel_total",
    "_delins_het", "_delins_hom", "_delins_total",
]

KEY_METRIC_TYPES = ["snp_total", "indel_total", "titv_total", "delins_total"]

METRIC_DISPLAY_NAMES = {
    "snp_total": "SNP Count",
    "indel_total": "Indel Count",
    "titv_total": "Ti/Tv Ratio",
    "delins_total": "Del/Ins Ratio",
}

SEVERITY_COLORS = {
    "not appreciable": "#d4edda",  # green
    "appreciable": "#f8d7da",      # red/pink
}

REPORT_CSS = """
body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    max-width: 1200px; margin: 0 auto; padding: 20px;
    color: #333; background: #fafafa; line-height: 1.5;
}
h1 { color: #1a237e; border-bottom: 3px solid #1a237e; padding-bottom: 8px; }
h2 { color: #283593; margin-top: 40px; border-bottom: 1px solid #ccc; padding-bottom: 4px; }
h3 { color: #3949ab; }
.card {
    background: white; border: 1px solid #ddd; border-radius: 8px;
    padding: 16px; margin: 12px 0; box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}
.summary-stat { display: inline-block; margin: 8px 16px 8px 0; }
.summary-stat .label { font-size: 0.85em; color: #666; }
.summary-stat .value { font-size: 1.4em; font-weight: 600; }
table {
    border-collapse: collapse; width: 100%; margin: 12px 0;
    font-size: 0.9em;
}
th, td { border: 1px solid #ddd; padding: 6px 10px; text-align: left; }
th { background: #f5f5f5; font-weight: 600; position: sticky; top: 0; }
tr:hover { background: #f9f9f9; }
.fig-container { text-align: center; margin: 16px 0; }
.fig-container img { max-width: 100%; height: auto; border: 1px solid #eee; border-radius: 4px; }
.caption { font-size: 0.85em; color: #666; margin-top: 4px; font-style: italic; }
.alert { padding: 12px 16px; border-radius: 6px; margin: 12px 0; }
.alert-info { background: #e3f2fd; border-left: 4px solid #2196f3; }
.alert-warn { background: #fff3e0; border-left: 4px solid #ff9800; }
.alert-error { background: #ffebee; border-left: 4px solid #f44336; }
details { margin: 8px 0; }
summary { cursor: pointer; font-weight: 600; color: #1565c0; }
.methods { font-size: 0.9em; color: #555; }
.grid-2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
@media (max-width: 768px) { .grid-2 { grid-template-columns: 1fr; } }
@media print {
    @page { margin: 1cm 1.2cm; }
    body { background: white; max-width: none; padding: 0; margin: 0; }
    .card { box-shadow: none; border: 1px solid #ccc; }
    h2 { break-after: avoid; }
    .fig-container { break-inside: avoid; }
    .fig-container img { max-width: 95%; }
    table { font-size: 0.8em; }
    th { position: static !important; background: #f5f5f5 !important;
         -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    tr[style] { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .summary-stat { display: inline-block !important; }
    details { display: block; }
    details > summary { display: none; }
    details > *:not(summary) { display: block; }
}
"""

# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------

@dataclass
class ReportData:
    """Container for all loaded report data."""
    comparisons: pd.DataFrame = field(default_factory=pd.DataFrame)
    group_summaries: pd.DataFrame = field(default_factory=pd.DataFrame)
    sample_stats: Optional[pd.DataFrame] = None
    config: Optional[Dict[str, Any]] = None
    timing: Optional[Dict[str, Any]] = None
    input_path: str = ""
    effect_threshold: float = 0.5
    comparison_name: str = "group"  # human label for comparison variable


def _open_file(path: str, mode: str = "r"):
    """Open a file from GCS or local filesystem."""
    if path.startswith("gs://"):
        import gcsfs
        fs = gcsfs.GCSFileSystem(requester_pays=True)
        return fs.open(path.replace("gs://", ""), mode)
    return open(path, mode)


def _file_exists(path: str) -> bool:
    """Check if a file exists on GCS or locally."""
    if path.startswith("gs://"):
        try:
            import gcsfs
            fs = gcsfs.GCSFileSystem(requester_pays=True)
            return fs.exists(path.replace("gs://", ""))
        except Exception:
            return False
    return os.path.isfile(path)


def _join_path(base: str, name: str) -> str:
    """Join a base path with a filename, handling GCS and local paths."""
    if base.startswith("gs://"):
        return base.rstrip("/") + "/" + name
    return os.path.join(base, name)


def _validate_comparisons(df: pd.DataFrame) -> pd.DataFrame:
    """Validate that comparisons TSV has required v2 columns."""
    required = {"group_x", "group_y", "comparison"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Comparisons TSV missing required columns: {missing}. "
            f"Expected: group_x, group_y, comparison"
        )
    return df


def load_report_data(
    input_path: str,
    load_sample_stats: bool = True,
    effect_threshold: float = 0.5,
) -> ReportData:
    """Load all available result files from a directory.

    Args:
        input_path: Local directory or GCS path containing result files.
        load_sample_stats: Whether to load the (potentially large) sample_stats.tsv.
        effect_threshold: Cohen's d threshold for flagging effects.

    Returns:
        ReportData with all loaded DataFrames and metadata.
    """
    data = ReportData(input_path=input_path, effect_threshold=effect_threshold)

    # Required files
    comp_path = _join_path(input_path, "pairwise_comparisons.tsv")
    summ_path = _join_path(input_path, "group_summaries.tsv")

    logger.info("Loading pairwise_comparisons.tsv ...")
    with _open_file(comp_path) as f:
        data.comparisons = pd.read_csv(f, sep="\t")
    logger.info(f"  {len(data.comparisons)} comparison rows loaded")

    # Validate schema
    data.comparisons = _validate_comparisons(data.comparisons)

    logger.info("Loading group_summaries.tsv ...")
    with _open_file(summ_path) as f:
        data.group_summaries = pd.read_csv(f, sep="\t")
    logger.info(f"  {len(data.group_summaries)} summary rows loaded")

    # Optional files
    config_path = _join_path(input_path, "config.json")
    if _file_exists(config_path):
        logger.info("Loading config.json ...")
        with _open_file(config_path) as f:
            data.config = json.load(f)

    timing_path = _join_path(input_path, "timing_metrics.json")
    if _file_exists(timing_path):
        logger.info("Loading timing_metrics.json ...")
        with _open_file(timing_path) as f:
            data.timing = json.load(f)

    # Determine comparison_name from config or from data
    if data.config and data.config.get("comparison_name"):
        data.comparison_name = data.config["comparison_name"]
    elif "comparison" in data.comparisons.columns:
        # Use the first value from the comparison column
        data.comparison_name = data.comparisons["comparison"].iloc[0]
    else:
        data.comparison_name = "group"

    if load_sample_stats:
        stats_path = _join_path(input_path, "sample_stats.tsv")
        if _file_exists(stats_path):
            logger.info("Loading sample_stats.tsv ...")
            with _open_file(stats_path) as f:
                data.sample_stats = pd.read_csv(f, sep="\t")
            logger.info(f"  {len(data.sample_stats)} sample rows loaded")

    return data


# ---------------------------------------------------------------------------
# Metric Parsing Helpers
# ---------------------------------------------------------------------------

def parse_metric_name(metric: str) -> Tuple[str, str]:
    """Split a compound metric name into (interval, metric_type).

    Handles underscore-heavy names like 'GC_gt_85_indel_total' by matching
    known suffixes from longest to shortest.

    Returns:
        Tuple of (interval_name, metric_type), e.g. ('Low_Mappability', 'snp_total').
        Falls back to ('unknown', metric) if no known suffix matches.
    """
    for suffix in sorted(KNOWN_METRIC_SUFFIXES, key=len, reverse=True):
        if metric.endswith(suffix):
            interval = metric[: -len(suffix)]
            metric_type = suffix.lstrip("_")
            if interval:
                return interval, metric_type
    return "unknown", metric


def severity_label(d: float) -> str:
    """Return a severity label for a Cohen's d value."""
    return "appreciable" if abs(d) >= APPRECIABLE_THRESHOLD else "not appreciable"


def display_metric(metric_type: str) -> str:
    """Return a human-readable display name for a metric type."""
    return METRIC_DISPLAY_NAMES.get(metric_type, metric_type)


def severity_color(d: float) -> str:
    """Return a background color for a Cohen's d value."""
    return SEVERITY_COLORS[severity_label(d)]


# ---------------------------------------------------------------------------
# Figure Utilities
# ---------------------------------------------------------------------------

def fig_to_base64(fig: plt.Figure, dpi: int = 130) -> str:
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def embed_figure(fig: plt.Figure, caption: str = "", dpi: int = 130) -> str:
    """Render a matplotlib figure as an embedded HTML img tag."""
    b64 = fig_to_base64(fig, dpi=dpi)
    html = f'<div class="fig-container">'
    html += f'<img src="data:image/png;base64,{b64}" alt="{caption}">'
    if caption:
        html += f'<div class="caption">{caption}</div>'
    html += "</div>"
    return html


def safe_render(fn, *args, **kwargs) -> str:
    """Call a chart function, catching exceptions and returning fallback HTML."""
    try:
        return fn(*args, **kwargs)
    except Exception as e:
        name = getattr(fn, "__name__", "chart")
        logger.warning(f"Failed to render {name}: {e}", exc_info=True)
        return (
            f'<div class="alert alert-error">'
            f"Could not render {name}: {e}</div>"
        )


def _get_group_colors(groups: List[str]) -> Dict[str, str]:
    """Assign consistent bold colors to groups."""
    BOLD_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                   "#a65628", "#f781bf", "#999999"]
    sorted_groups = sorted(groups)
    return {g: BOLD_COLORS[i % len(BOLD_COLORS)] for i, g in enumerate(sorted_groups)}


# ---------------------------------------------------------------------------
# Chart Functions
# ---------------------------------------------------------------------------

def chart_significance_table(data: ReportData) -> str:
    """Render the executive summary significance table."""
    df = data.comparisons.copy()
    threshold = data.effect_threshold

    flagged = df[
        (df["cohens_d"].abs() > threshold) | (df["p_value"] < 0.05)
    ].copy()
    flagged["abs_d"] = flagged["cohens_d"].abs()
    flagged = flagged.sort_values("abs_d", ascending=False)

    # Parse metric into interval + metric_type
    parsed = flagged["metric"].apply(parse_metric_name)
    flagged["interval"] = [p[0] for p in parsed]
    flagged["metric_type"] = [p[1] for p in parsed]
    flagged["severity"] = flagged["cohens_d"].apply(severity_label)
    flagged["group_pair"] = flagged["group_x"] + " vs " + flagged["group_y"]

    total = len(df)
    n_flagged = len(flagged)
    n_appreciable = (flagged["abs_d"] >= APPRECIABLE_THRESHOLD).sum()

    # Text summary
    html = '<div class="card">'
    html += f"<p><strong>{n_appreciable}</strong> of {total} comparisons show "
    html += f"appreciable effects (|d| &ge; {APPRECIABLE_THRESHOLD})."

    if n_flagged > 0:
        top_intervals = flagged["interval"].value_counts().head(3).index.tolist()
        html += f" Most affected intervals: {', '.join(top_intervals)}."

    html += "</p></div>"

    if n_flagged == 0:
        html += '<div class="alert alert-info">No comparisons exceed the effect size threshold.</div>'
        return html

    # Build table
    use_details = n_flagged > 20
    table_html = _build_flagged_table(flagged, data.comparison_name)

    if use_details:
        top_20 = _build_flagged_table(flagged.head(20), data.comparison_name)
        html += top_20
        html += f"<details><summary>Show all {n_flagged} flagged comparisons</summary>"
        html += table_html
        html += "</details>"
    else:
        html += table_html

    return html


def _build_flagged_table(df: pd.DataFrame, comparison_name: str) -> str:
    """Build an HTML table from flagged comparisons."""
    ancestry_col = "ancestry" if "ancestry" in df.columns else None
    pair_label = f"{comparison_name.replace('_', ' ').title()} Pair"

    html = "<table><thead><tr>"
    if ancestry_col:
        html += "<th>Ancestry</th>"
    html += f"<th>{pair_label}</th><th>Interval</th><th>Metric</th>"
    html += "<th>Cohen's d</th><th>p-value</th><th>Severity</th>"
    html += "</tr></thead><tbody>"

    for _, row in df.iterrows():
        bg = severity_color(row["cohens_d"])
        html += f'<tr style="background:{bg}">'
        if ancestry_col:
            anc = row.get("ancestry", "")
            html += f"<td>{anc}</td>"
        html += f"<td>{row['group_pair']}</td>"
        html += f"<td>{row['interval']}</td>"
        html += f"<td>{display_metric(row['metric_type'])}</td>"
        html += f"<td>{row['cohens_d']:.3f}</td>"
        html += f"<td>{row['p_value']:.2e}</td>"
        html += f"<td>{row['severity']}</td>"
        html += "</tr>"

    html += "</tbody></table>"
    return html


def chart_effect_heatmap(data: ReportData) -> str:
    """Render Cohen's d heatmaps, one per ancestry."""
    df = data.comparisons.copy()
    parsed = df["metric"].apply(parse_metric_name)
    df["interval"] = [p[0] for p in parsed]
    df["metric_type"] = [p[1] for p in parsed]
    df["group_pair"] = df["group_x"] + " vs " + df["group_y"]
    df["row_label"] = df["interval"] + " / " + df["metric_type"].map(display_metric)

    # Filter to key metrics for readability
    df = df[df["metric_type"].isin(KEY_METRIC_TYPES)]

    ancestries = sorted(df["ancestry"].dropna().unique()) if "ancestry" in df.columns else [None]
    html = ""

    for anc in ancestries:
        sub = df[df["ancestry"] == anc] if anc else df
        if sub.empty:
            continue

        pivot = sub.pivot_table(
            index="row_label", columns="group_pair",
            values="cohens_d", aggfunc="first",
        )
        if pivot.empty:
            continue

        # Annotation: d value + asterisk if significant
        pval_pivot = sub.pivot_table(
            index="row_label", columns="group_pair",
            values="p_value", aggfunc="first",
        )

        annot = pivot.copy().astype(str)
        for ri in pivot.index:
            for ci in pivot.columns:
                d_val = pivot.loc[ri, ci]
                p_val = pval_pivot.loc[ri, ci] if ri in pval_pivot.index and ci in pval_pivot.columns else 1.0
                if pd.isna(d_val):
                    annot.loc[ri, ci] = ""
                else:
                    star = "*" if (not pd.isna(p_val) and p_val < 0.05) else ""
                    annot.loc[ri, ci] = f"{d_val:.2f}{star}"

        n_rows = max(len(pivot), 1)
        n_cols = max(len(pivot.columns), 1)
        fig_h = max(4, n_rows * 0.5 + 1.5)
        fig_w = max(6, n_cols * 1.8 + 3)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        sns.heatmap(
            pivot.astype(float), annot=annot, fmt="",
            cmap="RdBu_r", norm=norm, ax=ax,
            linewidths=0.5, linecolor="white",
            cbar_kws={"label": "Cohen's d", "shrink": 0.8},
        )

        title = "Effect Size Heatmap"
        if anc:
            title += f" — {anc.upper()}"
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_ylabel("")
        ax.set_xlabel("")
        plt.xticks(rotation=30, ha="right", fontsize=9)
        plt.yticks(fontsize=9)

        caption = "Cohen's d values. * indicates p < 0.05. Color clamped to [-1, 1]."
        html += embed_figure(fig, caption)

    return html


def chart_volcano(data: ReportData) -> str:
    """Render volcano plots: Cohen's d vs -log10(p), one panel per ancestry."""
    df = data.comparisons.copy()
    parsed = df["metric"].apply(parse_metric_name)
    df["interval"] = [p[0] for p in parsed]
    df["metric_type"] = [p[1] for p in parsed]

    ancestries = sorted(df["ancestry"].dropna().unique()) if "ancestry" in df.columns else [None]
    n_anc = len(ancestries)
    fig, axes = plt.subplots(1, n_anc, figsize=(7 * n_anc, 6), squeeze=False)

    intervals = sorted(df["interval"].unique())
    BOLD_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                   "#a65628", "#f781bf", "#999999"]
    interval_palette = dict(zip(intervals, BOLD_COLORS[:len(intervals)]))

    # Bonferroni threshold
    n_tests = len(df)
    bonf_thresh = -np.log10(0.05 / n_tests) if n_tests > 0 else 5

    for i, anc in enumerate(ancestries):
        ax = axes[0, i]
        sub = df[df["ancestry"] == anc] if anc else df

        for interval in intervals:
            isub = sub[sub["interval"] == interval]
            if isub.empty:
                continue
            ax.scatter(
                isub["cohens_d"], isub["neg_log10_p"],
                c=interval_palette[interval],
                s=60, alpha=0.8, label=interval,
                edgecolors="black", linewidth=0.4, zorder=3,
            )

        # Label outlier points (|d| > 0.3 and above Bonferroni)
        outliers = sub[(sub["cohens_d"].abs() > 0.3) & (sub["neg_log10_p"] > bonf_thresh)]
        for _, row in outliers.iterrows():
            label = f"{row['interval']}\n{row['metric_type']}"
            ax.annotate(
                label, (row["cohens_d"], row["neg_log10_p"]),
                fontsize=7, ha="center", va="bottom",
                xytext=(0, 6), textcoords="offset points",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="grey", alpha=0.8, lw=0.5),
            )

        ax.axvline(data.effect_threshold, ls="--", c="grey", lw=1, alpha=0.5)
        ax.axvline(-data.effect_threshold, ls="--", c="grey", lw=1, alpha=0.5)
        ax.axhline(bonf_thresh, ls="--", c="red", lw=1, alpha=0.5,
                   label=f"Bonferroni (p={0.05/n_tests:.1e})")

        title = "Volcano Plot"
        if anc:
            title += f" — {anc.upper()}"
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("Cohen's d", fontsize=11)
        ax.set_ylabel("-log10(p)", fontsize=11)

        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(
            by_label.values(), by_label.keys(),
            fontsize=9, loc="upper right", framealpha=0.9,
        )

    fig.tight_layout()
    caption = (
        f"Each point is one (interval, metric, group-pair) comparison. "
        f"Vertical lines at d = +/-{data.effect_threshold}. "
        f"Points labeled when |d| > 0.3 and above Bonferroni threshold."
    )
    return embed_figure(fig, caption)


def chart_cross_ancestry(data: ReportData) -> str:
    """Scatter of Cohen's d in one ancestry vs another to show consistency."""
    df = data.comparisons.copy()
    if "ancestry" not in df.columns or df["ancestry"].nunique() < 2:
        return '<div class="alert alert-info">Cross-ancestry plot requires multiple ancestries.</div>'

    parsed = df["metric"].apply(parse_metric_name)
    df["interval"] = [p[0] for p in parsed]
    df["metric_type"] = [p[1] for p in parsed]
    df["group_pair"] = df["group_x"] + " vs " + df["group_y"]
    df["comparison_key"] = df["group_pair"] + " | " + df["metric"]

    # Focus on key metrics
    df = df[df["metric_type"].isin(KEY_METRIC_TYPES)]

    ancestries = sorted(df["ancestry"].unique())
    if len(ancestries) < 2:
        return '<div class="alert alert-info">Cross-ancestry plot requires multiple ancestries.</div>'

    # Use first two ancestries (typically the two largest groups)
    anc_x_name, anc_y_name = ancestries[0], ancestries[1]
    dx = df[df["ancestry"] == anc_x_name].set_index("comparison_key")["cohens_d"]
    dy = df[df["ancestry"] == anc_y_name].set_index("comparison_key")["cohens_d"]
    common = dx.index.intersection(dy.index)

    if len(common) == 0:
        return '<div class="alert alert-info">No matching comparisons across ancestries.</div>'

    # Build merged frame for plotting
    merged = pd.DataFrame({
        f"d_{anc_x_name}": dx.loc[common],
        f"d_{anc_y_name}": dy.loc[common],
    })

    # Add interval info for coloring
    key_to_interval = df.drop_duplicates("comparison_key").set_index("comparison_key")["interval"]
    key_to_pair = df.drop_duplicates("comparison_key").set_index("comparison_key")["group_pair"]
    merged["interval"] = key_to_interval.loc[common].values
    merged["group_pair"] = key_to_pair.loc[common].values

    intervals = sorted(merged["interval"].unique())
    BOLD_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                   "#a65628", "#f781bf", "#999999"]
    interval_palette = dict(zip(intervals, BOLD_COLORS[:len(intervals)]))

    group_pairs = sorted(merged["group_pair"].unique())
    MARKERS = ["o", "s", "D", "^", "v", "P", "X", "*"]
    pair_markers = dict(zip(group_pairs, MARKERS[:len(group_pairs)]))

    fig, ax = plt.subplots(figsize=(8, 7))

    for gp in group_pairs:
        for iv in intervals:
            sub = merged[(merged["group_pair"] == gp) & (merged["interval"] == iv)]
            if sub.empty:
                continue
            ax.scatter(
                sub[f"d_{anc_x_name}"], sub[f"d_{anc_y_name}"],
                c=interval_palette[iv], marker=pair_markers[gp],
                s=80, alpha=0.85, edgecolors="black", linewidth=0.4,
                zorder=3,
            )

    # Diagonal reference
    lim_min = min(merged[f"d_{anc_x_name}"].min(), merged[f"d_{anc_y_name}"].min()) - 0.1
    lim_max = max(merged[f"d_{anc_x_name}"].max(), merged[f"d_{anc_y_name}"].max()) + 0.1
    ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1, alpha=0.4, label="Perfect consistency")
    ax.axhspan(-0.2, 0.2, color="green", alpha=0.05)
    ax.axvspan(-0.2, 0.2, color="green", alpha=0.05)
    ax.axhline(0, color="grey", lw=0.5)
    ax.axvline(0, color="grey", lw=0.5)

    ax.set_xlabel(f"Cohen's d — {anc_x_name.upper()}", fontsize=12)
    ax.set_ylabel(f"Cohen's d — {anc_y_name.upper()}", fontsize=12)
    ax.set_title("Cross-Ancestry Consistency", fontsize=14, fontweight="bold")

    # Label outliers far from diagonal
    for idx, row in merged.iterrows():
        d_diff = abs(row[f"d_{anc_x_name}"] - row[f"d_{anc_y_name}"])
        d_max = max(abs(row[f"d_{anc_x_name}"]), abs(row[f"d_{anc_y_name}"]))
        if d_max > 0.3 or d_diff > 0.3:
            ax.annotate(
                f"{row['interval']}\n{row['group_pair']}",
                (row[f"d_{anc_x_name}"], row[f"d_{anc_y_name}"]),
                fontsize=6.5, ha="center", va="bottom",
                xytext=(0, 6), textcoords="offset points",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="grey", alpha=0.8, lw=0.5),
            )

    # Build legend: colors for intervals, markers for group pairs
    interval_handles = [plt.Line2D([0], [0], marker="o", color="w",
                        markerfacecolor=interval_palette[iv], markersize=9,
                        markeredgecolor="black", markeredgewidth=0.4, label=iv)
                        for iv in intervals]
    pair_handles = [plt.Line2D([0], [0], marker=pair_markers[gp], color="w",
                    markerfacecolor="grey", markersize=8,
                    markeredgecolor="black", markeredgewidth=0.4, label=gp)
                    for gp in group_pairs]
    ax.legend(handles=interval_handles + pair_handles, fontsize=8,
              loc="upper left", framealpha=0.9, ncol=1)

    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()

    # Correlation
    r = merged[f"d_{anc_x_name}"].corr(merged[f"d_{anc_y_name}"])
    caption = (
        f"Each point is one (interval, metric, group-pair) comparison. "
        f"Points on the diagonal indicate identical effects in both ancestries (r = {r:.2f}). "
        f"Color = interval, shape = group pair. Green bands mark |d| < 0.2."
    )
    return embed_figure(fig, caption)


def chart_grouped_bars(data: ReportData) -> str:
    """Grouped bar charts of mean metrics from group_summaries."""
    gs = data.group_summaries.copy()
    comp_name = data.comparison_name
    if gs.empty:
        return '<div class="alert alert-info">No group summary data available.</div>'

    # Discover intervals and metric types from columns
    mean_cols = [c for c in gs.columns if c.endswith("_mean")]
    parsed_cols = [(c, *parse_metric_name(c.replace("_mean", ""))) for c in mean_cols]

    # Filter to key metric types
    key_cols = [(c, iv, mt) for c, iv, mt in parsed_cols if mt in KEY_METRIC_TYPES]
    if not key_cols:
        return '<div class="alert alert-info">No key metrics found in group summaries.</div>'

    intervals = sorted(set(iv for _, iv, _ in key_cols))
    has_ancestry = "ancestry" in gs.columns
    groups = sorted(gs[comp_name].unique())
    group_colors = _get_group_colors(groups)

    html = ""
    for metric_type in KEY_METRIC_TYPES:
        relevant = [(c, iv) for c, iv, mt in key_cols if mt == metric_type]
        if not relevant:
            continue

        n_intervals = len(relevant)
        if has_ancestry:
            ancestries = sorted(gs["ancestry"].dropna().unique())
            n_anc = len(ancestries)
        else:
            ancestries = [None]
            n_anc = 1

        fig, axes = plt.subplots(
            n_anc, n_intervals,
            figsize=(3.5 * n_intervals, 3 * n_anc),
            squeeze=False,
        )

        for ai, anc in enumerate(ancestries):
            sub = gs[gs["ancestry"] == anc] if anc else gs
            for ii, (mean_col, interval) in enumerate(relevant):
                ax = axes[ai, ii]
                std_col = mean_col.replace("_mean", "_std")
                x = np.arange(len(groups))
                width = 0.7

                means = []
                stds = []
                colors = []
                for group in groups:
                    row = sub[sub[comp_name] == group]
                    if len(row) == 0:
                        means.append(0)
                        stds.append(0)
                    else:
                        means.append(row[mean_col].values[0])
                        stds.append(row[std_col].values[0] if std_col in row.columns else 0)
                    colors.append(group_colors[group])

                ax.bar(x, means, width, yerr=stds, color=colors,
                       capsize=3, edgecolor="white", linewidth=0.5, alpha=0.85)
                ax.set_xticks(x)
                ax.set_xticklabels(groups, fontsize=7, rotation=30, ha="right")

                if ai == 0:
                    ax.set_title(interval, fontsize=9, fontweight="bold")
                if ii == 0:
                    label = display_metric(metric_type)
                    if anc:
                        label = f"{anc.upper()}\n{display_metric(metric_type)}"
                    ax.set_ylabel(label, fontsize=9)

        comp_display = comp_name.replace("_", " ").title()
        fig.suptitle(f"Mean {display_metric(metric_type)} by {comp_display}", fontsize=12, fontweight="bold")
        fig.tight_layout()
        html += embed_figure(fig, f"Mean {display_metric(metric_type)} (+/- 1 SD) by {comp_display} and interval.")

    return html


def chart_distributions(data: ReportData) -> str:
    """Violin plots of per-sample metrics (requires sample_stats.tsv)."""
    if data.sample_stats is None:
        return '<div class="alert alert-info">sample_stats.tsv not loaded — distribution plots skipped.</div>'

    ss = data.sample_stats.copy()
    comp_name = data.comparison_name
    has_ancestry = "ancestry" in ss.columns

    # Subsample for rendering speed
    max_per_group = 5000
    if has_ancestry:
        ss = ss.groupby([comp_name, "ancestry"], group_keys=False).apply(
            lambda g: g.sample(n=min(len(g), max_per_group), random_state=42)
        )
    else:
        ss = ss.groupby(comp_name, group_keys=False).apply(
            lambda g: g.sample(n=min(len(g), max_per_group), random_state=42)
        )

    # Find metric columns for key metrics
    id_cols = {"s", comp_name, "ancestry", "ancestry_pred_other"}
    metric_cols = [c for c in ss.columns if c not in id_cols]
    parsed = [(c, *parse_metric_name(c)) for c in metric_cols]
    key = [(c, iv, mt) for c, iv, mt in parsed if mt in KEY_METRIC_TYPES]

    if not key:
        return '<div class="alert alert-info">No key metrics found in sample stats.</div>'

    intervals = sorted(set(iv for _, iv, _ in key))
    groups = sorted(ss[comp_name].unique())
    group_colors = _get_group_colors(groups)
    palette = {g: group_colors[g] for g in groups}

    html = ""
    for metric_type in KEY_METRIC_TYPES:
        cols_for_mt = [(c, iv) for c, iv, mt in key if mt == metric_type]
        if not cols_for_mt:
            continue

        n_iv = len(cols_for_mt)
        if has_ancestry:
            ancestries = sorted(ss["ancestry"].dropna().unique())
            n_anc = len(ancestries)
        else:
            ancestries = [None]
            n_anc = 1

        fig, axes = plt.subplots(
            n_anc, n_iv,
            figsize=(3.5 * n_iv, 3 * n_anc),
            squeeze=False,
        )

        for ai, anc in enumerate(ancestries):
            sub = ss[ss["ancestry"] == anc] if anc else ss
            for ii, (col, interval) in enumerate(cols_for_mt):
                ax = axes[ai, ii]
                plot_df = sub[[comp_name, col]].dropna()
                if plot_df.empty:
                    ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
                    continue

                sns.violinplot(
                    data=plot_df, x=comp_name, y=col, hue=comp_name,
                    palette=palette, ax=ax, inner="quartile",
                    linewidth=0.8, cut=0, density_norm="width", legend=False,
                )
                ax.set_xlabel("")
                ax.tick_params(axis="x", rotation=30, labelsize=7)

                if ai == 0:
                    ax.set_title(interval, fontsize=9, fontweight="bold")
                if ii == 0:
                    label = display_metric(metric_type)
                    if anc:
                        label = f"{anc.upper()}\n{display_metric(metric_type)}"
                    ax.set_ylabel(label, fontsize=9)
                else:
                    ax.set_ylabel("")

        comp_display = comp_name.replace("_", " ").title()
        fig.suptitle(f"Distribution of {display_metric(metric_type)} by {comp_display}", fontsize=12, fontweight="bold")
        fig.tight_layout()
        html += embed_figure(fig, f"Violin plots of per-sample {display_metric(metric_type)}. Inner lines show quartiles.")

    return html


def chart_pca(data: ReportData) -> str:
    """PCA scatter of sample metrics, colored by comparison group."""
    if data.sample_stats is None:
        return '<div class="alert alert-info">sample_stats.tsv not loaded — PCA plots skipped.</div>'

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    ss = data.sample_stats.copy()
    comp_name = data.comparison_name
    has_ancestry = "ancestry" in ss.columns

    id_cols = {"s", comp_name, "ancestry", "ancestry_pred_other"}
    metric_cols = [c for c in ss.columns if c not in id_cols]

    if len(metric_cols) < 3:
        return '<div class="alert alert-info">Too few metrics for PCA.</div>'

    groups = sorted(ss[comp_name].unique())
    group_colors = _get_group_colors(groups)

    if has_ancestry:
        ancestries = sorted(ss["ancestry"].dropna().unique())
    else:
        ancestries = [None]

    n_anc = len(ancestries)
    fig, axes = plt.subplots(1, n_anc, figsize=(6.5 * n_anc, 6), squeeze=False)

    for i, anc in enumerate(ancestries):
        ax = axes[0, i]
        sub = ss[ss["ancestry"] == anc].copy() if anc else ss.copy()

        # Subsample for speed
        if len(sub) > 5000:
            sub = sub.sample(n=5000, random_state=42)

        X = sub[metric_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        if X.shape[0] < 10:
            ax.text(0.5, 0.5, "Too few samples", ha="center", va="center", transform=ax.transAxes)
            continue

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)

        for group in groups:
            mask = sub[comp_name].values == group
            if mask.sum() == 0:
                continue
            ax.scatter(
                pcs[mask, 0], pcs[mask, 1],
                c=[group_colors[group]], label=group,
                s=18, alpha=0.45, edgecolors="black", linewidths=0.3,
            )

        var1 = pca.explained_variance_ratio_[0] * 100
        var2 = pca.explained_variance_ratio_[1] * 100
        ax.set_xlabel(f"PC1 ({var1:.1f}%)", fontsize=10)
        ax.set_ylabel(f"PC2 ({var2:.1f}%)", fontsize=10)

        title = "PCA of Sample Metrics"
        if anc:
            title += f" — {anc.upper()}"
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.legend(fontsize=9, framealpha=0.8, markerscale=3)

    fig.tight_layout()
    caption = "PCA on standardized per-sample metrics. Computed within ancestry to avoid population structure dominating PCs."
    return embed_figure(fig, caption)


def chart_timing(data: ReportData) -> str:
    """Horizontal bar chart of pipeline step durations."""
    if data.timing is None:
        return '<div class="alert alert-info">timing_metrics.json not available — performance section skipped.</div>'

    timings = data.timing.get("step_timings", {})
    if not timings:
        return '<div class="alert alert-info">No step timings found.</div>'

    # Exclude total from bar chart
    step_timings = {k: v for k, v in timings.items() if k != "total_pipeline"}
    if not step_timings:
        return ""

    labels = list(step_timings.keys())
    values = [step_timings[k] for k in labels]
    # Clean up step names
    display_labels = [l.replace("step", "Step ").replace("_", " ").strip() for l in labels]

    fig, ax = plt.subplots(figsize=(8, max(3, len(labels) * 0.5)))
    colors = sns.color_palette("Blues_r", len(labels))
    bars = ax.barh(range(len(labels)), values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(display_labels, fontsize=9)
    ax.set_xlabel("Duration (seconds)", fontsize=10)
    ax.set_title("Pipeline Step Durations", fontsize=12, fontweight="bold")
    ax.invert_yaxis()

    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:.0f}s", va="center", fontsize=8)

    fig.tight_layout()

    html = embed_figure(fig, "Wall-clock duration of each pipeline step.")

    # Cluster info card
    cluster = data.timing.get("cluster_info", {})
    if cluster and "error" not in cluster:
        html += '<div class="card"><h3>Cluster Configuration</h3><table>'
        for k, v in cluster.items():
            html += f"<tr><td><strong>{k}</strong></td><td>{v}</td></tr>"
        html += "</table></div>"

    total = timings.get("total_pipeline")
    if total:
        html += f'<div class="card"><strong>Total pipeline time:</strong> {total:.0f}s ({total/60:.1f} min)</div>'

    return html


# ---------------------------------------------------------------------------
# HTML Assembly
# ---------------------------------------------------------------------------

def _render_header(data: ReportData, title: str) -> str:
    """Render the report header and run metadata section."""
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    html = f"<h1>{title}</h1>"
    html += f'<p style="color:#666">Generated {now} from <code>{data.input_path}</code></p>'

    comp_display = data.comparison_name.replace("_", " ").title()

    # Config summary card
    if data.config:
        cfg = data.config
        html += '<div class="card"><h3>Run Configuration</h3>'
        html += '<div class="grid-2">'
        fields = [
            ("Comparison", f"{cfg.get('comparison_name', '')} (col: {cfg.get('comparison_col', '')})"),
            ("Samples/Group", cfg.get("samples_per_group")),
            ("Ancestries", ", ".join(cfg.get("ancestries", []) or ["all"])),
            ("Comparison Values", ", ".join(cfg.get("comparison_values", []) or ["all"])),
            ("Intervals", ", ".join((cfg.get("interval_dict") or {}).keys())),
            ("Run Name", cfg.get("run_name")),
            ("Filter to PASS", cfg.get("filter_to_pass")),
        ]
        for label, value in fields:
            if value is not None and value != "":
                html += f'<div class="summary-stat"><span class="label">{label}</span><br>'
                html += f'<span class="value">{value}</span></div>'
        html += "</div></div>"

    # Sample size table from comparisons
    df = data.comparisons
    if not df.empty:
        html += f'<div class="card"><h3>Sample Sizes</h3>'
        # Extract unique (group, ancestry) counts from n_x / n_y
        rows = []
        has_anc = "ancestry" in df.columns
        for _, r in df.iterrows():
            anc = r.get("ancestry", "") if has_anc else ""
            rows.append({comp_display: r["group_x"], "ancestry": anc, "n": r["n_x"]})
            rows.append({comp_display: r["group_y"], "ancestry": anc, "n": r["n_y"]})

        sizes = pd.DataFrame(rows).drop_duplicates()
        if has_anc:
            pivot = sizes.pivot_table(index=comp_display, columns="ancestry", values="n", aggfunc="first")
            pivot = pivot.reindex(sorted(pivot.columns), axis=1)
        else:
            pivot = sizes.groupby(comp_display)["n"].first().to_frame()

        html += pivot.to_html(classes="", na_rep="-", float_format=lambda x: f"{x:.0f}")
        html += "</div>"

    return html


def _render_methods(data: ReportData) -> str:
    """Render methods and interpretation guide."""
    html = '<div class="methods">'
    html += "<h2>Methods &amp; Interpretation Guide</h2>"

    html += "<h3>Statistical Methods</h3>"
    html += "<ul>"
    html += "<li><strong>Cohen's d</strong>: Standardized mean difference = (mean_x - mean_y) / pooled_SD, "
    html += "where pooled_SD = sqrt(((n_x-1)*s_x&sup2; + (n_y-1)*s_y&sup2;) / (n_x+n_y-2)).</li>"
    html += "<li><strong>Welch's t-test</strong>: Two-sample t-test not assuming equal variances.</li>"
    html += "<li><strong>Multiple testing</strong>: Bonferroni threshold shown on volcano plots. "
    html += "Individual p-values should be interpreted cautiously given the number of comparisons.</li>"
    html += "</ul>"

    html += "<h3>Effect Size Classification</h3>"
    html += "<table>"
    html += "<tr><th>|Cohen's d|</th><th>Interpretation</th><th>Color</th></tr>"
    html += f'<tr style="background:{SEVERITY_COLORS["not appreciable"]}">'
    html += "<td>&lt; 0.5</td><td>Not appreciable</td><td></td></tr>"
    html += f'<tr style="background:{SEVERITY_COLORS["appreciable"]}">'
    html += "<td>&ge; 0.5</td><td>Appreciable batch effect</td><td></td></tr>"
    html += "</table>"
    html += "<p>The |d| &ge; 0.5 threshold corresponds to Cohen's convention for a "
    html += "medium-sized effect and was chosen before looking at any data. "
    html += "Randomized control experiments show a noise floor of |d| &lt; 0.35.</p>"

    # Interval descriptions from config if available
    if data.config and data.config.get("interval_dict"):
        html += "<h3>Genomic Intervals</h3>"
        html += "<table><tr><th>Interval</th><th>Path</th></tr>"
        for name, path in data.config["interval_dict"].items():
            html += f"<tr><td><strong>{name}</strong></td><td><code>{path}</code></td></tr>"
        html += "</table>"

    html += "<h3>Caveats</h3>"
    html += "<ul>"
    html += "<li>Cohen's d uses pooled standard deviation as denominator.</li>"
    html += "<li>Large sample sizes inflate statistical significance — focus on effect sizes.</li>"
    html += "<li>Within-ancestry stratification controls for population structure but reduces power for small ancestry groups.</li>"
    html += "</ul>"
    html += "</div>"
    return html


def render_report(data: ReportData, title: str = "Batch Effect Report") -> str:
    """Assemble the full HTML report."""
    sections = []

    # Header
    sections.append(_render_header(data, title))

    # Executive Summary
    sections.append("<h2>Executive Summary</h2>")
    sections.append(safe_render(chart_significance_table, data))

    # Heatmap
    sections.append("<h2>Effect Size Heatmaps</h2>")
    sections.append(safe_render(chart_effect_heatmap, data))

    # Volcano
    sections.append("<h2>Volcano Plots</h2>")
    sections.append(safe_render(chart_volcano, data))

    # Cross-ancestry
    sections.append("<h2>Cross-Ancestry Consistency</h2>")
    sections.append(safe_render(chart_cross_ancestry, data))

    # Grouped bars
    comp_display = data.comparison_name.replace("_", " ").title()
    sections.append(f"<h2>Group Means by {comp_display}</h2>")
    sections.append(safe_render(chart_grouped_bars, data))

    # Distributions (conditional)
    sections.append("<h2>Per-Sample Distributions</h2>")
    sections.append(safe_render(chart_distributions, data))

    # PCA (conditional)
    sections.append("<h2>PCA of Sample Metrics</h2>")
    sections.append(safe_render(chart_pca, data))

    # Timing (conditional)
    sections.append("<h2>Pipeline Performance</h2>")
    sections.append(safe_render(chart_timing, data))

    # Methods
    sections.append(safe_render(_render_methods, data))

    body = "\n".join(sections)

    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>{REPORT_CSS}</style>
</head>
<body>
{body}
<hr>
<p style="text-align:center; color:#999; font-size:0.85em;">
    Generated by batch_e_reporter.py &mdash; {now}
</p>
</body>
</html>"""


def generate_report(
    input_path: str,
    output_path: Optional[str] = None,
    title: str = "Batch Effect Report",
    load_sample_stats: bool = True,
    effect_threshold: float = 0.5,
):
    """Load results and generate an HTML report.

    In a Jupyter notebook this returns an IPython HTML display object that
    renders inline. From the CLI or plain Python it returns the HTML string.
    """
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="All-NaN")

    logger.info("=" * 60)
    logger.info("Batch Effect Reporter")
    logger.info("=" * 60)

    data = load_report_data(
        input_path=input_path,
        load_sample_stats=load_sample_stats,
        effect_threshold=effect_threshold,
    )

    logger.info("Rendering report ...")
    html = render_report(data, title=title)
    logger.info(f"Report size: {len(html) / 1024:.0f} KB")

    # Save a copy if requested
    if output_path:
        if output_path.startswith("gs://"):
            with _open_file(output_path, "w") as f:
                f.write(html)
        else:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            with open(output_path, "w") as f:
                f.write(html)
        logger.info(f"Report saved to {output_path}")

    try:
        from IPython.display import HTML as IPyHTML
        return IPyHTML(_wrap_with_download_link(html, title))
    except ImportError:
        return html


def _wrap_with_download_link(html: str, title: str) -> str:
    """Wrap report HTML with a base64 data-URI download link at the top."""
    b64 = base64.b64encode(html.encode("utf-8")).decode("utf-8")
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    filename = f"batch_effect_report_{timestamp}.html"

    download_link = (
        f'<p style="text-align:center; margin-top:12px; font-style:italic; font-size:0.9em; color:#666;">'
        f'<a download="{filename}" '
        f'href="data:text/html;base64,{b64}" '
        f'style="color:#1565c0; text-decoration:none; border-bottom:1px dotted #1565c0;">'
        f'Download this report as HTML</a></p>'
        f'<p style="text-align:center; margin-top:4px; font-size:0.8em; color:#999;">'
        f'Tip: open the downloaded file in Chrome and use Print \u2192 Save as PDF for best results.</p>'
    )

    return html + download_link


# ---------------------------------------------------------------------------
# CLI Entry Point
# ---------------------------------------------------------------------------

def cli():
    """Command-line interface for the batch effect reporter."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate HTML report from batch effect pipeline results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python batch_e_reporter.py ./results/\n"
            "  python batch_e_reporter.py gs://bucket/batch_effect_results/run/ -o report.html\n"
            "  python batch_e_reporter.py ./results/ --no-sample-stats --effect-threshold 0.5\n"
        ),
    )
    parser.add_argument(
        "input_path",
        help="Local or GCS directory containing pipeline output files",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output HTML path (default: report.html in current directory)",
    )
    parser.add_argument(
        "--title",
        default="Batch Effect Report",
        help="Report title",
    )
    parser.add_argument(
        "--no-sample-stats",
        action="store_true",
        help="Skip loading sample_stats.tsv (faster, skips violins and PCA)",
    )
    parser.add_argument(
        "--effect-threshold",
        type=float,
        default=0.5,
        help="Cohen's d threshold for flagging (default: 0.5)",
    )

    args = parser.parse_args()

    output = args.output or "report.html"

    generate_report(
        input_path=args.input_path,
        output_path=output,
        title=args.title,
        load_sample_stats=not args.no_sample_stats,
        effect_threshold=args.effect_threshold,
    )
    print(f"Report written to: {output}")


if __name__ == "__main__":
    cli()
