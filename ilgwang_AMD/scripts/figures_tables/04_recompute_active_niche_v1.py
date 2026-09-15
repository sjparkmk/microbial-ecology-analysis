\
\
\
\
\


from __future__ import annotations

import argparse
import csv
import hashlib
import platform
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Circle


MM = 1 / 25.4
A4_LANDSCAPE = (297 * MM, 210 * MM)
A4_PORTRAIT = (210 * MM, 297 * MM)

SAMPLES = [
    "1st_BFP",
    "1st_BFW1",
    "1st_BFW2",
    "2nd_BFP",
    "2nd_BFW",
    "2nd_wt",
    "3rd_wt",
]
TYPE_SAMPLES = {
    "P_biofilm": ["1st_BFP", "2nd_BFP"],
    "W_biofilm": ["1st_BFW1", "1st_BFW2", "2nd_BFW"],
    "water": ["2nd_wt", "3rd_wt"],
}
TYPE_LABELS = {
    "P_biofilm": "AMD precipitate-associated\nbiofilm",
    "W_biofilm": "Wall-attached\nbiofilm",
    "water": "Water",
}
TYPE_COLORS = {
    "P_biofilm": "#1BAF7A",
    "W_biofilm": "#EB6834",
    "water": "#2A78D6",
}
INK = "#111111"
SECONDARY = "#5E6462"
GRID = "#E4E7E5"
BASELINE = "#A8AAA7"
TEAL = "#1F6F5F"


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.labelsize": 10.5,
        "axes.labelweight": "bold",
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "axes.edgecolor": "#000000",
        "axes.linewidth": 1.15,
        "xtick.color": "#000000",
        "ytick.color": "#000000",
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)


def read_rows(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_rows(path: Path, rows: list[dict], fields: list[str], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_axes(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")
    ax.spines["left"].set_linewidth(1.15)
    ax.spines["bottom"].set_linewidth(1.15)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.65, zorder=0)
        ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for suffix, kwargs in (
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": 320}),
        ("tiff", {"dpi": 400, "pil_kwargs": {"compression": "tiff_lzw"}}),
    ):
        path = out_dir / f"{stem}.{suffix}"
        fig.savefig(path, bbox_inches=None, facecolor="white", **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def levins(values: np.ndarray) -> tuple[float, float]:
    total = float(values.sum())
    if total <= 0:
        raise ValueError("Levins input has no abundance")
    proportions = values / total
    b_value = float(1.0 / np.square(proportions).sum())
    ba_value = float((b_value - 1.0) / (len(values) - 1.0))
    return b_value, ba_value


def build_niche_rows(
    seqtab_path: Path,
    taxonomy_path: Path,
) -> tuple[list[dict], dict[str, float], list[dict]]:
    seqtab = read_rows(seqtab_path)
    if len(seqtab) != 663:
        raise ValueError(f"Expected 663 ASVs; observed {len(seqtab)}")

    taxonomy_rows = read_rows(taxonomy_path)
    taxonomy = {row["sequence"]: row for row in taxonomy_rows}
    if len(taxonomy) != len(taxonomy_rows):
        raise ValueError("Taxonomy sequence keys are not unique")

    library_depths = np.zeros(len(SAMPLES), dtype=np.int64)
    asv_records: list[dict] = []
    genus_counts: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(len(SAMPLES), dtype=np.int64))
    genus_kingdom: dict[str, set[str]] = defaultdict(set)
    genus_phylum: dict[str, set[str]] = defaultdict(set)

    for row in seqtab:
        counts = np.array([int(row[sample]) for sample in SAMPLES], dtype=np.int64)
        library_depths += counts
        sequence = row["sequence"]
        if sequence not in taxonomy:
            raise ValueError(f"ASV sequence missing from taxonomy: {row['ASV_ID']}")
        taxon = taxonomy[sequence]
        genus_raw = taxon.get("Genus", "").strip()
        genus = "" if genus_raw.upper() in {"", "NA", "NAN", "NULL"} else genus_raw
        asv_records.append(
            {
                "ASV_ID": row["ASV_ID"],
                "counts": counts,
                "total_reads": int(counts.sum()),
            }
        )
        if genus:
            genus_counts[genus] += counts
            if taxon.get("Kingdom", "").strip():
                genus_kingdom[genus].add(taxon["Kingdom"].strip())
            if taxon.get("Phylum", "").strip():
                genus_phylum[genus].add(taxon["Phylum"].strip())

    if int(library_depths.sum()) != 716_915:
        raise ValueError(f"Expected 716,915 reads; observed {int(library_depths.sum())}")
    if len(genus_counts) != 129:
        raise ValueError(f"Expected 129 named genera; observed {len(genus_counts)}")

    sample_index = {sample: idx for idx, sample in enumerate(SAMPLES)}
    output: list[dict] = []
    for genus in sorted(genus_counts):
        counts = genus_counts[genus]
        sample_ra = counts / library_depths * 100.0
        type_means = np.array(
            [
                np.mean([sample_ra[sample_index[s]] for s in TYPE_SAMPLES[code]])
                for code in TYPE_SAMPLES
            ],
            dtype=float,
        )
        pooled_counts = np.array(
            [
                np.sum([counts[sample_index[s]] for s in TYPE_SAMPLES[code]])
                for code in TYPE_SAMPLES
            ],
            dtype=float,
        )
        primary_b, primary_ba = levins(type_means)
        pooled_b, pooled_ba = levins(pooled_counts)
        output.append(
            {
                "Genus": genus,
                "Kingdom": "; ".join(sorted(genus_kingdom[genus])),
                "Phylum": "; ".join(sorted(genus_phylum[genus])),
                "total_reads": int(counts.sum()),
                "n_samples_present": int(np.sum(counts > 0)),
                "mean_relative_abundance_percent": float(sample_ra.mean()),
                "type_mean_percent_P_biofilm": float(type_means[0]),
                "type_mean_percent_W_biofilm": float(type_means[1]),
                "type_mean_percent_water": float(type_means[2]),
                "Levins_B_primary_equal_sample_weight": primary_b,
                "Levins_BA_primary_equal_sample_weight": primary_ba,
                "Levins_B_pooled_read_sensitivity": pooled_b,
                "Levins_BA_pooled_read_sensitivity": pooled_ba,
                "absolute_BA_difference": abs(primary_ba - pooled_ba),
            }
        )
    output.sort(key=lambda row: (-row["mean_relative_abundance_percent"], row["Genus"]))

    deltas = np.array([row["absolute_BA_difference"] for row in output])
    audit = {
        "n_genera": float(len(output)),
        "median_absolute_BA_difference": float(np.median(deltas)),
        "maximum_absolute_BA_difference": float(np.max(deltas)),
    }
    return output, audit, asv_records


def build_overlap_rows(asv_records: list[dict]) -> tuple[list[dict], dict[str, int]]:
    sample_index = {sample: idx for idx, sample in enumerate(SAMPLES)}
    values: dict[str, int] = defaultdict(int)
    type_totals = {code: 0 for code in TYPE_SAMPLES}
    for record in asv_records:
        counts = record["counts"]
        present = {
            code: bool(sum(counts[sample_index[s]] for s in samples) > 0)
            for code, samples in TYPE_SAMPLES.items()
        }
        for code, status in present.items():
            type_totals[code] += int(status)
        key = (present["P_biofilm"], present["W_biofilm"], present["water"])
        region_map = {
            (True, False, False): "P_biofilm_only",
            (False, True, False): "W_biofilm_only",
            (False, False, True): "water_only",
            (True, True, False): "P_biofilm_and_W_biofilm_only",
            (True, False, True): "P_biofilm_and_water_only",
            (False, True, True): "W_biofilm_and_water_only",
            (True, True, True): "all_three",
        }
        if key not in region_map:
            raise ValueError(f"ASV absent from all Types: {record['ASV_ID']}")
        values[region_map[key]] += 1

    order = [
        "P_biofilm_only",
        "W_biofilm_only",
        "water_only",
        "P_biofilm_and_W_biofilm_only",
        "P_biofilm_and_water_only",
        "W_biofilm_and_water_only",
        "all_three",
    ]
    expected = [102, 127, 133, 84, 24, 31, 162]
    observed = [values[key] for key in order]
    if observed != expected or sum(observed) != 663:
        raise ValueError(f"Unexpected ASV overlap values: {observed}")
    if type_totals != {"P_biofilm": 372, "W_biofilm": 404, "water": 350}:
        raise ValueError(f"Unexpected Type detection totals: {type_totals}")
    rows = [{"region": key, "ASV_count": values[key]} for key in order]
    totals = {"any_P_biofilm": 372, "any_W_biofilm": 404, "any_water": 350, "all_ASVs": 663}
    return rows, totals


def build_figure4_audit(asv_records: list[dict], tips_path: Path) -> list[dict]:
    total_reads = sum(record["total_reads"] for record in asv_records)
    displayed = [
        {
            "ASV_ID": record["ASV_ID"],
            "total_reads": record["total_reads"],
            "pooled_relative_abundance": record["total_reads"] / total_reads,
        }
        for record in asv_records
        if record["total_reads"] / total_reads >= 0.002
    ]
    displayed.sort(key=lambda row: -row["pooled_relative_abundance"])
    tips = read_rows(tips_path, "\t")
    if len(displayed) != 56 or {row["ASV_ID"] for row in displayed} != {row["ASV_ID"] for row in tips}:
        raise ValueError("Figure 4 tip set does not match the canonical 0.2% display filter")
    return displayed


def figure_continuous_niche(rows: list[dict], figure_dir: Path) -> list[Path]:
    top_rows = sorted(rows, key=lambda row: -row["mean_relative_abundance_percent"])[:15]
    top_names = {row["Genus"] for row in top_rows}
    other_rows = [row for row in rows if row["Genus"] not in top_names]

    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(
        1, 2, width_ratios=[1.18, 1.0], left=0.085, right=0.975, bottom=0.17, top=0.82, wspace=0.42
    )
    ax1 = fig.add_subplot(grid[0, 0])
    ax2 = fig.add_subplot(grid[0, 1])

    ax1.scatter(
        [row["mean_relative_abundance_percent"] for row in other_rows],
        [row["Levins_BA_primary_equal_sample_weight"] for row in other_rows],
        s=22, color=BASELINE, edgecolor="white", linewidth=0.4, alpha=0.78, zorder=2,
    )
    ax1.scatter(
        [row["mean_relative_abundance_percent"] for row in top_rows],
        [row["Levins_BA_primary_equal_sample_weight"] for row in top_rows],
        s=36, color=TEAL, edgecolor="white", linewidth=0.55, alpha=0.95, zorder=3,
    )
    ax1.set_xscale("log")
    ax1.set_ylim(-0.02, 1.02)
    ax1.set_xlabel("Mean relative sequence abundance across samples (%)")
    ax1.set_ylabel("Standardized Levins’ niche breadth (Bₐ)")
    ax1.set_title("a  All named genera (n = 129)", loc="left", pad=10)
    clean_axes(ax1, grid_axis="y")

    ranked = sorted(top_rows, key=lambda row: row["Levins_BA_primary_equal_sample_weight"])
    y_values = np.arange(len(ranked))
    b_values = np.array([row["Levins_BA_primary_equal_sample_weight"] for row in ranked])
    ax2.hlines(y_values, 0, b_values, color="#D8D8D4", linewidth=1.0, zorder=1)
    ax2.scatter(b_values, y_values, s=38, color=TEAL, edgecolor="white", linewidth=0.55, zorder=3)
    ax2.set_yticks(y_values, [row["Genus"] for row in ranked])
    for tick in ax2.get_yticklabels():
        tick.set_fontstyle("italic")
    ax2.set_xlim(0, 1.05)
    ax2.set_xlabel("Standardized Levins’ niche breadth (Bₐ)")
    ax2.set_title("b  Fifteen most abundant named genera", loc="left", pad=10)
    clean_axes(ax2, grid_axis="x")
    for y_value, b_value in zip(y_values, b_values):
        if b_value > 0.94:
            ax2.text(b_value - 0.025, y_value, f"{b_value:.2f}", ha="right", va="center", fontsize=7.7)
        else:
            ax2.text(b_value + 0.022, y_value, f"{b_value:.2f}", ha="left", va="center", fontsize=7.7)

    fig.suptitle(
        "Genus-level niche breadth is reported as a continuous descriptor",
        x=0.085, y=0.935, ha="left", fontsize=12.0, fontweight="bold", color=INK,
    )
    fig.text(
        0.085, 0.885,
        "Primary calculation: equal-weighted mean of sample relative abundances within each of three microhabitat Types",
        ha="left", va="center", fontsize=8.8, color=SECONDARY,
    )
    fig.text(
        0.085, 0.055,
        "Teal marks the 15 genera with the highest mean relative sequence abundance. No generalist or specialist threshold is applied.",
        ha="left", va="center", fontsize=8.2, color=SECONDARY,
    )
    return save_figure(fig, figure_dir, "SFig5_Continuous_Levins_BA_A4_v4")


def figure_aggregation_sensitivity(rows: list[dict], figure_dir: Path) -> list[Path]:
    top_rows = sorted(rows, key=lambda row: -row["mean_relative_abundance_percent"])[:15]
    top_names = {row["Genus"] for row in top_rows}
    other_rows = [row for row in rows if row["Genus"] not in top_names]
    median_delta = float(np.median([row["absolute_BA_difference"] for row in rows]))

    fig = plt.figure(figsize=A4_PORTRAIT)
    ax = fig.add_axes([0.16, 0.24, 0.78, 0.55])
    ax.plot([0, 1], [0, 1], linestyle="--", color="#6F6F6F", linewidth=1.0, zorder=1)
    ax.scatter(
        [row["Levins_BA_primary_equal_sample_weight"] for row in other_rows],
        [row["Levins_BA_pooled_read_sensitivity"] for row in other_rows],
        s=26, color=BASELINE, edgecolor="white", linewidth=0.45, alpha=0.78, zorder=2,
    )
    ax.scatter(
        [row["Levins_BA_primary_equal_sample_weight"] for row in top_rows],
        [row["Levins_BA_pooled_read_sensitivity"] for row in top_rows],
        s=44, color=TEAL, edgecolor="white", linewidth=0.6, alpha=0.96, zorder=3,
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Bₐ from equal-weighted Type means (primary)")
    ax.set_ylabel("Bₐ from pooled reads within Type (sensitivity)")
    clean_axes(ax, grid_axis="both")


    fig.suptitle(
        "Niche-breadth values change with Type aggregation",
        x=0.12, y=0.91, ha="left", fontsize=12.0, fontweight="bold", color=INK,
    )
    fig.text(
        0.12, 0.865,
        f"Same 663-ASV input and 129 named genera; median |ΔBₐ| = {median_delta:.3f}",
        ha="left", va="center", fontsize=8.8, color=SECONDARY,
    )
    legend = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=TEAL, markeredgecolor="white", markersize=7, label="15 most abundant named genera"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=BASELINE, markeredgecolor="white", markersize=6, label="Other named genera"),
        Line2D([0], [0], color="#6F6F6F", linestyle="--", linewidth=1.0, label="Equal values"),
    ]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.53, 0.135), frameon=False, ncol=1, fontsize=8.2)
    fig.text(
        0.50, 0.060,
        "Pooled-read aggregation gives more weight to Types with more samples and to samples with greater read depth.\n"
        "It is shown only as a sensitivity analysis.",
        ha="center", va="center", fontsize=8.0, color=SECONDARY, linespacing=1.35,
    )
    return save_figure(fig, figure_dir, "SFig6_Levins_BA_Aggregation_Sensitivity_A4_v4")


def figure_asv_overlap(overlap_rows: list[dict], figure_dir: Path) -> list[Path]:
    values = {row["region"]: int(row["ASV_count"]) for row in overlap_rows}
    fig = plt.figure(figsize=A4_PORTRAIT)
    ax = fig.add_axes([0.055, 0.15, 0.89, 0.70])
    ax.set_aspect("equal")
    ax.axis("off")
    radius = 1.45
    centers = {"P_biofilm": (-0.72, 0.35), "W_biofilm": (0.72, 0.35), "water": (0.0, -0.70)}
    for code in ("P_biofilm", "W_biofilm", "water"):
        ax.add_patch(
            Circle(
                centers[code], radius, facecolor=TYPE_COLORS[code], edgecolor=TYPE_COLORS[code], linewidth=2.0, alpha=0.44
            )
        )
    number_positions = {
        "P_biofilm_only": (-1.30, 0.50),
        "W_biofilm_only": (1.30, 0.50),
        "water_only": (0.0, -1.62),
        "P_biofilm_and_W_biofilm_only": (0.0, 1.18),
        "P_biofilm_and_water_only": (-0.58, -0.30),
        "W_biofilm_and_water_only": (0.58, -0.30),
        "all_three": (0.0, 0.22),
    }
    for region, position in number_positions.items():
        ax.text(*position, str(values[region]), ha="center", va="center", fontsize=12.0, color=INK)
    ax.text(-1.05, 1.92, TYPE_LABELS["P_biofilm"], ha="center", va="bottom", color=TYPE_COLORS["P_biofilm"], fontsize=9.5, fontweight="bold", linespacing=1.08)
    ax.text(1.05, 1.92, TYPE_LABELS["W_biofilm"], ha="center", va="bottom", color=TYPE_COLORS["W_biofilm"], fontsize=9.5, fontweight="bold", linespacing=1.08)
    ax.text(0.0, -2.25, TYPE_LABELS["water"], ha="center", va="top", color=TYPE_COLORS["water"], fontsize=10.0, fontweight="bold")
    ax.set_xlim(-2.25, 2.25)
    ax.set_ylim(-2.55, 2.30)
    fig.suptitle(
        "ASV detection overlaps across the three microhabitat Types",
        x=0.08, y=0.93, ha="left", fontsize=12.0, fontweight="bold", color=INK,
    )
    fig.text(
        0.08, 0.075,
        "Presence was defined as ≥1 read in ≥1 sample of a Type (n = 663 ASVs). Circle areas are schematic and are not scaled to counts.",
        ha="left", va="center", fontsize=8.2, color=SECONDARY,
    )
    return save_figure(fig, figure_dir, "SFig7_ASV_Overlap_A4_v4")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path("/input/youre/data/"))
    parser.add_argument("--output-root", type=Path, default=Path("/output/your/data"))
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    package_root = args.output_root.resolve()
    output_root = package_root / "build" / "niche_stage"
    figure_dir = output_root / "figures"
    data_dir = package_root / "source_data"
    log_dir = package_root / "logs"
    for directory in (figure_dir, data_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    seqtab_path = source_root / "output_data" / "seqtab_final_long.csv"
    taxonomy_path = source_root / "output_data" / "taxa_silva144.csv"
    inputs = [seqtab_path, taxonomy_path]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    niche_rows, niche_audit, asv_records = build_niche_rows(seqtab_path, taxonomy_path)
    overlap_rows, overlap_totals = build_overlap_rows(asv_records)

    niche_path = data_dir / "SFig5_SFig6_continuous_niche_values_v4.tsv"
    overlap_path = data_dir / "SFig7_ASV_overlap_regions_v4.tsv"
    niche_fields = list(niche_rows[0].keys())
    write_rows(niche_path, niche_rows, niche_fields)
    write_rows(overlap_path, overlap_rows, ["region", "ASV_count"])
    outputs: list[Path] = [niche_path, overlap_path]
    outputs += figure_continuous_niche(niche_rows, figure_dir)
    outputs += figure_aggregation_sensitivity(niche_rows, figure_dir)
    outputs += figure_asv_overlap(overlap_rows, figure_dir)

    log_path = log_dir / "78_niche_overlap_generation_v8.log"
    lines = [
        "Continuous niche-breadth and ASV-overlap figure generation",
        f"UTC timestamp: {datetime.now(timezone.utc).isoformat()}",
        f"Python: {platform.python_version()}",
        f"numpy: {np.__version__}",
        f"matplotlib: {mpl.__version__}",
        "",
        "INPUTS",
    ]
    lines.extend(f"{path}\tSHA256={sha256(path)}" for path in inputs)
    lines.extend(
        [
            "",
            "VALIDATION",
            "Canonical table: 716,915 reads; 663 ASVs; 7 samples",
            f"Named genera: {int(niche_audit['n_genera'])}",
            f"Median absolute BA aggregation difference: {niche_audit['median_absolute_BA_difference']:.10g}",
            f"Maximum absolute BA aggregation difference: {niche_audit['maximum_absolute_BA_difference']:.10g}",
            f"ASV overlap totals: {overlap_totals}",
            "All genus labels and abundance values recomputed directly from taxa_silva144.csv plus seqtab_final_long.csv",
            "Categorical generalist/specialist assignments: not generated",
            "",
            "OUTPUTS",
        ]
    )
    lines.extend(f"{path}\tSHA256={sha256(path)}" for path in outputs)
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"WROTE {log_path}")
    for path in outputs:
        print(f"WROTE {path}")


if __name__ == "__main__":
    main()
