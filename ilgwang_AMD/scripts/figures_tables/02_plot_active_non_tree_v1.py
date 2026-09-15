\
\
\
\
\
\
\


from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch, Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform


MM = 1 / 25.4
A4_LANDSCAPE = (297 * MM, 210 * MM)
A4_PORTRAIT = (210 * MM, 297 * MM)

TYPE_ORDER = ["P_biofilm", "W_biofilm", "water"]
TYPE_LABELS = {
    "P_biofilm": "AMD precipitate-associated biofilm",
    "W_biofilm": "wall-attached biofilm",
    "water": "water",
}
TYPE_COLORS = {
    "P_biofilm": "#1BAF7A",
    "W_biofilm": "#EB6834",
    "water": "#2A78D6",
}
TYPE_MARKERS = {"P_biofilm": "o", "W_biofilm": "s", "water": "^"}

PHYLUM_COLORS = {
    "Pseudomonadota": "#2F78C9",
    "Thermoplasmatota": "#EF6A32",
    "Nitrospirota": "#21AE83",
    "Bacillota--other": "#EBA500",
    "Actinomycetota": "#DE6FA1",
    "Planctomycetota": "#168A17",
    "Thermoproteota": "#5141A5",
    "Acidobacteriota": "#E34A4A",
    "Unclassified": "#95521B",
    "Others": "#918F8A",
}
NICHE_COLORS = {
    "Generalist": "#1BAF7A",
    "Specialist": "#982568",
    "Unclassified": "#8E8E8E",
    "Insufficient_data": "#D5D2CA",
}

SAMPLE_ORDER = [
    "1st_BFP",
    "1st_BFW1",
    "1st_BFW2",
    "2nd_BFP",
    "2nd_BFW",
    "2nd_wt",
    "3rd_wt",
]


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.0,
        "axes.labelsize": 11.0,
        "axes.labelweight": "bold",
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
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
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def serialize_tsv(rows: list[list[object]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def write_if_changed(path: Path, payload: bytes) -> str:
    if path.exists() and path.read_bytes() == payload:
        return "UNCHANGED"
    status = "UPDATED" if path.exists() else "CREATED"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return status


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_axes(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")
    ax.spines["left"].set_linewidth(1.15)
    ax.spines["bottom"].set_linewidth(1.15)
    if grid:
        ax.grid(axis="y", color="#E3E1DB", linewidth=0.65, zorder=0)
    ax.tick_params(colors="#000000", width=1.0)
    ax.xaxis.label.set_fontsize(11.0)
    ax.yaxis.label.set_fontsize(11.0)
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontweight("bold")
    ax.xaxis.label.set_color("#000000")
    ax.yaxis.label.set_color("#000000")


def save_figure(fig: plt.Figure, fig_dir: Path, stem: str) -> list[Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    settings = [
        ("pdf", None, {}),
        ("svg", None, {}),
        ("png", 300, {}),
        ("tiff", 400, {"pil_kwargs": {"compression": "tiff_lzw"}}),
    ]
    for extension, dpi, extra in settings:
        path = fig_dir / f"{stem}.{extension}"
        fig.savefig(
            path,
            format=extension,
            dpi=dpi,
            facecolor="white",
            edgecolor="none",
            transparent=False,
            bbox_inches=None,
            pad_inches=0,
            metadata={"Creator": "plot_a4_figures_v2.py"} if extension == "pdf" else None,
            **extra,
        )
        outputs.append(path)
    width, height = fig.get_size_inches()
    print(f"SAVED {stem}: {width:.5f} x {height:.5f} in; " + ", ".join(p.name for p in outputs))
    plt.close(fig)
    return outputs


def figure_cpcoa(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "Fig2_CPCoA_source_A4_v1.tsv")
    fig = plt.figure(figsize=A4_LANDSCAPE)
    ax = fig.add_axes([0.065, 0.115, 0.625, 0.81])

    label_offsets = {
        "1st_BFP": (8, 8),
        "2nd_BFP": (8, 7),
        "1st_BFW1": (8, 8),
        "1st_BFW2": (8, 8),
        "2nd_BFW": (8, 7),
        "2nd_wt": (-10, -16),
        "3rd_wt": (-18, 11),
    }
    for type_code in TYPE_ORDER:
        subset = [row for row in rows if row["Type_code"] == type_code]
        xs = [float(row["CAP1"]) for row in subset]
        ys = [float(row["CAP2"]) for row in subset]
        ax.scatter(
            xs,
            ys,
            s=82,
            marker=TYPE_MARKERS[type_code],
            color=TYPE_COLORS[type_code],
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        for row, x, y in zip(subset, xs, ys):
            ax.annotate(
                row["Sample"],
                (x, y),
                xytext=label_offsets[row["Sample"]],
                textcoords="offset points",
                ha="left" if label_offsets[row["Sample"]][0] >= 0 else "right",
                va="center",
                fontsize=8.1,
                color="#303030",
                arrowprops={"arrowstyle": "-", "color": "#8C8C8C", "lw": 0.65},
                zorder=4,
            )

    ax.axhline(0, color="#D4D1C9", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#D4D1C9", linewidth=0.8, zorder=1)
    ax.set_xlim(-0.92, 1.23)
    ax.set_ylim(-1.25, 1.58)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("CAP1", fontweight="bold")
    ax.set_ylabel("CAP2", fontweight="bold")
    clean_axes(ax, grid=False)

    legend_ax = fig.add_axes([0.715, 0.55, 0.27, 0.36])
    legend_ax.axis("off")
    legend_ax.text(0, 1.0, "Type", fontsize=11.0, fontweight="bold", va="top")
    y_positions = [0.74, 0.49, 0.24]
    for type_code, y in zip(TYPE_ORDER, y_positions):
        legend_ax.plot(
            [0.04], [y], marker=TYPE_MARKERS[type_code], markersize=9.8,
            color=TYPE_COLORS[type_code], markeredgecolor="white", markeredgewidth=0.7,
            linestyle="none", transform=legend_ax.transAxes, clip_on=False,
        )
        legend_ax.text(
            0.12, y,
            f"{TYPE_LABELS[type_code]}\n({type_code})" if type_code != "water" else "water",
            transform=legend_ax.transAxes, va="center", fontsize=9.2, linespacing=1.06,
        )
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    return save_figure(fig, fig_dir, "Fig2_CPCoA_A4_v2")


def figure_genus_bubble(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "Fig3_Genus_Bubble_source_A4_v1.tsv")
    taxa = [
        row["Taxon"]
        for row in sorted(rows, key=lambda row: int(row["rank_order"]))
        if row["sample"] == SAMPLE_ORDER[0]
    ]
    status_by_taxon = {row["Taxon"]: row["Taxon_status"] for row in rows}
    samples_by_type = {
        type_code: [s for s in SAMPLE_ORDER if any(r["sample"] == s and r["Type_code"] == type_code for r in rows)]
        for type_code in TYPE_ORDER
    }

    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=[2, 3, 2, 1.50],
        left=0.285,
        right=0.975,
        bottom=0.135,
        top=0.90,
        wspace=0.035,
    )
    y = np.arange(len(taxa))
    title_map = {
        "P_biofilm": "AMD precipitate-associated\nbiofilm (P_biofilm)",
        "W_biofilm": "wall-attached biofilm\n(W_biofilm)",
        "water": "water",
    }
    for panel, type_code in enumerate(TYPE_ORDER):
        ax = fig.add_subplot(grid[0, panel])
        samples = samples_by_type[type_code]
        for x, sample in enumerate(samples):
            sample_rows = {row["Taxon"]: row for row in rows if row["sample"] == sample}
            values = np.array([float(sample_rows[t]["relative_abundance_pct"]) for t in taxa])
            visible = values > 0
            sizes = np.where(visible, np.maximum(4.0, values * 8.5), 0.0)
            ax.scatter(
                np.full(visible.sum(), x), y[visible], s=sizes[visible],
                color=TYPE_COLORS[type_code], edgecolor="white", linewidth=0.45,
                alpha=0.95, zorder=3,
            )
        ax.set_xlim(-0.60, len(samples) - 0.40)
        ax.set_ylim(len(taxa) - 0.30, -0.70)
        ax.set_xticks(range(len(samples)), samples, rotation=43, ha="right", rotation_mode="anchor")
        ax.set_yticks(y)
        if panel == 0:
            ax.set_yticklabels(taxa)
            for tick, taxon in zip(ax.get_yticklabels(), taxa):
                status = status_by_taxon[taxon]
                tick.set_fontsize(7.6)
                if status == "Assigned genus":
                    tick.set_fontstyle("italic")
                    tick.set_color("#303030")
                elif taxon == "Others":
                    tick.set_color("#555555")
                else:
                    tick.set_color("#8A4716")
        else:
            ax.set_yticklabels([])
        ax.set_title(title_map[type_code], fontsize=10.0, fontweight="bold", pad=8, linespacing=1.05)
        ax.grid(axis="y", color="#E0DED8", linewidth=0.55, zorder=0)
        ax.spines["left"].set_color("#000000")
        ax.spines["bottom"].set_color("#000000")
        ax.spines["left"].set_linewidth(1.05)
        ax.spines["bottom"].set_linewidth(1.05)
        ax.spines["top"].set_color("#D2CFC8")
        ax.spines["right"].set_color("#D2CFC8")
        ax.spines["top"].set_linewidth(0.55)
        ax.spines["right"].set_linewidth(0.55)
        ax.tick_params(axis="x", labelsize=7.5, pad=2)
        ax.tick_params(axis="y", length=3)

    key_ax = fig.add_subplot(grid[0, 3])
    key_ax.axis("off")
    key_ax.text(0.03, 0.98, "Relative abundance\n(%)", fontsize=9.8, fontweight="bold", va="top", linespacing=1.02)
    for value, y_key in zip([5, 20, 40], [0.80, 0.65, 0.47]):
        key_ax.scatter([0.20], [y_key], s=max(4.0, value * 8.5), color="#5A5955", edgecolor="white", linewidth=0.45)
        key_ax.text(0.47, y_key, str(value), va="center", fontsize=8.7)
    key_ax.text(
        0.03,
        0.25,
        "Brown labels:\ngenus unclassified;\nnearest assigned\nhigher rank shown",
        fontsize=8.2,
        color="#8A4716",
        va="top",
        linespacing=1.12,
    )
    key_ax.set_xlim(0, 1)
    key_ax.set_ylim(0, 1)
    return save_figure(fig, fig_dir, "Fig3_Genus_Bubble_A4_v2")


def figure_ferrovum(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "Fig4_Ferrovum_top15_source_A4_v1.tsv")
    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(
        1, 4, width_ratios=[0.95, 3.5, 2.45, 3.15],
        left=0.035, right=0.985, bottom=0.12, top=0.90, wspace=0.06,
    )
    axes = [fig.add_subplot(grid[0, i]) for i in range(4)]
    y = np.arange(len(rows))
    for ax in axes:


        ax.set_ylim(len(rows) - 0.45, -1.35)
        for separator in np.arange(0.5, len(rows), 1):
            ax.axhline(separator, color="#ECEAE5", linewidth=0.55, zorder=0)

    axes[0].axis("off")
    axes[0].set_xlim(0, 1)
    axes[0].text(0.0, -0.90, "ASV", fontsize=10.4, fontweight="bold", transform=axes[0].transData)
    for index, row in enumerate(rows):
        label = row["ASV_ID"] + (" *" if row["ASV_ID"] == "ASV_1" else "")
        axes[0].text(0.0, index, label, va="center", fontsize=8.1, fontweight="bold", color="#66645F")

    axes[1].axis("off")
    axes[1].set_xlim(0, 1)
    axes[1].text(0.0, -0.90, "Genus-level label", fontsize=10.4, fontweight="bold", transform=axes[1].transData)
    for index, row in enumerate(rows):
        assigned = row["label_status"] == "Assigned genus"
        axes[1].text(
            0.0, index, row["display_label"], va="center", fontsize=8.1,
            fontstyle="italic" if assigned else "normal",
            color="#252525" if assigned else "#8A4716",
        )

    abundance = np.array([float(row["pct_of_total"]) for row in rows])
    axes[2].barh(y, abundance, height=0.62, color="#5A5955", edgecolor="white", linewidth=0.4, zorder=2)
    axes[2].set_xlim(0, 25.5)
    axes[2].set_xlabel("Pooled relative abundance (%)", fontweight="bold")
    axes[2].set_title("Abundance", fontsize=10.4, fontweight="bold", pad=9)
    axes[2].set_yticks([])
    for index, value in enumerate(abundance):
        axes[2].text(value + 0.35, index, f"{value:.1f}%", va="center", fontsize=7.7, color="#454545")
    clean_axes(axes[2], grid=False)
    axes[2].grid(axis="x", color="#E0DED8", linewidth=0.6, zorder=0)

    bootstrap = np.array([float(row["Genus_boot"]) for row in rows])
    meets = np.array([row["meets_minBoot80"].lower() == "true" for row in rows])
    bar_colors = np.where(meets, "#2A64AD", "#98501C")
    axes[3].barh(y, bootstrap, height=0.62, color=bar_colors, edgecolor="white", linewidth=0.4, zorder=2)
    axes[3].axvline(80, color="#8A4716", linestyle=(0, (3, 2)), linewidth=1.2, zorder=3)
    axes[3].text(
        82.0, -0.83, "minBoot = 80",
        transform=axes[3].transData, ha="left", va="center",
        fontsize=8.8, fontweight="bold", color="#8A4716",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.94},
        zorder=5,
    )
    axes[3].set_xlim(0, 118)
    axes[3].set_xlabel("Genus-level bootstrap support (%)", fontweight="bold")
    axes[3].set_title("Taxonomic support", fontsize=10.4, fontweight="bold", pad=9)
    axes[3].set_yticks([])
    for index, value in enumerate(bootstrap):
        axes[3].text(value + 2.2, index, f"{value:.0f}", va="center", fontsize=7.7, color=bar_colors[index], fontweight="bold")
    clean_axes(axes[3], grid=False)
    axes[3].grid(axis="x", color="#E0DED8", linewidth=0.6, zorder=0)
    return save_figure(fig, fig_dir, "SFig11_Ferrovum_Taxonomic_Resolution_A4_v2")


def figure_alpha(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "SFig1_AlphaDiversity_source_A4_v1.tsv")
    metrics = [
        ("Observed", "Observed ASVs"),
        ("Chao1", "Chao1"),
        ("ACE", "ACE"),
        ("Shannon", "Shannon (H')"),
        ("Simpson", "Simpson (1-D)"),
        ("InvSimpson", "Inverse Simpson"),
        ("Pielou_evenness", "Pielou's evenness (J')"),
    ]
    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(2, 4, left=0.055, right=0.985, bottom=0.115, top=0.95, wspace=0.30, hspace=0.42)
    count_by_type = Counter(row["Type_code"] for row in rows)
    for panel, (column, title) in enumerate(metrics):
        ax = fig.add_subplot(grid[panel // 4, panel % 4])
        for x, type_code in enumerate(TYPE_ORDER):
            subset = [row for row in rows if row["Type_code"] == type_code]
            values = np.array([float(row[column]) for row in subset])
            offsets = np.linspace(-0.08, 0.08, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(
                x + offsets, values, s=42, marker=TYPE_MARKERS[type_code],
                color=TYPE_COLORS[type_code], edgecolor="white", linewidth=0.6, zorder=3,
            )
            median = float(np.median(values))
            ax.plot([x - 0.16, x + 0.16], [median, median], color=TYPE_COLORS[type_code], linewidth=2.2, solid_capstyle="round", zorder=4)
        ax.set_xlim(-0.52, 2.52)
        ax.set_xticks(
            range(3),
            [f"{t}\n(n={count_by_type[t]})" for t in TYPE_ORDER],
            fontsize=7.5,
        )
        ax.set_title(title, fontsize=10.2, fontweight="bold", pad=7)
        clean_axes(ax, grid=True)

    legend_ax = fig.add_subplot(grid[1, 3])
    legend_ax.axis("off")
    legend_ax.text(0.0, 0.98, "Type", fontsize=10.4, fontweight="bold", va="top")
    for type_code, y in zip(TYPE_ORDER, [0.74, 0.46, 0.20]):
        legend_ax.plot(
            [0.08], [y], marker=TYPE_MARKERS[type_code], markersize=8.5,
            color=TYPE_COLORS[type_code], markeredgecolor="white", linestyle="none",
            transform=legend_ax.transAxes,
        )
        legend_ax.text(
            0.20, y,
            f"{type_code}\n{TYPE_LABELS[type_code]}" if type_code != "water" else "water",
            transform=legend_ax.transAxes, va="center", fontsize=8.6, linespacing=1.06,
        )
    return save_figure(fig, fig_dir, "SFig1_AlphaDiversity_A4_v2")


def figure_bray_curtis(root: Path, source_dir: Path, fig_dir: Path) -> list[Path]:
    seqtab_path = root / "output_data" / "seqtab_final_long.csv"
    sample_sheet_path = root / "input_data" / "sample_sheet.tsv"
    seq_rows = read_csv(seqtab_path)
    sample_meta_rows = read_tsv(sample_sheet_path)
    sample_meta = {row["analysis_name"]: row for row in sample_meta_rows}
    matrix = np.array([[float(row[sample]) for row in seq_rows] for sample in SAMPLE_ORDER])
    distances = squareform(pdist(matrix, metric="braycurtis"))
    condensed = squareform(distances, checks=True)
    cluster = linkage(condensed, method="average")
    order_index = leaves_list(cluster).tolist()
    ordered_samples = [SAMPLE_ORDER[index] for index in order_index]
    ordered = distances[np.ix_(order_index, order_index)]

    expected_order = ["2nd_BFW", "2nd_wt", "3rd_wt", "2nd_BFP", "1st_BFW1", "1st_BFP", "1st_BFW2"]
    if ordered_samples != expected_order:
        raise ValueError(f"Unexpected Bray-Curtis display order: {ordered_samples}")
    matrix_rows: list[list[object]] = [["sample"] + SAMPLE_ORDER]
    for sample, row in zip(SAMPLE_ORDER, distances):
        matrix_rows.append([sample] + [f"{value:.12f}" for value in row])
    matrix_path = source_dir / "SFig2_BrayCurtis_matrix_A4_v2.tsv"
    matrix_status = write_if_changed(matrix_path, serialize_tsv(matrix_rows))
    order_rows = [["display_order", "sample", "Type_code"]] + [
        [index + 1, sample, sample_meta[sample]["Type"]] for index, sample in enumerate(ordered_samples)
    ]
    order_path = source_dir / "SFig2_BrayCurtis_display_order_A4_v2.tsv"
    order_status = write_if_changed(order_path, serialize_tsv(order_rows))
    print(f"HEATMAP SOURCE: matrix={matrix_status}, order={order_status}; order={ordered_samples}")

    fig = plt.figure(figsize=A4_LANDSCAPE)


    heat_left, heat_bottom, heat_height = 0.305, 0.145, 0.615
    heat_width = heat_height * (210.0 / 297.0)
    top_dend_ax = fig.add_axes([heat_left, 0.815, heat_width, 0.12])
    dendrogram(cluster, ax=top_dend_ax, color_threshold=0, above_threshold_color="#4A4946", no_labels=True)
    top_dend_ax.axis("off")
    top_strip_ax = fig.add_axes([heat_left, 0.775, heat_width, 0.030])
    top_rgba = np.array([[to_rgba(TYPE_COLORS[sample_meta[s]["Type"]]) for s in ordered_samples]])
    top_strip_ax.imshow(top_rgba, aspect="auto", interpolation="nearest", extent=(-0.5, 6.5, -0.5, 0.5))
    top_strip_ax.set_xlim(-0.5, 6.5)
    top_strip_ax.set_ylim(-0.5, 0.5)
    top_strip_ax.set_xticks([])
    top_strip_ax.set_yticks([])
    for boundary in np.arange(0.5, 6.5, 1.0):
        top_strip_ax.axvline(boundary, color="white", linewidth=0.7)
    for spine in top_strip_ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(0.85)

    left_dend_ax = fig.add_axes([0.135, heat_bottom, 0.120, heat_height])
    dendrogram(cluster, ax=left_dend_ax, orientation="left", color_threshold=0, above_threshold_color="#4A4946", no_labels=True)
    left_dend_ax.axis("off")
    left_strip_ax = fig.add_axes([0.270, heat_bottom, 0.023, heat_height])
    left_rgba = np.array([[to_rgba(TYPE_COLORS[sample_meta[s]["Type"]])] for s in ordered_samples])
    left_strip_ax.imshow(left_rgba, aspect="auto", interpolation="nearest", extent=(-0.5, 0.5, 6.5, -0.5))
    left_strip_ax.set_xlim(-0.5, 0.5)
    left_strip_ax.set_ylim(6.5, -0.5)
    left_strip_ax.set_xticks([])
    left_strip_ax.set_yticks([])
    for boundary in np.arange(0.5, 6.5, 1.0):
        left_strip_ax.axhline(boundary, color="white", linewidth=0.7)
    for spine in left_strip_ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(0.85)

    heat_ax = fig.add_axes([heat_left, heat_bottom, heat_width, heat_height])
    cmap = LinearSegmentedColormap.from_list("bray", ["#F7F7F5", "#B8D0EA", "#4A8BDA", "#174A7E"])
    image = heat_ax.imshow(ordered, vmin=0, vmax=1, cmap=cmap, aspect="equal", interpolation="nearest")
    heat_ax.set_xticks(range(7), ordered_samples, rotation=43, ha="right", rotation_mode="anchor")
    heat_ax.set_yticks(range(7), ordered_samples)
    heat_ax.yaxis.tick_right()
    heat_ax.tick_params(axis="x", colors="#000000", width=1.0, labelsize=8.2)
    heat_ax.tick_params(axis="y", colors="#000000", width=1.0, labelsize=8.2, pad=5)
    for spine in heat_ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(1.05)
    for i in range(7):
        for j in range(7):
            value = ordered[i, j]
            heat_ax.text(
                j, i, f"{value:.2f}", ha="center", va="center", fontsize=7.3,
                color="white" if value >= 0.57 else "#3F3F3F",
            )

    legend_ax = fig.add_axes([heat_left + heat_width + 0.060, 0.58, 0.19, 0.32])
    legend_ax.axis("off")
    legend_ax.text(0.0, 0.99, "Type", fontsize=10.5, fontweight="bold", va="top")
    for type_code, y in zip(TYPE_ORDER, [0.74, 0.47, 0.20]):
        legend_ax.add_patch(Rectangle((0.01, y - 0.06), 0.09, 0.12, transform=legend_ax.transAxes, color=TYPE_COLORS[type_code], clip_on=False))
        legend_ax.text(
            0.15, y,
            f"{type_code}\n{TYPE_LABELS[type_code]}" if type_code != "water" else "water",
            transform=legend_ax.transAxes, va="center", fontsize=8.6, linespacing=1.05,
        )

    color_ax = fig.add_axes([0.855, 0.17, 0.024, 0.27])
    cbar = fig.colorbar(image, cax=color_ax)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.tick_params(labelsize=8.2, width=1.0, colors="#000000")
    cbar.set_label("Bray-Curtis dissimilarity", fontsize=10.0, fontweight="bold", labelpad=7)
    cbar.outline.set_edgecolor("#000000")
    cbar.outline.set_linewidth(1.0)
    return save_figure(fig, fig_dir, "SFig2_BrayCurtis_Heatmap_A4_v2")


def figure_phylum(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "SFig3_Phylum_Composition_display_A4_v1.tsv")
    phyla = [
        row["Phylum"]
        for row in sorted(rows, key=lambda row: int(row["rank_order"]))
        if row["sample"] == SAMPLE_ORDER[0]
    ]
    samples_by_type = {
        t: [s for s in SAMPLE_ORDER if any(r["sample"] == s and r["Type_code"] == t for r in rows)]
        for t in TYPE_ORDER
    }
    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(
        1, 4, width_ratios=[2, 3, 2, 1.35], left=0.07, right=0.99, bottom=0.145, top=0.90, wspace=0.05
    )
    title_map = {
        "P_biofilm": "AMD precipitate-associated\nbiofilm (P_biofilm)",
        "W_biofilm": "wall-attached biofilm\n(W_biofilm)",
        "water": "water",
    }
    for panel, type_code in enumerate(TYPE_ORDER):
        ax = fig.add_subplot(grid[0, panel])
        samples = samples_by_type[type_code]
        bottom = np.zeros(len(samples))
        for phylum in phyla:
            values = np.array([
                float(next(r["relative_abundance_pct"] for r in rows if r["sample"] == sample and r["Phylum"] == phylum))
                for sample in samples
            ])
            ax.bar(
                range(len(samples)), values, bottom=bottom, width=0.72,
                color=PHYLUM_COLORS[phylum], edgecolor="white", linewidth=0.45, label=phylum, zorder=2,
            )
            bottom += values
        if not np.allclose(bottom, 100.0, atol=1e-7):
            raise ValueError(f"Phylum bars do not sum to 100 for {type_code}: {bottom}")
        ax.set_ylim(0, 100)
        ax.set_xlim(-0.58, len(samples) - 0.42)
        ax.set_xticks(range(len(samples)), samples, rotation=42, ha="right", rotation_mode="anchor")
        ax.set_title(title_map[type_code], fontsize=10.0, fontweight="bold", pad=8, linespacing=1.05)
        if panel == 0:
            ax.set_ylabel("Relative abundance (%)", fontweight="bold")
        else:
            ax.set_yticklabels([])
        clean_axes(ax, grid=True)
        ax.set_axisbelow(True)

    legend_ax = fig.add_subplot(grid[0, 3])
    legend_ax.axis("off")
    legend_ax.text(0.02, 0.98, "Phylum", fontsize=10.4, fontweight="bold", va="top")
    for phylum, y in zip(phyla, np.linspace(0.87, 0.10, len(phyla))):
        legend_ax.add_patch(Rectangle((0.02, y - 0.018), 0.08, 0.036, color=PHYLUM_COLORS[phylum], transform=legend_ax.transAxes, clip_on=False))
        legend_ax.text(0.15, y, phylum, va="center", fontsize=8.4, transform=legend_ax.transAxes)
    return save_figure(fig, fig_dir, "SFig3_Phylum_Composition_A4_v2")


def figure_niche_counts(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "TableS5_Niche_Breadth_primary_A4_v1.tsv")
    classes = ["Generalist", "Specialist", "Unclassified", "Insufficient_data"]
    all_counts = Counter(row["classification"] for row in rows)
    eligible_rows = [row for row in rows if row["eligible (baseline gate)"].lower() == "true"]
    eligible_counts = Counter(row["classification"] for row in eligible_rows)
    panels = [
        (all_counts, "All genera", "n = 129"),
        (eligible_counts, "Eligible genera (baseline gate)", "n = 44; present in >=2/7 samples and mean >=0.01%"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=A4_LANDSCAPE, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.86, wspace=0.10)
    for ax, (counts, title, subtitle) in zip(axes, panels):
        values = [counts.get(cls, 0) for cls in classes]
        bars = ax.bar(
            range(4), values, width=0.62, color=[NICHE_COLORS[cls] for cls in classes],
            edgecolor="white", linewidth=0.6, zorder=2,
        )
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 2.0, str(value), ha="center", va="bottom", fontsize=9.0, fontweight="bold", color="#4A4946")
        ax.set_xticks(range(4), ["Generalist", "Nominal\nspecialist", "Unclassified", "Insufficient\ndata"])
        ax.set_ylim(0, 104)
        ax.set_title(f"{title}\n{subtitle}", fontsize=10.1, fontweight="bold", pad=9, linespacing=1.15)
        clean_axes(ax, grid=True)
    axes[0].set_ylabel("Number of genera", fontweight="bold")
    fig.text(
        0.985, 0.025, "Nominal specialist candidates are sensitivity-unstable; stable calls = 0.",
        ha="right", va="bottom", fontsize=8.4, fontweight="bold", color="#70204E",
    )
    return save_figure(fig, fig_dir, "SFig4_NicheBreadth_Counts_A4_v2")


def figure_niche_b(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "TableS5_Niche_Breadth_primary_A4_v1.tsv")
    eligible = [row for row in rows if row["eligible (baseline gate)"].lower() == "true"]
    class_order = ["Specialist", "Unclassified"]
    values_by_class = {
        cls: np.array([float(row["Levins_B_standardized"]) for row in eligible if row["classification"] == cls])
        for cls in class_order
    }
    if [len(values_by_class[cls]) for cls in class_order] != [2, 42]:
        raise ValueError("Unexpected eligible niche class counts")

    fig = plt.figure(figsize=A4_PORTRAIT)
    ax = fig.add_axes([0.145, 0.19, 0.81, 0.67])
    box = ax.boxplot(
        [values_by_class[cls] for cls in class_order], positions=[0, 1], widths=0.46,
        patch_artist=True, showfliers=False,
        boxprops={"linewidth": 1.0, "edgecolor": "#4D4C49"},
        whiskerprops={"linewidth": 1.0, "color": "#4D4C49"},
        capprops={"linewidth": 1.0, "color": "#4D4C49"},
        medianprops={"linewidth": 2.0, "color": "white"},
    )
    for patch, cls in zip(box["boxes"], class_order):
        patch.set_facecolor(NICHE_COLORS[cls])
        patch.set_alpha(0.90)
    rng = np.random.default_rng(100)
    for x, cls in enumerate(class_order):
        values = values_by_class[cls]
        jitter = rng.uniform(-0.095, 0.095, len(values))
        ax.scatter(
            x + jitter, values, s=30, color=NICHE_COLORS[cls],
            edgecolor="white", linewidth=0.55, alpha=0.85, zorder=3,
        )
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 1], ["Nominal specialist\n(n=2; unstable)", "Unclassified\n(n=42)"], fontsize=9.2)
    ax.set_ylabel("Levins' B (standardized; three Type states)", fontweight="bold")
    ax.set_title("Nominal, sensitivity-unstable classification among eligible genera (n = 44)", fontsize=10.8, fontweight="bold", pad=12)
    clean_axes(ax, grid=True)
    return save_figure(fig, fig_dir, "SFig5_NicheBreadth_LevinsB_A4_v2")


def figure_venn(source_dir: Path, fig_dir: Path) -> list[Path]:
    rows = read_tsv(source_dir / "SFig10_Venn_regions_A4_v1.tsv")
    values = {row["region"]: int(row["ASV_count"]) for row in rows}
    fig = plt.figure(figsize=A4_PORTRAIT)
    ax = fig.add_axes([0.055, 0.12, 0.89, 0.78])
    ax.set_aspect("equal")
    ax.axis("off")
    radius = 1.45
    centers = {"P_biofilm": (-0.72, 0.35), "W_biofilm": (0.72, 0.35), "water": (0.0, -0.70)}
    for type_code in TYPE_ORDER:
        ax.add_patch(
            Circle(
                centers[type_code], radius,
                facecolor=TYPE_COLORS[type_code], edgecolor=TYPE_COLORS[type_code],
                linewidth=2.2, alpha=0.46,
            )
        )

    number_positions = {
        "P_biofilm_only": (-1.30, 0.50),
        "W_biofilm_only": (1.30, 0.50),
        "water_only": (0.0, -1.62),
        "P_biofilm_and_W_biofilm_only": (0.0, 1.18),
        "P_biofilm_and_water_only": (-0.62, -0.38),
        "W_biofilm_and_water_only": (0.62, -0.38),
        "all_three_Types": (0.0, 0.22),
    }
    for region, position in number_positions.items():
        ax.text(*position, str(values[region]), ha="center", va="center", fontsize=13.0, color="#111111")

    ax.text(
        -1.58, 1.95, "P_biofilm\nAMD precipitate-associated biofilm",
        ha="center", va="bottom", color=TYPE_COLORS["P_biofilm"], fontsize=9.0, fontweight="bold", linespacing=1.05,
    )
    ax.text(
        1.58, 1.95, "W_biofilm\nwall-attached biofilm",
        ha="center", va="bottom", color=TYPE_COLORS["W_biofilm"], fontsize=9.0, fontweight="bold", linespacing=1.05,
    )
    ax.text(0.0, -2.30, "water", ha="center", va="top", color=TYPE_COLORS["water"], fontsize=10.0, fontweight="bold")
    ax.set_xlim(-2.45, 2.45)
    ax.set_ylim(-2.45, 2.28)
    return save_figure(fig, fig_dir, "SFig10_Venn_ASV_byType_A4_v2")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path("/input/youre/data/"))
    parser.add_argument("--output-root", type=Path, default=Path("/output/your/data"))
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    source_dir = source_root / "source_data"
    v2_source_dir = source_root / "source_data"
    v2_source_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_root / "build" / "non_tree_stage" / "figures"
    required_sources = [
        "Fig2_CPCoA_source_A4_v1.tsv",
        "Fig3_Genus_Bubble_source_A4_v1.tsv",
        "Fig4_Ferrovum_top15_source_A4_v1.tsv",
        "SFig1_AlphaDiversity_source_A4_v1.tsv",
        "SFig3_Phylum_Composition_display_A4_v1.tsv",
    ]
    missing = [name for name in required_sources if not (source_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing prepared source files: {missing}")

    outputs: list[Path] = []
    outputs += figure_cpcoa(source_dir, fig_dir)
    outputs += figure_genus_bubble(source_dir, fig_dir)
    outputs += figure_ferrovum(source_dir, fig_dir)
    outputs += figure_alpha(source_dir, fig_dir)
    outputs += figure_phylum(source_dir, fig_dir)

    register_rows: list[list[object]] = [["output", "bytes", "sha256"]]
    for path in outputs:
        register_rows.append([str(path.relative_to(output_root)).replace("\\", "/"), path.stat().st_size, sha256(path)])
    register_path = output_root / "logs" / "03_python_active_figure_stage_v1.tsv"
    status = write_if_changed(register_path, serialize_tsv(register_rows))
    print(f"OUTPUT REGISTER: {status}; {len(outputs)} files; {register_path}")
    print("PYTHON FIGURE GENERATION: PASS")


if __name__ == "__main__":
    main()
