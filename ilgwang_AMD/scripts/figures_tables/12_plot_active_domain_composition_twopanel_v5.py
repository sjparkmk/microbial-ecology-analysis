\
\
\
\
\

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from PIL import Image

MM = 1.0 / 25.4
A4_LANDSCAPE = (297 * MM, 210 * MM)
STEM = "SFig3_Domain_Read_Fractions_A4_v5"
SAMPLES = ["1st_BFP", "1st_BFW1", "1st_BFW2", "2nd_BFP", "2nd_BFW", "2nd_wt", "3rd_wt"]
DOMAINS = ["Bacteria", "Archaea"]
DOMAIN_COLORS = {"Bacteria": "#A8ACB0", "Archaea": "#5B3A86"}
TYPE_ORDER = ["P_biofilm", "W_biofilm", "water"]
TYPE_COLORS = {"P_biofilm": "#1BAF7A", "W_biofilm": "#EB6834", "water": "#2A78D6"}
TYPE_MARKERS = {"P_biofilm": "o", "W_biofilm": "s", "water": "^"}
TYPE_LABELS = {
    "P_biofilm": "P_biofilm (AMD precipitate-associated biofilm)",
    "W_biofilm": "W_biofilm (wall-attached biofilm)",
    "water": "water",
}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.5, "axes.labelsize": 10.5, "axes.labelweight": "bold",
    "axes.titlesize": 10.8, "axes.titleweight": "bold", "xtick.labelsize": 8.2, "ytick.labelsize": 8.2,
    "axes.edgecolor": "#000000", "axes.linewidth": 1.1, "xtick.color": "#000000", "ytick.color": "#000000",
    "xtick.major.width": 0.9, "ytick.major.width": 0.9, "xtick.major.size": 3.2, "ytick.major.size": 3.2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "figure.facecolor": "white",
    "axes.facecolor": "white", "savefig.facecolor": "white",
})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def check(checks: list[dict[str, str]], name: str, passed: bool, observed: object, expected: object) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": str(observed), "expected": str(expected)})
    if not passed:
        raise ValueError(f"{name}: observed={observed!r}; expected={expected!r}")


def build_data(source_root: Path) -> tuple[list[dict[str, object]], list[dict[str, str]], list[Path]]:
    data_root = source_root
    tax_path = data_root / "output_data" / "taxa_silva144.csv"
    seq_path = data_root / "output_data" / "seqtab_final_long.csv"
    meta_path = data_root / "input_data" / "sample_sheet.tsv"
    taxonomy_rows = read_csv(tax_path)
    seq_rows = read_csv(seq_path)
    metadata_rows = read_csv(meta_path, "\t")
    checks: list[dict[str, str]] = []
    check(checks, "canonical_taxonomy_rows", len(taxonomy_rows) == 663, len(taxonomy_rows), 663)
    check(checks, "canonical_seqtab_rows", len(seq_rows) == 663, len(seq_rows), 663)
    check(checks, "metadata_rows", len(metadata_rows) == 7, len(metadata_rows), 7)
    tax_map = {row["sequence"]: row.get("Kingdom", "").strip() for row in taxonomy_rows}
    check(checks, "taxonomy_unique_sequence", len(tax_map) == 663, len(tax_map), 663)
    seq_sequences = [row["sequence"] for row in seq_rows]
    check(checks, "seqtab_unique_sequence", len(set(seq_sequences)) == 663, len(set(seq_sequences)), 663)
    check(checks, "exact_sequence_join", set(seq_sequences) == set(tax_map), len(set(seq_sequences) ^ set(tax_map)), 0)
    meta = {row["analysis_name"]: row for row in metadata_rows}
    check(checks, "sample_set", set(meta) == set(SAMPLES), sorted(meta), SAMPLES)

    counts = {sample: {domain: 0 for domain in DOMAINS} for sample in SAMPLES}
    for row in seq_rows:
        kingdom = tax_map[row["sequence"]]
        check(checks, f"kingdom_{row['ASV_ID']}", kingdom in DOMAINS, kingdom, DOMAINS)
        for sample in SAMPLES:
            counts[sample][kingdom] += int(row[sample])
    bacteria_total = sum(counts[s]["Bacteria"] for s in SAMPLES)
    archaea_total = sum(counts[s]["Archaea"] for s in SAMPLES)
    check(checks, "bacterial_total", bacteria_total == 610601, bacteria_total, 610601)
    check(checks, "archaeal_total", archaea_total == 106314, archaea_total, 106314)
    check(checks, "total_reads", bacteria_total + archaea_total == 716915, bacteria_total + archaea_total, 716915)

    output: list[dict[str, object]] = []
    for idx, sample in enumerate(SAMPLES, start=1):
        b = counts[sample]["Bacteria"]
        a = counts[sample]["Archaea"]
        total = a + b
        row = {
            "sample_order": idx, "sample": sample, "Group": meta[sample]["Group"], "Type_code": meta[sample]["Type"],
            "Type_label": meta[sample]["Type_label"], "Archaea_reads": a, "Bacteria_reads": b, "total_reads": total,
            "within_sample_archaea_pct": 100.0 * a / total, "within_sample_bacteria_pct": 100.0 * b / total,
            "archaea_share_across_samples_pct": 100.0 * a / archaea_total,
            "bacteria_share_across_samples_pct": 100.0 * b / bacteria_total,
            "archaea_total_reads": archaea_total, "bacteria_total_reads": bacteria_total,
        }
        output.append(row)
    check(checks, "panel_a_sums", all(abs(r["within_sample_archaea_pct"] + r["within_sample_bacteria_pct"] - 100.0) < 1e-10 for r in output), "all seven sums", "100%")
    check(checks, "panel_b_archaea_sum", math.isclose(sum(r["archaea_share_across_samples_pct"] for r in output), 100.0, abs_tol=1e-10), sum(r["archaea_share_across_samples_pct"] for r in output), "100%")
    check(checks, "panel_b_bacteria_sum", math.isclose(sum(r["bacteria_share_across_samples_pct"] for r in output), 100.0, abs_tol=1e-10), sum(r["bacteria_share_across_samples_pct"] for r in output), "100%")
    check(checks, "type_codes", {r["Type_code"] for r in output} == set(TYPE_ORDER), sorted({r["Type_code"] for r in output}), TYPE_ORDER)
    return output, checks, [tax_path, seq_path, meta_path]


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#000000")
    ax.spines["bottom"].set_color("#000000")
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.tick_params(colors="#000000", width=0.9)
    ax.grid(axis="y", color="#E4E7E5", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)


def add_type_markers(ax: plt.Axes, data: list[dict[str, object]]) -> None:
    for index, row in enumerate(data):
        code = str(row["Type_code"])
        ax.plot([index], [-0.035], marker=TYPE_MARKERS[code], markersize=5.3, markerfacecolor=TYPE_COLORS[code],
                markeredgecolor="#FFFFFF", markeredgewidth=0.45, linestyle="none", transform=ax.get_xaxis_transform(),
                clip_on=False, zorder=5)
    for label, row in zip(ax.get_xticklabels(), data):
        label.set_color(TYPE_COLORS[str(row["Type_code"])])
        label.set_fontweight("bold")


def plot(data: list[dict[str, object]]) -> plt.Figure:
    x = np.arange(len(data))
    fig = plt.figure(figsize=A4_LANDSCAPE)
    grid = fig.add_gridspec(1, 2, left=0.070, right=0.985, bottom=0.245, top=0.770, wspace=0.215)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    bacteria = np.array([float(row["within_sample_bacteria_pct"]) for row in data])
    archaea = np.array([float(row["within_sample_archaea_pct"]) for row in data])
    ax_a.bar(x, bacteria, width=0.68, color=DOMAIN_COLORS["Bacteria"], edgecolor="#FFFFFF", linewidth=0.65, zorder=2)
    ax_a.bar(x, archaea, bottom=bacteria, width=0.68, color=DOMAIN_COLORS["Archaea"], edgecolor="#FFFFFF", linewidth=0.65, zorder=2)
    for index, value in enumerate(archaea):
        ax_a.text(index, 102.1, f"{value:.2f}%", ha="center", va="bottom", fontsize=7.8, color=DOMAIN_COLORS["Archaea"], fontweight="bold")
    ax_a.set_ylim(0, 108)
    ax_a.set_yticks([0, 20, 40, 60, 80, 100])
    ax_a.set_ylabel("Reads within sample (%)", fontweight="bold")
    ax_a.set_xlabel("Sample", fontweight="bold", labelpad=10)
    ax_a.set_xticks(x, [str(row["sample"]) for row in data], rotation=34, ha="right", rotation_mode="anchor")
    ax_a.tick_params(axis="x", pad=11)
    ax_a.set_title("Within-sample domain composition\nLabels show Archaea (%)", pad=10, linespacing=1.25)
    clean_axes(ax_a); add_type_markers(ax_a, data)

    width = 0.34
    bact_share = np.array([float(row["bacteria_share_across_samples_pct"]) for row in data])
    arch_share = np.array([float(row["archaea_share_across_samples_pct"]) for row in data])
    bars_b = ax_b.bar(x - width / 2, bact_share, width=width, color=DOMAIN_COLORS["Bacteria"], edgecolor="#FFFFFF", linewidth=0.55, zorder=2)
    bars_a = ax_b.bar(x + width / 2, arch_share, width=width, color=DOMAIN_COLORS["Archaea"], edgecolor="#FFFFFF", linewidth=0.55, zorder=2)
    ymax = max(55.0, math.ceil((max(float(arch_share.max()), float(bact_share.max())) + 5.0) / 5.0) * 5.0)
    ax_b.set_ylim(0, ymax)
    ax_b.set_ylabel("Share of domain-assigned reads (%)", fontweight="bold")
    ax_b.set_xlabel("Sample", fontweight="bold", labelpad=10)
    ax_b.set_xticks(x, [str(row["sample"]) for row in data], rotation=34, ha="right", rotation_mode="anchor")
    ax_b.tick_params(axis="x", pad=11)
    ax_b.set_title("Distribution of each domain among samples\nEach domain sums to 100% across samples", pad=10, linespacing=1.25)
    clean_axes(ax_b); add_type_markers(ax_b, data)
    for bars, values, color in [(bars_b, bact_share, "#5F6468"), (bars_a, arch_share, DOMAIN_COLORS["Archaea"])]:
        for bar, value in zip(bars, values):
            ax_b.text(bar.get_x() + bar.get_width() / 2, value + 0.85, f"{value:.1f}", ha="center", va="bottom", fontsize=7.6, color=color, fontweight="bold")

    ax_a.text(-0.125, 1.10, "(a)", transform=ax_a.transAxes, ha="left", va="top", fontsize=10.5, fontweight="bold")
    ax_b.text(-0.125, 1.10, "(b)", transform=ax_b.transAxes, ha="left", va="top", fontsize=10.5, fontweight="bold")
    fig.suptitle("Bacterial and archaeal read fractions across Ilgwang AMD samples", x=0.5, y=0.965, fontsize=12.2, fontweight="bold")
    domains = [Patch(facecolor=DOMAIN_COLORS[n], edgecolor="none", label=n) for n in DOMAINS]
    dleg = fig.legend(handles=domains, loc="upper center", bbox_to_anchor=(0.5, 0.912), ncol=2, frameon=False, title="Domain", fontsize=9.0, title_fontsize=9.3, handlelength=1.35, handletextpad=0.45, columnspacing=1.25)
    dleg.get_title().set_fontweight("bold")
    types = [Line2D([0], [0], marker=TYPE_MARKERS[c], linestyle="none", markerfacecolor=TYPE_COLORS[c], markeredgecolor="#FFFFFF", markeredgewidth=0.5, markersize=7.2, label=TYPE_LABELS[c]) for c in TYPE_ORDER]
    tleg = fig.legend(handles=types, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=3, frameon=False, title="Microhabitat Type", fontsize=8.25, title_fontsize=9.2, handletextpad=0.50, columnspacing=1.05)
    tleg.get_title().set_fontweight("bold")
    return fig


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path("/input/youre/data/"))
    parser.add_argument("--output-root", type=Path, default=Path("/output/your/data"))
    args = parser.parse_args()
    source_root, output_root = args.source_root.resolve(), args.output_root.resolve()
    figure_dir = output_root / "build" / "domain_stage" / "figures"
    source_dir, log_dir, qa_dir = output_root / "source_data", output_root / "logs", output_root / "qa"
    for directory in [figure_dir, source_dir, log_dir, qa_dir]: directory.mkdir(parents=True, exist_ok=True)
    paths = {
        "pdf": figure_dir / f"{STEM}.pdf", "svg": figure_dir / f"{STEM}.svg", "png": figure_dir / f"{STEM}.png", "tiff": figure_dir / f"{STEM}.tiff",
        "source": source_dir / f"{STEM}_source.tsv", "validation": log_dir / f"{STEM}_validation.tsv", "geometry": log_dir / f"{STEM}_geometry.tsv",
        "runlog": log_dir / f"{STEM}_run.log", "manifest": log_dir / f"{STEM}_manifest.tsv", "gray": qa_dir / f"{STEM}_grayscale_QA.png",
    }
    existing = [str(p) for p in paths.values() if p.exists()]
    if existing: raise FileExistsError("Refusing to overwrite: " + "; ".join(existing))
    data, checks, inputs = build_data(source_root)
    write_tsv(paths["source"], data)
    write_tsv(paths["validation"], checks)
    fig = plot(data)
    fig.savefig(paths["pdf"], format="pdf")
    fig.savefig(paths["svg"], format="svg")
    fig.savefig(paths["png"], format="png", dpi=320)
    fig.savefig(paths["tiff"], format="tiff", dpi=400, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    with Image.open(paths["png"]) as image: image.convert("L").save(paths["gray"])
    geometry = [
        {"parameter":"figure_width_mm","value":"297"}, {"parameter":"figure_height_mm","value":"210"}, {"parameter":"png_dpi","value":"320"}, {"parameter":"tiff_dpi","value":"400"},
        {"parameter":"minimum_intended_text_pt","value":"7.6"}, {"parameter":"axis_line_width_pt","value":"1.1"},
        {"parameter":"panel_a_denominator","value":"Bacteria + Archaea reads within each sample"}, {"parameter":"panel_b_denominator","value":"all reads assigned to the respective domain across seven samples"},
        {"parameter":"statistical_test","value":"none; descriptive only"}, {"parameter":"python_version","value":platform.python_version()}, {"parameter":"matplotlib_version","value":mpl.__version__},
        {"parameter":"taxonomy","value":"DADA2 assignTaxonomy; SILVA r144; minBoot=80"},
    ]
    write_tsv(paths["geometry"], geometry)
    summary = [f"PASS: {len(checks)} validation checks", "Canonical taxonomy: taxa_silva144.csv joined by exact sequence to seqtab_final_long.csv", "Archaea reads: 106,314; Bacteria reads: 610,601; total reads: 716,915", "Panel (a) denominator: within-sample Bacteria + Archaea reads", "Panel (b) denominator: all reads assigned to each domain across seven samples", "Statistical test: none; descriptive only"]
    paths["runlog"].write_text("\n".join(summary)+"\n", encoding="utf-8")
    manifest = [{"role":"plotting_script","path":str(Path(__file__).resolve()),"bytes":Path(__file__).stat().st_size,"sha256":sha256(Path(__file__))}]
    for role, p in zip(["input_taxonomy","input_seqtab","input_metadata"], inputs): manifest.append({"role":role,"path":str(p),"bytes":p.stat().st_size,"sha256":sha256(p)})
    for role, p in paths.items():
        if role != "manifest":
            manifest.append({"role":"output_"+role,"path":str(p),"bytes":p.stat().st_size,"sha256":sha256(p)})
    write_tsv(paths["manifest"], manifest)
    print("\n".join(summary))

if __name__ == "__main__": main()
