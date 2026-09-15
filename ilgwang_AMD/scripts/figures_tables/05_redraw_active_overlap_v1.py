from __future__ import annotations

import argparse
import csv
import hashlib
import math
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


A4_PORTRAIT = (8.2677165354, 11.6929133858)
INK = "#111111"
SECONDARY = "#5E6462"
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
CENTERS = {
    "P_biofilm": (-0.72, 0.35),
    "W_biofilm": (0.72, 0.35),
    "water": (0.0, -0.70),
}
RADIUS = 1.45


NUMBER_POSITIONS = {
    "P_biofilm_only": (-1.35, 0.52),
    "W_biofilm_only": (1.35, 0.52),
    "water_only": (0.0, -1.62),
    "P_biofilm_and_W_biofilm_only": (0.0, 1.22),
    "P_biofilm_and_water_only": (-0.85, -0.38),
    "W_biofilm_and_water_only": (0.85, -0.38),
    "all_three": (0.0, 0.10),
}
EXPECTED_MEMBERSHIP = {
    "P_biofilm_only": {"P_biofilm"},
    "W_biofilm_only": {"W_biofilm"},
    "water_only": {"water"},
    "P_biofilm_and_W_biofilm_only": {"P_biofilm", "W_biofilm"},
    "P_biofilm_and_water_only": {"P_biofilm", "water"},
    "W_biofilm_and_water_only": {"W_biofilm", "water"},
    "all_three": {"P_biofilm", "W_biofilm", "water"},
}
EXPECTED_COUNTS = {
    "P_biofilm_only": 102,
    "W_biofilm_only": 127,
    "water_only": 133,
    "P_biofilm_and_W_biofilm_only": 84,
    "P_biofilm_and_water_only": 24,
    "W_biofilm_and_water_only": 31,
    "all_three": 162,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_values(path: Path) -> dict[str, int]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        values = {row["region"]: int(row["ASV_count"]) for row in csv.DictReader(handle, delimiter="\t")}
    if values != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected overlap data: {values}")
    if sum(values.values()) != 663:
        raise ValueError("Overlap regions do not sum to 663 ASVs")
    return values


def membership_margin(region: str, x: float, y: float) -> tuple[float, str]:
    expected_inside = EXPECTED_MEMBERSHIP[region]
    margins: list[tuple[float, str]] = []
    for type_code, (cx, cy) in CENTERS.items():
        distance = math.hypot(x - cx, y - cy)
        if type_code in expected_inside:
            margins.append((RADIUS - distance, f"inside:{type_code}"))
        else:
            margins.append((distance - RADIUS, f"outside:{type_code}"))
    return min(margins, key=lambda item: item[0])


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = [output_dir / f"{stem}.{suffix}" for suffix in ("pdf", "svg", "png", "tiff")]
    fig.savefig(outputs[0], facecolor="white")
    fig.savefig(outputs[1], facecolor="white")
    fig.savefig(outputs[2], dpi=320, facecolor="white")
    fig.savefig(outputs[3], dpi=400, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return outputs


def draw(values: dict[str, int], output_dir: Path) -> list[Path]:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "axes.linewidth": 1.1,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=A4_PORTRAIT)
    ax = fig.add_axes([0.055, 0.15, 0.89, 0.70])
    ax.set_aspect("equal")
    ax.axis("off")

    for code in ("P_biofilm", "W_biofilm", "water"):
        ax.add_patch(
            Circle(
                CENTERS[code],
                RADIUS,
                facecolor=TYPE_COLORS[code],
                edgecolor=TYPE_COLORS[code],
                linewidth=1.8,
                alpha=0.43,
            )
        )

    for region, position in NUMBER_POSITIONS.items():
        ax.text(
            *position,
            str(values[region]),
            ha="center",
            va="center",
            fontsize=12.0,
            fontweight="medium",
            color=INK,
            zorder=5,
        )

    label_fontsize = 12.5
    ax.text(
        -1.08,
        1.94,
        TYPE_LABELS["P_biofilm"],
        ha="center",
        va="bottom",
        color=TYPE_COLORS["P_biofilm"],
        fontsize=label_fontsize,
        fontweight="bold",
        linespacing=1.05,
    )
    ax.text(
        1.08,
        1.94,
        TYPE_LABELS["W_biofilm"],
        ha="center",
        va="bottom",
        color=TYPE_COLORS["W_biofilm"],
        fontsize=label_fontsize,
        fontweight="bold",
        linespacing=1.05,
    )
    ax.text(
        0.0,
        -2.25,
        TYPE_LABELS["water"],
        ha="center",
        va="top",
        color=TYPE_COLORS["water"],
        fontsize=label_fontsize,
        fontweight="bold",
    )
    ax.set_xlim(-2.25, 2.25)
    ax.set_ylim(-2.55, 2.30)

    fig.suptitle(
        "ASV detection overlaps across the three microhabitat Types",
        x=0.08,
        y=0.93,
        ha="left",
        fontsize=12.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.08,
        0.075,
        "Presence was defined as ≥1 read in ≥1 sample of a Type (n = 663 ASVs). Circle areas are schematic and are not scaled to counts.",
        ha="left",
        va="center",
        fontsize=8.2,
        color=SECONDARY,
    )
    return save_figure(fig, output_dir, "SFig7_ASV_Overlap_A4_v5")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=Path("/input/youre/data/"))
    parser.add_argument("--output-root", type=Path, default=Path("/output/your/data"))
    args = parser.parse_args()
    final_dir = args.output_root.resolve()
    input_path = args.input_root.resolve() / "source_data" / "SFig7_ASV_overlap_regions_v4.tsv"
    output_dir = final_dir / "build" / "overlap_stage" / "figures"
    log_dir = final_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    values = read_values(input_path)
    validation_rows = []
    for region, (x, y) in NUMBER_POSITIONS.items():
        margin, limiting_boundary = membership_margin(region, x, y)
        if margin <= 0:
            raise ValueError(f"Number position is outside its intended region: {region}")
        validation_rows.append({
            "region": region,
            "ASV_count": values[region],
            "x": f"{x:.3f}",
            "y": f"{y:.3f}",
            "minimum_center_to_boundary_clearance": f"{margin:.6f}",
            "limiting_boundary": limiting_boundary,
            "status": "PASS",
        })

    outputs = draw(values, output_dir)
    validation_path = log_dir / "86_sfig7_layout_QA_v1.tsv"
    with validation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(validation_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(validation_rows)

    run_log = log_dir / "85_sfig7_layout_update_v1.log"
    lines = [
        "Supplementary Figure S7 layout update",
        f"UTC timestamp: {datetime.now(timezone.utc).isoformat()}",
        f"Python: {platform.python_version()}",
        f"matplotlib: {mpl.__version__}",
        "",
        f"INPUT\t{input_path}\tSHA256={sha256(input_path)}",
        "CHANGE\tMoved the 24 and 31 labels away from the opposing circle boundaries and recentered all region counts.",
        "CHANGE\tIncreased all three sample-Type labels to 12.5 pt bold.",
        f"VALIDATION\tAll seven number anchors lie in their intended mutually exclusive regions; minimum anchor clearance={min(float(row['minimum_center_to_boundary_clearance']) for row in validation_rows):.6f} data units.",
        "VALIDATION\tOverlap counts unchanged and sum to 663.",
        "",
        "OUTPUTS",
    ]
    lines.extend(f"{path}\tSHA256={sha256(path)}" for path in outputs)
    lines.append(f"{validation_path}\tSHA256={sha256(validation_path)}")
    run_log.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"WROTE {run_log}")
    for path in outputs:
        print(f"WROTE {path}")
    print(f"WROTE {validation_path}")


if __name__ == "__main__":
    main()
