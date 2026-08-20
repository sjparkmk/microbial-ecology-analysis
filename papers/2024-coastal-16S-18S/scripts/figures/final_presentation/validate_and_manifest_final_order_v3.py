#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

from PIL import Image
from pypdf import PdfReader


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[4]
PKG = ROOT / "revision" / "re_analysis" / "63_figure_presentation_refinement_R_v1"
V2 = ROOT / "revision" / "rev_ms" / "fig_tables_rev" / "final_order_v2"
V3 = ROOT / "revision" / "rev_ms" / "fig_tables_rev" / "final_order_v3"
V2_MANIFEST = V2 / "00_final_order_manifest_v2.tsv"
V3_MANIFEST = V3 / "00_final_order_manifest_v3.tsv"
QA_JSON = PKG / "qc" / "final_order_v3_validation_v1.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def pdf_text(path: Path) -> str:
    reader = PdfReader(str(path))
    if len(reader.pages) != 1:
        raise AssertionError(f"PDF is not one page: {path}")
    return "\n".join(page.extract_text() or "" for page in reader.pages)


def add_pair(mapping: dict[str, Path], final_stem: str, source_stem: Path) -> None:
    for suffix in (".pdf", ".tiff"):
        mapping[final_stem + suffix] = Path(str(source_stem) + suffix)


modified: dict[str, Path] = {}
add_pair(
    modified,
    "01_main_figures/Figure_1_16S_18S_taxonomic_composition_v1",
    PKG / "generated" / "figure1" / "delivery" / "Figure_1_16S_18S_taxonomic_composition_v1",
)
add_pair(
    modified,
    "01_main_figures/Figure_2_16S_18S_alpha_effect_sizes_v1",
    PKG / "generated" / "hierarchy" / "figures" / "Figure_3_16S_18S_alpha_effect_sizes_v1",
)
add_pair(
    modified,
    "01_main_figures/Figure_3_16S_18S_beta_diversity_v1",
    PKG / "generated" / "hierarchy" / "figures" / "Figure_2_16S_18S_beta_diversity_v1",
)
add_pair(
    modified,
    "01_main_figures/Figure_4_paired_sample_type_differential_abundance_v1",
    PKG / "generated" / "daa" / "artifacts_v1" / "figures" / "Figure_4_paired_sample_type_differential_abundance_v1",
)
add_pair(
    modified,
    "01_main_figures/Figure_5_16S_18S_cross_marker_concordance_v1",
    PKG / "generated" / "f5" / "delivery" / "Supplementary_Figure_S13_cross_marker_concordance_R_v1",
)
add_pair(
    modified,
    "02_supplementary_figures/Supplementary_Figure_S2_taxonomic_composition_Other_decomposition_v2",
    PKG / "generated" / "s2" / "delivery" / "Supplementary_Figure_S3_taxonomic_bubble_R_v2",
)
add_pair(
    modified,
    "02_supplementary_figures/Supplementary_Figure_S7_16S_core_ASV_patterns_v1",
    PKG / "generated" / "core_delivery" / "Supplementary_Figure_S5_16S_core_ASV_patterns_R_v5",
)
add_pair(
    modified,
    "02_supplementary_figures/Supplementary_Figure_S8_18S_core_ASV_patterns_v1",
    PKG / "generated" / "core_delivery" / "Supplementary_Figure_S6_18S_core_ASV_patterns_R_v5",
)
add_pair(
    modified,
    "02_supplementary_figures/Supplementary_Figure_S9_core_ASV_phylogeny_abundance_v1",
    PKG / "generated" / "core_delivery" / "Supplementary_Figure_S7_core_ASV_phylogeny_abundance_R_v5",
)


def main() -> None:
    if not V2_MANIFEST.exists() or not V3.exists():
        raise AssertionError("Required final-order package is missing")
    if V3_MANIFEST.exists():
        raise AssertionError(f"Refusing to overwrite existing manifest: {V3_MANIFEST}")

    with V2_MANIFEST.open("r", encoding="utf-8-sig", newline="") as handle:
        prior_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(prior_rows) != 41:
        raise AssertionError(f"Expected 41 prior manifest rows, found {len(prior_rows)}")

    payload = sorted(
        path.relative_to(V3).as_posix()
        for path in V3.rglob("*")
        if path.is_file() and not path.name.startswith("00_")
    )
    if len(payload) != 41:
        raise AssertionError(f"Expected 41 payload files, found {len(payload)}")

    prior_by_rel: dict[str, dict[str, str]] = {}
    marker = "/final_order_v2/"
    for row in prior_rows:
        final_path = row["FinalPath"].replace("\\", "/")
        if marker not in final_path:
            raise AssertionError(f"Unexpected v2 FinalPath: {final_path}")
        prior_by_rel[final_path.split(marker, 1)[1]] = row
    if set(payload) != set(prior_by_rel):
        raise AssertionError("v3 payload paths do not match the v2 order/package contract")
    if set(modified) - set(payload):
        raise AssertionError("A modified-output mapping is absent from the v3 payload")

    manifest_rows: list[dict[str, object]] = []
    modified_matches = 0
    unchanged_matches = 0
    for row in prior_rows:
        final_path_v2 = row["FinalPath"].replace("\\", "/")
        rel = final_path_v2.split(marker, 1)[1]
        prior = V2 / rel
        final = V3 / rel
        generation = modified.get(rel, prior)
        for path in (prior, final, generation):
            if not path.exists():
                raise AssertionError(f"Missing manifest target: {path}")

        prior_hash = sha256(prior)
        final_hash = sha256(final)
        generation_hash = sha256(generation)
        if rel in modified:
            if final_hash != generation_hash or final.stat().st_size != generation.stat().st_size:
                raise AssertionError(f"Modified delivery does not match its generated source: {rel}")
            modified_matches += 1
            disposition = "PRESENTATION_REFINED_R_V1"
            validation = "SOURCE_HASH_AND_FORMAT_PASS"
        else:
            if final_hash != prior_hash or final.stat().st_size != prior.stat().st_size:
                raise AssertionError(f"Unrequested payload changed from v2: {rel}")
            unchanged_matches += 1
            disposition = "BYTE_IDENTICAL_FROM_V2"
            validation = "BYTE_IDENTICAL_PASS"

        manifest_rows.append(
            {
                "ArtifactType": row["ArtifactType"],
                "FinalOrder": row["FinalOrder"],
                "Format": row["Format"],
                "ExpectedPanelLabels": row["ExpectedPanelLabels"],
                "PriorPackagePath": relative(prior),
                "PriorBytes": prior.stat().st_size,
                "PriorSHA256": prior_hash,
                "GenerationSourcePath": relative(generation),
                "GenerationSourceBytes": generation.stat().st_size,
                "GenerationSourceSHA256": generation_hash,
                "FinalPath": relative(final),
                "FinalBytes": final.stat().st_size,
                "FinalSHA256": final_hash,
                "Disposition": disposition,
                "Validation": validation,
            }
        )

    pdfs = sorted((V3 / "01_main_figures").glob("*.pdf")) + sorted(
        (V3 / "02_supplementary_figures").glob("*.pdf")
    )
    tiffs = sorted((V3 / "01_main_figures").glob("*.tiff")) + sorted(
        (V3 / "02_supplementary_figures").glob("*.tiff")
    )
    if len(pdfs) != 17 or len(tiffs) != 17:
        raise AssertionError("Expected 17 PDF/TIFF figure pairs")

    pdf_texts = {path.name: pdf_text(path) for path in pdfs}
    tiff_records = []
    for path in tiffs:
        with Image.open(path) as image:
            dpi = image.info.get("dpi", (0, 0))
            compression = str(image.info.get("compression", ""))
            record = {
                "file": path.name,
                "width_px": image.width,
                "height_px": image.height,
                "mode": image.mode,
                "dpi_x": float(dpi[0]),
                "dpi_y": float(dpi[1]),
                "compression": compression,
            }
            if image.mode != "RGB":
                raise AssertionError(f"TIFF is not RGB: {path}")
            if abs(float(dpi[0]) - 600) > 1 or abs(float(dpi[1]) - 600) > 1:
                raise AssertionError(f"TIFF is not 600 dpi: {path}")
            if "lzw" not in compression.lower():
                raise AssertionError(f"TIFF is not LZW-compressed: {path}")
            tiff_records.append(record)

    f2 = pdf_texts["Figure_2_16S_18S_alpha_effect_sizes_v1.pdf"]
    q_labels = re.findall(r"q\s*=\s*0\.\d{3}", f2)
    if len(q_labels) != 12 or re.search(r"q\s*=\s*\.\d", f2):
        raise AssertionError("Main Figure 2 q-label formatting check failed")

    f3 = pdf_texts["Figure_3_16S_18S_beta_diversity_v1.pdf"]
    if "A 16S" not in f3 or "B 18S" not in f3:
        raise AssertionError("Main Figure 3 panel labels are missing")

    s2 = pdf_texts["Supplementary_Figure_S2_taxonomic_composition_Other_decomposition_v2.pdf"]
    if s2.count("Ulleungdo") < 2:
        raise AssertionError("Supplementary Figure S2 does not retain both Ulleungdo labels")
    s2_reader = PdfReader(str(V3 / "02_supplementary_figures" /
                              "Supplementary_Figure_S2_taxonomic_composition_Other_decomposition_v2.pdf"))
    s2_width_mm = float(s2_reader.pages[0].mediabox.width) * 25.4 / 72.0
    if abs(s2_width_mm - 210.0) > 0.5:
        raise AssertionError(f"Supplementary Figure S2 width drifted: {s2_width_mm}")

    for name in (
        "Supplementary_Figure_S7_16S_core_ASV_patterns_v1.pdf",
        "Supplementary_Figure_S8_18S_core_ASV_patterns_v1.pdf",
    ):
        text = " ".join(pdf_texts[name].split())
        required = ("Mean relative abundance (%)", "Prevalence (count >0, %)", "25", "50", "75", "100")
        if not all(token in text for token in required):
            raise AssertionError(f"Prevalence/fill legend check failed: {name}")

    s9 = " ".join(pdf_texts["Supplementary_Figure_S9_core_ASV_phylogeny_abundance_v1.pdf"].split())
    if "Core status" not in s9 or "prevalence >=80%; count >0" not in s9:
        raise AssertionError("Supplementary Figure S9 core-status legend check failed")

    figure_rel = [p for p in payload if p.startswith(("01_main_figures/", "02_supplementary_figures/"))]
    forbidden_names = [p for p in figure_rel if re.search(r"random[_ -]?forest|ncm|tnst", p, re.I)]
    if forbidden_names:
        raise AssertionError(f"Excluded figure type present: {forbidden_names}")

    fieldnames = list(manifest_rows[0])
    with V3_MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)

    qa = {
        "status": "PASS",
        "blockers": 0,
        "payload_files": len(payload),
        "main_figure_files": 10,
        "supplementary_figure_files": 24,
        "supplementary_table_files": 7,
        "modified_source_matches": modified_matches,
        "unchanged_v2_byte_matches": unchanged_matches,
        "pdf_one_page_pass": len(pdfs),
        "tiff_rgb_600dpi_lzw_pass": len(tiffs),
        "figure2_q_labels_0xxx": len(q_labels),
        "figure3_shared_coordinate_window": {
            "x": [-1.240169, 1.178654],
            "y": [-1.717716, 1.489935],
            "geometry": "coord_equal",
        },
        "supplementary_figure_s2_width_mm": s2_width_mm,
        "supplementary_figure_s2_ulleungdo_text_count": s2.count("Ulleungdo"),
        "supplementary_figure_s7_s8_prevalence_breaks_percent": [25, 50, 75, 100],
        "supplementary_figure_s9_marker_semantics": "binary core status; prevalence >=80% at count >0",
        "excluded_rf_ncm_tnst_figure_hits": forbidden_names,
        "tiff_records": tiff_records,
        "manifest": relative(V3_MANIFEST),
    }
    QA_JSON.parent.mkdir(parents=True, exist_ok=True)
    with QA_JSON.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(qa, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    print("FINAL_ORDER_V3_VALIDATION_PASS")


if __name__ == "__main__":
    main()
