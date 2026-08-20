#!/usr/bin/env python3
"""Create the versioned PR2 lineage-flag correction supplement.

The v1 parallel taxonomy values are not changed.  This supplement fixes only
the PR2 Metazoa/Chloroplastida flags so that the audited PR2 subdivision
crosswalk is considered in addition to literal SINTAX rank values.  It also
adds selected QC flags to a v2 copy of the Additional-Data taxonomy table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Iterable, Sequence


VERSION = "v2"
CANONICAL_PR2_COLUMNS = (
    "pr2_canonical_domain",
    "pr2_canonical_kingdom",
    "pr2_canonical_division_or_phylum",
    "pr2_canonical_class",
    "pr2_canonical_order",
    "pr2_canonical_family",
    "pr2_canonical_genus",
    "pr2_canonical_species",
    "pr2_canonical_subdivision",
)
APPENDED_LEDGER_COLUMNS = (
    "analysis_keep",
    "exclusion_reason",
    "flag_pr2_metazoa",
    "flag_pr2_chloroplastida",
    "flag_pr2_nucleomorph",
    "flag_pr2_non_eukaryota",
    "flag_silva_metazoa",
    "flag_silva_chloroplastida",
    "flag_crossdb_discordant",
    "flag_zero_total_count",
)
OUTPUTS = {
    "parallel_v2": "results/18s_taxonomy_parallel_pr2_silva123_all_v2.tsv",
    "ledger_v2": "qc/18s_asv_flags_ledger_v2.tsv",
    "summary_v2": "provenance/package_summary_18s_parallel_taxonomy_v2.json",
    "readme_v2": "provenance/README_18s_parallel_taxonomy_v2.md",
    "manifest_v2": "provenance/sha256_manifest_18s_parallel_taxonomy_v2.tsv",
}


class ValidationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    require(path.is_file(), f"Missing TSV: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, f"Missing header: {path}")
        rows = list(reader)
    return list(reader.fieldnames), rows


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[3]


def package_paths(project_root: Path, output_root: Path) -> dict[str, Path]:
    return {
        "parallel_v1": output_root / "results" / "18s_taxonomy_parallel_pr2_silva123_all_v1.tsv",
        "ledger_v1": output_root / "qc" / "18s_asv_flags_ledger_v1.tsv",
        "summary_v1": output_root / "provenance" / "package_summary_18s_parallel_taxonomy_v1.json",
        "helper_v1": project_root / "revision" / "re_analysis" / "scripts" / "build_18s_parallel_taxonomy_v1.py",
    }


def corrected_rows(
    paths: dict[str, Path],
) -> tuple[list[str], list[dict[str, str]], list[str], list[dict[str, str]], dict[str, object]]:
    parallel_header, parallel_rows = read_tsv(paths["parallel_v1"])
    ledger_header, ledger_rows = read_tsv(paths["ledger_v1"])
    require(parallel_header and parallel_header[0] == "ASV_ID", "Unexpected v1 parallel header")
    require(ledger_header and ledger_header[0] == "ASV_ID", "Unexpected v1 ledger header")
    for column in CANONICAL_PR2_COLUMNS:
        require(column in parallel_header, f"Missing PR2 canonical column: {column}")
    for column in APPENDED_LEDGER_COLUMNS:
        require(column in ledger_header, f"Missing v1 ledger column: {column}")

    parallel_ids = [row["ASV_ID"] for row in parallel_rows]
    ledger_ids = [row["ASV_ID"] for row in ledger_rows]
    require(len(parallel_ids) == len(set(parallel_ids)), "Duplicate ASV IDs in v1 parallel table")
    require(len(ledger_ids) == len(set(ledger_ids)), "Duplicate ASV IDs in v1 ledger")
    require(set(parallel_ids) == set(ledger_ids), "v1 parallel/ledger ASV sets differ")
    ledger_by_id = {row["ASV_ID"]: dict(row) for row in ledger_rows}

    corrected_parallel: list[dict[str, str]] = []
    corrected_ledger: list[dict[str, str]] = []
    metazoa_count = 0
    chloroplastida_count = 0
    for parallel in parallel_rows:
        feature_id = parallel["ASV_ID"]
        ledger = ledger_by_id[feature_id]
        canonical_values = {
            parallel[column].strip().casefold()
            for column in CANONICAL_PR2_COLUMNS
            if parallel[column].strip() and parallel[column].strip() != "Unassigned"
        }
        is_metazoa = "metazoa" in canonical_values
        is_chloroplastida = "chloroplastida" in canonical_values
        metazoa_count += is_metazoa
        chloroplastida_count += is_chloroplastida
        ledger["flag_pr2_metazoa"] = bool_text(is_metazoa)
        ledger["flag_pr2_chloroplastida"] = bool_text(is_chloroplastida)
        corrected_ledger.append(ledger)

        enriched = dict(parallel)
        for column in APPENDED_LEDGER_COLUMNS:
            enriched[column] = ledger[column]
        corrected_parallel.append(enriched)

    with paths["summary_v1"].open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    require(summary["feature_contract"]["all_features"] == len(parallel_rows), "v1 summary feature count differs")
    summary["package_version"] = VERSION
    summary["flag_counts"]["flag_pr2_metazoa"] = metazoa_count
    summary["flag_counts"]["flag_pr2_chloroplastida"] = chloroplastida_count
    summary["v2_correction"] = {
        "scope": "lineage flags and Additional-Data flag columns only; counts and taxonomy assignments unchanged",
        "reason": "v1 did not include the PR2 division-subdivision crosswalk when setting Metazoa/Chloroplastida flags",
        "pr2_metazoa_asvs_after_correction": metazoa_count,
        "pr2_chloroplastida_asvs_after_correction": chloroplastida_count,
        "v1_parallel_taxonomy_sha256": sha256_file(paths["parallel_v1"]),
        "v1_ledger_sha256": sha256_file(paths["ledger_v1"]),
    }
    return (
        [*parallel_header, *APPENDED_LEDGER_COLUMNS],
        corrected_parallel,
        ledger_header,
        corrected_ledger,
        summary,
    )


def readme_text(summary: dict[str, object]) -> str:
    correction = summary["v2_correction"]
    return f"""# 18S PR2/SILVA parallel taxonomy flag supplement ({VERSION})

Use `results/18s_taxonomy_parallel_pr2_silva123_all_v2.tsv` as the
Additional-Data parallel taxonomy table. It preserves every v1 PR2/SILVA
taxonomy cell and appends selected analysis/QC flags.

The v2 correction changes no sequence, count, taxonomy assignment, confidence,
agreement call, or analysis-retention decision. It corrects PR2 lineage flags
by considering the audited `division-subdivision` crosswalk. Corrected counts:
PR2 Metazoa = {correction['pr2_metazoa_asvs_after_correction']:,} ASVs;
PR2 Chloroplastida = {correction['pr2_chloroplastida_asvs_after_correction']:,}
ASVs. Existing v1 files remain unchanged for provenance.
"""


def expected_content(project_root: Path, output_root: Path):
    paths = package_paths(project_root, output_root)
    for path in paths.values():
        require(path.is_file(), f"Missing source: {path}")
    return paths, corrected_rows(paths)


def write_package(stage_root: Path, project_root: Path, output_root: Path) -> None:
    paths, content = expected_content(project_root, output_root)
    parallel_header, parallel_rows, ledger_header, ledger_rows, summary = content
    write_tsv(stage_root / OUTPUTS["parallel_v2"], parallel_header, parallel_rows)
    write_tsv(stage_root / OUTPUTS["ledger_v2"], ledger_header, ledger_rows)
    summary_path = stage_root / OUTPUTS["summary_v2"]
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    readme_path = stage_root / OUTPUTS["readme_v2"]
    readme_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path.write_text(readme_text(summary), encoding="utf-8", newline="\n")

    rows: list[tuple[object, ...]] = []
    for label, path in sorted(paths.items()):
        rows.append(("source", label, str(path.resolve()), path.stat().st_size, sha256_file(path)))
    for label in ("parallel_v2", "ledger_v2", "summary_v2", "readme_v2"):
        path = stage_root / OUTPUTS[label]
        rows.append(("output", label, OUTPUTS[label], path.stat().st_size, sha256_file(path)))
    manifest = stage_root / OUTPUTS["manifest_v2"]
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("category", "label", "path", "bytes", "sha256"))
        writer.writerows(rows)


def check_package(
    project_root: Path,
    output_root: Path,
    *,
    source_root: Path | None = None,
) -> dict[str, object]:
    paths, content = expected_content(project_root, source_root or output_root)
    parallel_header, parallel_rows, ledger_header, ledger_rows, summary = content
    observed_parallel_header, observed_parallel_rows = read_tsv(output_root / OUTPUTS["parallel_v2"])
    observed_ledger_header, observed_ledger_rows = read_tsv(output_root / OUTPUTS["ledger_v2"])
    require(observed_parallel_header == parallel_header, "v2 parallel header mismatch")
    require(observed_parallel_rows == parallel_rows, "v2 parallel rows differ from recomputed correction")
    require(observed_ledger_header == ledger_header, "v2 ledger header mismatch")
    require(observed_ledger_rows == ledger_rows, "v2 ledger rows differ from recomputed correction")
    observed_summary = json.loads((output_root / OUTPUTS["summary_v2"]).read_text(encoding="utf-8"))
    require(observed_summary == summary, "v2 summary differs from recomputed correction")
    require((output_root / OUTPUTS["readme_v2"]).read_text(encoding="utf-8") == readme_text(summary), "v2 README mismatch")

    manifest_header, manifest_rows = read_tsv(output_root / OUTPUTS["manifest_v2"])
    require(manifest_header == ["category", "label", "path", "bytes", "sha256"], "v2 manifest header mismatch")
    require(len(manifest_rows) == len(paths) + 4, "v2 manifest row count mismatch")
    for row in manifest_rows:
        path = Path(row["path"]) if row["category"] == "source" else output_root / row["path"]
        require(path.is_file(), f"Manifested v2 path missing: {path}")
        require(path.stat().st_size == int(row["bytes"]), f"Manifested v2 size drift: {path}")
        require(sha256_file(path) == row["sha256"], f"Manifested v2 hash drift: {path}")
    return {
        "status": "CHECK_OK",
        "features": len(parallel_rows),
        "pr2_metazoa_asvs": summary["flag_counts"]["flag_pr2_metazoa"],
        "pr2_chloroplastida_asvs": summary["flag_counts"]["flag_pr2_chloroplastida"],
        "taxonomy_or_counts_changed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build/check the versioned 18S PR2 lineage-flag correction supplement.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true", help="Create v2 supplement; refuses existing targets")
    mode.add_argument("--check", action="store_true", help="Recompute and validate existing v2 supplement")
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    output_root = (args.output_dir or project_root / "revision" / "re_analysis" / "02_18s_unoise_pr2_v1").resolve()
    require(output_root.is_dir(), f"Missing completed 18S package: {output_root}")
    if args.check:
        print(json.dumps(check_package(project_root, output_root), ensure_ascii=False, indent=2))
        return 0
    targets = [output_root / relative for relative in OUTPUTS.values()]
    require(not any(path.exists() for path in targets), "Refusing to overwrite an existing v2 supplement target")
    stage_root = Path(tempfile.mkdtemp(prefix=".taxonomy_flag_v2_", dir=output_root))
    try:
        write_package(stage_root, project_root, output_root)
        check_package(project_root, stage_root, source_root=output_root)
        for relative in OUTPUTS.values():
            source = stage_root / relative
            target = output_root / relative
            require(source.is_file() and not target.exists(), f"Unsafe v2 install target: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            source.replace(target)
    finally:
        if stage_root.exists() and stage_root.parent == output_root and stage_root.name.startswith(".taxonomy_flag_v2_"):
            shutil.rmtree(stage_root)
    print(json.dumps(check_package(project_root, output_root), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(2)
