#!/usr/bin/env python3
"""Build and verify the 75-sample 16S rarefaction package at 20,000 reads.

The non-rarefied 16S analysis package is immutable input. Creation requires
``--run`` and uses a sibling staging directory followed by an atomic rename.
``--check`` re-reads every artifact, verifies source/output manifests and
linked feature order, and asks the isolated R worker to recompute the seed-1
rarefaction twice. This module intentionally performs no diversity analysis,
statistical testing, taxonomic aggregation, table generation, or plotting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


VERSION = "1.0.0"
PACKAGE_ID = "04_input_16s_rarefied_20000_v1"
SOURCE_PACKAGE_ID = "01_input_16s_75_v1"
DEPTH = 20_000
SEED = 1
EXPECTED_FEATURES = 17_562
EXPECTED_SAMPLES = 75
EXPECTED_SOURCE_TOTAL = 3_475_142
EXPECTED_SOURCE_MIN = 20_837
EXPECTED_SOURCE_MAX = 89_683
EXPECTED_TOTAL = DEPTH * EXPECTED_SAMPLES
RNG_CONTRACT = {
    "kind": "Mersenne-Twister",
    "normal_kind": "Inversion",
    "sample_kind": "Rejection",
}

SOURCE_FILES = {
    "counts": "16s_asv_counts_analysis_75_v1.tsv",
    "taxonomy": "16s_taxonomy_analysis_75_v1.tsv",
    "fasta": "16s_rep_seqs_analysis_75_v1.fasta",
    "sintax": "16s_sintax_analysis_75_v1.tsv",
    "metadata": "metadata_16s_75_v1.tsv",
    "summary": "package_summary_16s_75_v1.json",
    "manifest": "sha256_manifest_16s_75_v1.tsv",
}

OUTPUT_FILES = {
    "counts_all": "results/16s_asv_counts_rarefied_20000_all_v1.tsv",
    "counts_analysis": "results/16s_asv_counts_rarefied_20000_analysis_v1.tsv",
    "taxonomy": "results/16s_taxonomy_rarefied_20000_analysis_v1.tsv",
    "fasta": "results/16s_rep_seqs_rarefied_20000_analysis_v1.fasta",
    "sintax": "results/16s_sintax_rarefied_20000_analysis_v1.tsv",
    "sample_ledger": "qc/16s_rarefaction_sample_ledger_20000_v1.tsv",
    "feature_ledger": "qc/16s_rarefaction_feature_ledger_20000_v1.tsv",
    "validation": "qc/16s_rarefaction_validation_20000_v1.json",
    "metadata": "provenance/metadata_16s_75_v1.tsv",
    "config": "provenance/run_config_16s_rarefied_20000_v1.json",
    "summary": "provenance/package_summary_16s_rarefied_20000_v1.json",
    "readme": "provenance/README_16s_rarefied_20000_v1.md",
    "manifest": "provenance/sha256_manifest_16s_rarefied_20000_v1.tsv",
    "r_stdout": "logs/rarefy_16s_20000_v1.stdout.log",
    "r_stderr": "logs/rarefy_16s_20000_v1.stderr.log",
}

SAMPLE_LEDGER_HEADER = (
    "SampleID",
    "source_depth",
    "rarefied_depth",
    "removed_reads",
    "retention_fraction",
    "rarefaction_depth",
    "seed",
    "site",
    "Group",
    "type",
    "subgroup",
    "subID",
)
FEATURE_LEDGER_HEADER = (
    "ASV_ID",
    "source_total",
    "rarefied_total",
    "removed_reads",
    "retention_fraction",
    "source_positive_sample_count",
    "rarefied_positive_sample_count",
    "retained_in_analysis",
    "post_rarefaction_status",
)
MANIFEST_HEADER = ("scope", "label", "path", "bytes", "sha256")

SOURCE_MANIFEST_CONTRACT = {
    "counts": ("output", "counts_analysis"),
    "taxonomy": ("output", "taxonomy_analysis"),
    "fasta": ("output", "fasta_analysis"),
    "sintax": ("output", "sintax_analysis"),
    "metadata": ("output", "metadata"),
    "summary": ("output", "summary"),
}


class ValidationError(RuntimeError):
    """Raised when a source or output violates the frozen contract."""


@dataclass(frozen=True)
class CountTable:
    id_header: str
    samples: tuple[str, ...]
    feature_ids: tuple[str, ...]
    rows: tuple[tuple[int, ...], ...]

    @property
    def sample_totals(self) -> tuple[int, ...]:
        totals = [0] * len(self.samples)
        for row in self.rows:
            for index, value in enumerate(row):
                totals[index] += value
        return tuple(totals)

    @property
    def feature_totals(self) -> tuple[int, ...]:
        return tuple(sum(row) for row in self.rows)


@dataclass(frozen=True)
class KeyedTable:
    header: tuple[str, ...]
    order: tuple[str, ...]
    rows: dict[str, tuple[str, ...]]


@dataclass(frozen=True)
class HeaderlessKeyedTable:
    width: int
    order: tuple[str, ...]
    rows: dict[str, tuple[str, ...]]


@dataclass(frozen=True)
class FastaData:
    order: tuple[str, ...]
    sequences: dict[str, str]


@dataclass(frozen=True)
class SourceData:
    paths: dict[str, Path]
    hashes: dict[str, str]
    sizes: dict[str, int]
    manifest_validation: tuple[dict[str, object], ...]
    counts: CountTable
    taxonomy: KeyedTable
    sintax: HeaderlessKeyedTable
    fasta: FastaData
    metadata: KeyedTable


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_path(path: Path) -> str:
    return os.path.normcase(str(path.resolve()))


def ensure_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    require(not duplicates, f"Duplicate {label}: {sorted(duplicates)[:10]}")


def read_count_table(path: Path, expected_features: int | None = EXPECTED_FEATURES) -> CountTable:
    require(path.is_file(), f"Missing count table: {path}")
    feature_ids: list[str] = []
    rows: list[tuple[int, ...]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValidationError(f"Empty count table: {path}") from exc
        require(len(header) == EXPECTED_SAMPLES + 1, f"Unexpected count-table width: {path}")
        require(header[0] == "#OTUID", f"Expected #OTUID as first header: {path}")
        samples = tuple(header[1:])
        require(all(samples), f"Empty sample ID: {path}")
        ensure_unique(samples, f"sample IDs in {path}")
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            require(len(row) == len(header), f"{path}:{line_number} column-count mismatch")
            feature_id = row[0]
            require(re.fullmatch(r"ASV_[0-9]+", feature_id) is not None, f"{path}:{line_number} invalid ASV ID")
            values: list[int] = []
            for column_number, cell in enumerate(row[1:], start=2):
                require(re.fullmatch(r"\d+", cell) is not None, f"{path}:{line_number}:{column_number} invalid count {cell!r}")
                values.append(int(cell))
            feature_ids.append(feature_id)
            rows.append(tuple(values))
    ensure_unique(feature_ids, f"feature IDs in {path}")
    if expected_features is not None:
        require(len(feature_ids) == expected_features, f"Expected {expected_features} ASVs in {path}; observed {len(feature_ids)}")
    require(len(samples) == EXPECTED_SAMPLES, f"Expected {EXPECTED_SAMPLES} samples in {path}; observed {len(samples)}")
    require(not any(re.match(r"^SS_BS(?:_|$)", sample) for sample in samples), f"SS_BS present in {path}")
    return CountTable(header[0], samples, tuple(feature_ids), tuple(rows))


def read_keyed_tsv(path: Path, id_header: str) -> KeyedTable:
    require(path.is_file(), f"Missing TSV: {path}")
    order: list[str] = []
    rows: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise ValidationError(f"Empty TSV: {path}") from exc
        require(header and header[0] == id_header, f"Expected first column {id_header!r} in {path}")
        for line_number, row_list in enumerate(reader, start=2):
            if not row_list or all(cell == "" for cell in row_list):
                continue
            require(len(row_list) == len(header), f"{path}:{line_number} column-count mismatch")
            row = tuple(row_list)
            key = row[0]
            require(bool(key), f"{path}:{line_number} empty ID")
            require(key not in rows, f"Duplicate ID {key} in {path}")
            order.append(key)
            rows[key] = row
    return KeyedTable(header, tuple(order), rows)


def read_headerless_keyed_tsv(path: Path, expected_width: int = 4) -> HeaderlessKeyedTable:
    require(path.is_file(), f"Missing headerless TSV: {path}")
    order: list[str] = []
    rows: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row_list in enumerate(reader, start=1):
            if not row_list or all(cell == "" for cell in row_list):
                continue
            require(len(row_list) == expected_width, f"{path}:{line_number} expected {expected_width} columns")
            row = tuple(row_list)
            key = row[0]
            require(re.fullmatch(r"ASV_[0-9]+", key) is not None, f"{path}:{line_number} invalid ASV ID")
            require(key not in rows, f"Duplicate ID {key} in {path}")
            order.append(key)
            rows[key] = row
    require(bool(order), f"Empty headerless TSV: {path}")
    return HeaderlessKeyedTable(expected_width, tuple(order), rows)


def read_fasta(path: Path) -> FastaData:
    require(path.is_file(), f"Missing FASTA: {path}")
    order: list[str] = []
    sequences: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal current_id, chunks
        if current_id is None:
            return
        sequence = "".join(chunks).upper()
        require(bool(sequence), f"Empty FASTA sequence for {current_id}: {path}")
        require(current_id not in sequences, f"Duplicate FASTA ID {current_id}: {path}")
        order.append(current_id)
        sequences[current_id] = sequence
        current_id = None
        chunks = []

    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                require(re.fullmatch(r"ASV_[0-9]+", current_id) is not None, f"Invalid FASTA ID at {path}:{line_number}")
            else:
                require(current_id is not None, f"Sequence before FASTA header at {path}:{line_number}")
                require(re.fullmatch(r"[ACGTUN.-]+", line.upper()) is not None, f"Invalid sequence at {path}:{line_number}")
                chunks.append(line)
    flush()
    return FastaData(tuple(order), sequences)


def read_source_manifest(path: Path) -> dict[tuple[str, str], tuple[str, int, str]]:
    require(path.is_file(), f"Missing source manifest: {path}")
    records: dict[tuple[str, str], tuple[str, int, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["category", "label", "path", "bytes", "sha256"], f"Unexpected source-manifest header: {path}")
        for line_number, row in enumerate(reader, start=2):
            key = (row["category"], row["label"])
            require(key not in records, f"Duplicate source-manifest record {key} at {path}:{line_number}")
            require(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None, f"Invalid SHA-256 at {path}:{line_number}")
            records[key] = (row["path"], int(row["bytes"]), row["sha256"])
    return records


def resolve_source_manifest_record(manifest_path: Path, recorded_path: str) -> Path:
    candidate = Path(recorded_path)
    if candidate.is_absolute():
        return candidate.resolve()
    return (manifest_path.parent / candidate).resolve()


def validate_source_manifest_contract(paths: Mapping[str, Path]) -> tuple[dict[str, str], dict[str, int], tuple[dict[str, object], ...]]:
    manifest_path = paths["manifest"]
    records = read_source_manifest(manifest_path)
    hashes: dict[str, str] = {}
    sizes: dict[str, int] = {}
    validations: list[dict[str, object]] = []
    for source_label, (category, record_label) in SOURCE_MANIFEST_CONTRACT.items():
        key = (category, record_label)
        require(key in records, f"Missing {key} from {manifest_path}")
        recorded_path, recorded_bytes, recorded_hash = records[key]
        expected_path = paths[source_label].resolve()
        resolved_recorded_path = resolve_source_manifest_record(manifest_path, recorded_path)
        require(normalized_path(resolved_recorded_path) == normalized_path(expected_path), f"Source-manifest path mismatch for {source_label}")
        require(expected_path.is_file(), f"Missing protected source: {expected_path}")
        actual_bytes = expected_path.stat().st_size
        actual_hash = sha256_file(expected_path)
        require(actual_bytes == recorded_bytes, f"Source-manifest byte mismatch for {source_label}")
        require(actual_hash == recorded_hash, f"Source-manifest SHA-256 mismatch for {source_label}")
        hashes[source_label] = actual_hash
        sizes[source_label] = actual_bytes
        validations.append({
            "source_label": source_label,
            "manifest": str(manifest_path.resolve()),
            "manifest_category": category,
            "manifest_label": record_label,
            "path": str(expected_path),
            "bytes": actual_bytes,
            "sha256": actual_hash,
            "validated": True,
        })
    hashes["manifest"] = sha256_file(manifest_path)
    sizes["manifest"] = manifest_path.stat().st_size
    return hashes, sizes, tuple(validations)


def validate_source_data(source_root: Path) -> SourceData:
    require(source_root.resolve().name == SOURCE_PACKAGE_ID, f"Source package must be named {SOURCE_PACKAGE_ID}")
    paths = {label: (source_root / relative).resolve() for label, relative in SOURCE_FILES.items()}
    hashes, sizes, validations = validate_source_manifest_contract(paths)
    counts = read_count_table(paths["counts"])
    depths = counts.sample_totals
    require(sum(depths) == EXPECTED_SOURCE_TOTAL, f"Expected source total {EXPECTED_SOURCE_TOTAL}; observed {sum(depths)}")
    require(min(depths) == EXPECTED_SOURCE_MIN, f"Expected source minimum {EXPECTED_SOURCE_MIN}; observed {min(depths)}")
    require(max(depths) == EXPECTED_SOURCE_MAX, f"Expected source maximum {EXPECTED_SOURCE_MAX}; observed {max(depths)}")
    require(min(depths) >= DEPTH, f"At least one sample is below {DEPTH} reads")

    taxonomy = read_keyed_tsv(paths["taxonomy"], "OTUID")
    sintax = read_headerless_keyed_tsv(paths["sintax"])
    fasta = read_fasta(paths["fasta"])
    metadata = read_keyed_tsv(paths["metadata"], "SampleID")
    require(taxonomy.header == ("OTUID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "Unexpected 16S taxonomy schema")
    require(metadata.header == ("SampleID", "site", "Group", "type", "subgroup", "subID"), "Unexpected metadata schema")
    require(taxonomy.order == counts.feature_ids, "Taxonomy IDs/order differ from count table")
    require(sintax.order == counts.feature_ids, "SINTAX IDs/order differ from count table")
    require(fasta.order == counts.feature_ids, "Representative FASTA IDs/order differ from count table")
    require(metadata.order == counts.samples, "Metadata sample IDs/order differ from count table")

    summary = json.loads(paths["summary"].read_text(encoding="utf-8-sig"))
    require(summary["feature_counts"]["analysis_nonzero_tier"] == EXPECTED_FEATURES, "Source summary feature-count mismatch")
    require(summary["sample_counts"]["retained"] == EXPECTED_SAMPLES, "Source summary sample-count mismatch")
    require(summary["read_counts"]["retained"] == EXPECTED_SOURCE_TOTAL, "Source summary read-total mismatch")
    require(summary["retained_library_size"]["minimum"] == EXPECTED_SOURCE_MIN, "Source summary minimum-depth mismatch")
    require(summary["retained_library_size"]["maximum"] == EXPECTED_SOURCE_MAX, "Source summary maximum-depth mismatch")
    require(summary["excluded_samples"] == ["SS_BS_S1", "SS_BS_S2", "SS_BS_S3"], "Source summary exclusion contract mismatch")
    return SourceData(paths, hashes, sizes, validations, counts, taxonomy, sintax, fasta, metadata)


def find_rscript(explicit: str | None) -> Path:
    if explicit:
        candidate = Path(explicit).expanduser().resolve()
        require(candidate.is_file(), f"Rscript not found: {candidate}")
        return candidate
    discovered = shutil.which("Rscript")
    if discovered:
        return Path(discovered).resolve()
    require(False, "Rscript was not found on PATH; provide --rscript explicitly")
    raise AssertionError("unreachable")


def run_command(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        list(command), check=False, text=True, encoding="utf-8", errors="replace",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise ValidationError(
            f"Command failed ({completed.returncode}): {' '.join(command)}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    return completed


def collect_software(rscript: Path) -> tuple[dict[str, str], dict[str, object]]:
    expression = (
        'if (!requireNamespace("vegan", quietly=TRUE)) quit(status=2); '
        'cat(R.version.string, "\\t", as.character(utils::packageVersion("vegan")), '
        '"\\t", utils::packageDescription("vegan")$Built, sep="")'
    )
    completed = run_command((str(rscript), "--vanilla", "-e", expression))
    fields = completed.stdout.strip().split("\t", 2)
    require(len(fields) == 3, f"Could not parse R/vegan version output: {completed.stdout!r}")
    runtime, vegan_version, vegan_built = fields
    runtime_match = re.search(r"R version (\d+)\.(\d+)\.(\d+)", runtime)
    built_match = re.search(r"R (\d+)\.(\d+)\.(\d+)", vegan_built)
    patch_mismatch = bool(
        runtime_match and built_match
        and runtime_match.group(1, 2) == built_match.group(1, 2)
        and runtime_match.group(3) != built_match.group(3)
    )
    software = {
        "python": sys.version.replace("\n", " "),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "r": runtime,
        "vegan": vegan_version,
        "vegan_built": vegan_built,
    }
    compatibility = {
        "runtime_r_version": runtime,
        "vegan_built_under": vegan_built,
        "build_runtime_patch_mismatch": patch_mismatch,
        "status": "warning_nonblocking" if patch_mismatch else "matched_or_not_applicable",
        "deterministic_oracle_recomputations": 2,
        "message": (
            "vegan was built under a different R patch release; treat as non-blocking only after two cell-identical seed-1 recomputations and all count invariants pass."
            if patch_mismatch else "No R patch-level build/runtime mismatch was detected."
        ),
    }
    return software, compatibility


def r_worker_command(rscript: Path, worker: Path, mode: str, source_counts: Path, other: Path | None = None) -> tuple[str, ...]:
    command = [
        str(rscript), "--vanilla", str(worker), mode, "--input", str(source_counts),
        "--depth", str(DEPTH), "--seed", str(SEED),
    ]
    if mode == "--run":
        require(other is not None, "R --run requires an output path")
        command.extend(("--output", str(other)))
    elif mode == "--check":
        require(other is not None, "R --check requires an expected path")
        command.extend(("--expected", str(other)))
    return tuple(command)


def write_text_no_overwrite(path: Path, text: str) -> None:
    require(not path.exists(), f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def write_json_no_overwrite(path: Path, payload: object) -> None:
    write_text_no_overwrite(path, json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def write_tsv_no_overwrite(path: Path, header: Sequence[object], rows: Iterable[Sequence[object]]) -> None:
    require(not path.exists(), f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_headerless_tsv_no_overwrite(path: Path, rows: Iterable[Sequence[object]]) -> None:
    require(not path.exists(), f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def copy_no_overwrite(source: Path, target: Path) -> None:
    require(not target.exists(), f"Refusing to overwrite: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def write_count_table(path: Path, table: CountTable, retained_indices: Sequence[int] | None = None) -> None:
    indices = tuple(range(len(table.feature_ids))) if retained_indices is None else tuple(retained_indices)
    write_tsv_no_overwrite(path, (table.id_header, *table.samples), ((table.feature_ids[i], *table.rows[i]) for i in indices))


def write_keyed_subset(path: Path, table: KeyedTable, feature_ids: Sequence[str]) -> None:
    require(all(item in table.rows for item in feature_ids), f"Missing linked ID for {path}")
    write_tsv_no_overwrite(path, table.header, (table.rows[item] for item in feature_ids))


def write_headerless_subset(path: Path, table: HeaderlessKeyedTable, feature_ids: Sequence[str]) -> None:
    require(all(item in table.rows for item in feature_ids), f"Missing linked SINTAX ID for {path}")
    write_headerless_tsv_no_overwrite(path, (table.rows[item] for item in feature_ids))


def write_fasta_subset(path: Path, fasta: FastaData, feature_ids: Sequence[str]) -> None:
    require(not path.exists(), f"Refusing to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        for feature_id in feature_ids:
            require(feature_id in fasta.sequences, f"Missing FASTA ID: {feature_id}")
            sequence = fasta.sequences[feature_id]
            handle.write(f">{feature_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def fraction_text(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator:.10f}" if denominator else "0.0000000000"


def build_sample_ledger_rows(source: SourceData, rarefied: CountTable) -> Iterator[tuple[object, ...]]:
    source_depths = source.counts.sample_totals
    rare_depths = rarefied.sample_totals
    for index, sample_id in enumerate(source.counts.samples):
        metadata = source.metadata.rows[sample_id]
        yield (
            sample_id, source_depths[index], rare_depths[index], source_depths[index] - rare_depths[index],
            fraction_text(rare_depths[index], source_depths[index]), DEPTH, SEED, *metadata[1:],
        )


def build_feature_ledger_rows(source: SourceData, rarefied: CountTable) -> Iterator[tuple[object, ...]]:
    source_totals = source.counts.feature_totals
    rare_totals = rarefied.feature_totals
    for index, feature_id in enumerate(source.counts.feature_ids):
        source_total = source_totals[index]
        rare_total = rare_totals[index]
        retained = rare_total > 0
        yield (
            feature_id, source_total, rare_total, source_total - rare_total,
            fraction_text(rare_total, source_total),
            sum(value > 0 for value in source.counts.rows[index]),
            sum(value > 0 for value in rarefied.rows[index]),
            str(retained).lower(),
            "retained" if retained else "pruned_zero_after_rarefaction",
        )


def assert_rarefaction_invariants(source: CountTable, rarefied: CountTable) -> tuple[int, ...]:
    require(rarefied.id_header == source.id_header, "Rarefied count ID header differs")
    require(rarefied.samples == source.samples, "Rarefied sample IDs/order differ")
    require(rarefied.feature_ids == source.feature_ids, "Rarefied all-feature ASV IDs/order differ")
    for feature_index, (source_row, rare_row) in enumerate(zip(source.rows, rarefied.rows)):
        for sample_index, (source_value, rare_value) in enumerate(zip(source_row, rare_row)):
            require(0 <= rare_value <= source_value, f"Invalid rarefied cell at feature {feature_index}, sample {sample_index}")
    rare_depths = rarefied.sample_totals
    require(all(value == DEPTH for value in rare_depths), "Rarefied sample total differs from 20,000")
    require(sum(rare_depths) == EXPECTED_TOTAL, "Rarefied package total differs from 1,500,000")
    return tuple(index for index, total in enumerate(rarefied.feature_totals) if total > 0)


def validate_linked_outputs(root: Path, source: SourceData, rarefied_all: CountTable, retained_indices: Sequence[int]) -> dict[str, object]:
    retained_ids = tuple(source.counts.feature_ids[index] for index in retained_indices)
    analysis = read_count_table(root / OUTPUT_FILES["counts_analysis"], expected_features=None)
    require(analysis.samples == source.counts.samples, "Analysis count sample order mismatch")
    require(analysis.feature_ids == retained_ids, "Analysis count retained-ASV order mismatch")
    require(all(total > 0 for total in analysis.feature_totals), "Analysis tier contains a zero-total ASV")
    require(analysis.rows == tuple(rarefied_all.rows[index] for index in retained_indices), "Analysis counts differ from all-feature subset")

    taxonomy = read_keyed_tsv(root / OUTPUT_FILES["taxonomy"], "OTUID")
    sintax = read_headerless_keyed_tsv(root / OUTPUT_FILES["sintax"])
    fasta = read_fasta(root / OUTPUT_FILES["fasta"])
    for label, observed in (("taxonomy", taxonomy.order), ("SINTAX", sintax.order), ("FASTA", fasta.order)):
        require(observed == retained_ids, f"{label} retained-ASV linkage/order mismatch")
    for feature_id in retained_ids:
        require(taxonomy.rows[feature_id] == source.taxonomy.rows[feature_id], f"Taxonomy cell drift: {feature_id}")
        require(sintax.rows[feature_id] == source.sintax.rows[feature_id], f"SINTAX cell drift: {feature_id}")
        require(fasta.sequences[feature_id] == source.fasta.sequences[feature_id], f"FASTA sequence drift: {feature_id}")
    require(sha256_file(root / OUTPUT_FILES["metadata"]) == source.hashes["metadata"], "Copied metadata hash mismatch")
    return {
        "retained_ids_exact": True,
        "taxonomy_exact": True,
        "sintax_exact": True,
        "fasta_sequences_exact": True,
        "metadata_exact": True,
    }


def validate_ledgers(root: Path, source: SourceData, rarefied: CountTable, retained_indices: Sequence[int]) -> None:
    sample_ledger = read_keyed_tsv(root / OUTPUT_FILES["sample_ledger"], "SampleID")
    feature_ledger = read_keyed_tsv(root / OUTPUT_FILES["feature_ledger"], "ASV_ID")
    require(sample_ledger.header == SAMPLE_LEDGER_HEADER, "Sample-ledger header mismatch")
    require(feature_ledger.header == FEATURE_LEDGER_HEADER, "Feature-ledger header mismatch")
    require(sample_ledger.order == source.counts.samples, "Sample-ledger order mismatch")
    require(feature_ledger.order == source.counts.feature_ids, "Feature-ledger order mismatch")
    expected_sample_rows = tuple(tuple(str(value) for value in row) for row in build_sample_ledger_rows(source, rarefied))
    expected_feature_rows = tuple(tuple(str(value) for value in row) for row in build_feature_ledger_rows(source, rarefied))
    require(tuple(sample_ledger.rows[item] for item in sample_ledger.order) == expected_sample_rows, "Sample-ledger values mismatch")
    require(tuple(feature_ledger.rows[item] for item in feature_ledger.order) == expected_feature_rows, "Feature-ledger values mismatch")
    retained_set = {source.counts.feature_ids[index] for index in retained_indices}
    ledger_retained = {item for item in feature_ledger.order if feature_ledger.rows[item][7] == "true"}
    require(ledger_retained == retained_set, "Feature-ledger retained status mismatch")


def write_manifest(root: Path, source: SourceData, worker: Path, orchestrator: Path) -> None:
    records: list[tuple[object, ...]] = []
    for label in ("counts", "taxonomy", "fasta", "sintax", "metadata", "summary"):
        path = source.paths[label]
        records.append(("source_input", label, str(path), path.stat().st_size, sha256_file(path)))
    path = source.paths["manifest"]
    records.append(("source_manifest", "manifest", str(path), path.stat().st_size, sha256_file(path)))
    for label, code_path in (("r_worker", worker), ("python_orchestrator", orchestrator)):
        records.append(("code", label, str(code_path.resolve()), code_path.stat().st_size, sha256_file(code_path)))
    for label, relative in OUTPUT_FILES.items():
        if label == "manifest":
            continue
        output_path = root / relative
        require(output_path.is_file(), f"Cannot manifest missing output: {output_path}")
        records.append(("output", label, relative, output_path.stat().st_size, sha256_file(output_path)))
    write_tsv_no_overwrite(root / OUTPUT_FILES["manifest"], MANIFEST_HEADER, records)


def validate_manifest(root: Path, source: SourceData, worker: Path, orchestrator: Path) -> None:
    path = root / OUTPUT_FILES["manifest"]
    require(path.is_file(), f"Missing output manifest: {path}")
    records: dict[tuple[str, str], tuple[str, int, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(tuple(reader.fieldnames or ()) == MANIFEST_HEADER, "Output-manifest header mismatch")
        for line_number, row in enumerate(reader, start=2):
            key = (row["scope"], row["label"])
            require(key not in records, f"Duplicate output-manifest record at line {line_number}")
            require(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None, f"Invalid output-manifest SHA-256 at line {line_number}")
            records[key] = (row["path"], int(row["bytes"]), row["sha256"])

    expected_keys = {
        *(("source_input", label) for label in ("counts", "taxonomy", "fasta", "sintax", "metadata", "summary")),
        ("source_manifest", "manifest"),
        ("code", "r_worker"), ("code", "python_orchestrator"),
        *(("output", label) for label in OUTPUT_FILES if label != "manifest"),
    }
    require(set(records) == expected_keys, "Output-manifest record set mismatch")
    fixed_paths: dict[tuple[str, str], Path] = {
        **{("source_input", label): source.paths[label] for label in ("counts", "taxonomy", "fasta", "sintax", "metadata", "summary")},
        ("source_manifest", "manifest"): source.paths["manifest"],
        ("code", "r_worker"): worker,
        ("code", "python_orchestrator"): orchestrator,
        **{("output", label): root / relative for label, relative in OUTPUT_FILES.items() if label != "manifest"},
    }
    for key, (recorded_path, recorded_bytes, recorded_hash) in records.items():
        actual_path = fixed_paths[key].resolve()
        if key[0] == "output":
            require(recorded_path == OUTPUT_FILES[key[1]], f"Output-manifest relative path mismatch: {key}")
        else:
            require(normalized_path(Path(recorded_path)) == normalized_path(actual_path), f"Manifest absolute path mismatch: {key}")
        require(actual_path.is_file(), f"Manifest target missing: {actual_path}")
        require(actual_path.stat().st_size == recorded_bytes, f"Manifest byte mismatch: {key}")
        require(sha256_file(actual_path) == recorded_hash, f"Manifest SHA-256 mismatch: {key}")


def read_json(path: Path) -> dict[str, object]:
    require(path.is_file(), f"Missing JSON: {path}")
    with path.open("r", encoding="utf-8-sig") as handle:
        payload = json.load(handle)
    require(isinstance(payload, dict), f"Expected JSON object: {path}")
    return payload


def validate_json_contract(root: Path, source: SourceData, retained_indices: Sequence[int]) -> None:
    config = read_json(root / OUTPUT_FILES["config"])
    summary = read_json(root / OUTPUT_FILES["summary"])
    validation = read_json(root / OUTPUT_FILES["validation"])
    require(config.get("package_id") == PACKAGE_ID, "Config package ID mismatch")
    parameters = config.get("parameters")
    require(isinstance(parameters, dict), "Config parameters missing")
    require(parameters.get("depth") == DEPTH and parameters.get("seed") == SEED, "Config depth/seed mismatch")
    require(parameters.get("algorithm") == "vegan::rrarefy" and parameters.get("replacement") is False, "Config algorithm mismatch")
    require(parameters.get("rng") == RNG_CONTRACT, "Config RNG mismatch")
    require(summary.get("status") == "validated", "Summary is not validated")
    dimensions = summary.get("dimensions")
    require(isinstance(dimensions, dict), "Summary dimensions missing")
    require(dimensions.get("samples") == EXPECTED_SAMPLES, "Summary sample count mismatch")
    require(dimensions.get("source_feature_tier") == EXPECTED_FEATURES, "Summary source-feature count mismatch")
    require(dimensions.get("all_feature_rarefied_tier") == EXPECTED_FEATURES, "Summary all-feature count mismatch")
    require(dimensions.get("analysis_feature_tier") == len(retained_indices), "Summary analysis-feature count mismatch")
    require(dimensions.get("post_rarefaction_zero_pruned") == EXPECTED_FEATURES - len(retained_indices), "Summary pruned count mismatch")
    reads = summary.get("reads")
    require(isinstance(reads, dict) and reads.get("rarefied_total") == EXPECTED_TOTAL, "Summary rarefied total mismatch")
    require(validation.get("passed") is True and validation.get("blockers") == [], "Validation JSON reports failure")
    source_files = config.get("source_files")
    require(isinstance(source_files, dict), "Config source files missing")
    for label in source.hashes:
        require(label in source_files, f"Config source record missing: {label}")
        record = source_files[label]
        require(isinstance(record, dict), f"Malformed config source record: {label}")
        require(record.get("sha256") == source.hashes[label] and record.get("bytes") == source.sizes[label], f"Config source hash/size mismatch: {label}")


def expected_output_file_set(root: Path) -> set[str]:
    return {str((root / relative).resolve()) for relative in OUTPUT_FILES.values()}


def verify_exact_file_set(root: Path) -> None:
    observed = {str(path.resolve()) for path in root.rglob("*") if path.is_file()}
    expected = expected_output_file_set(root)
    require(observed == expected, f"Unexpected package file set; missing={sorted(expected-observed)}, extra={sorted(observed-expected)}")
    forbidden = [
        path for path in root.rglob("*") if path.is_file()
        and re.search(r"(?i)(alpha|beta|divers|ordin|pcoa|nmds|permanova|anosim|figure|plot|\.png$|\.pdf$|\.svg$|\.tiff?$)", path.name)
    ]
    require(not forbidden, f"Forbidden downstream artifact(s): {forbidden}")


def deterministic_oracle(rscript: Path, worker: Path, source_counts: Path, expected: Path, runs: int = 2) -> None:
    for _ in range(runs):
        run_command(r_worker_command(rscript, worker, "--check", source_counts, expected))


def validate_package(root: Path, source: SourceData, rscript: Path, worker: Path, orchestrator: Path, run_oracle: bool) -> dict[str, object]:
    rarefied_all = read_count_table(root / OUTPUT_FILES["counts_all"])
    retained_indices = assert_rarefaction_invariants(source.counts, rarefied_all)
    linkage = validate_linked_outputs(root, source, rarefied_all, retained_indices)
    validate_ledgers(root, source, rarefied_all, retained_indices)
    validate_json_contract(root, source, retained_indices)
    validate_manifest(root, source, worker, orchestrator)
    verify_exact_file_set(root)
    if run_oracle:
        deterministic_oracle(rscript, worker, source.paths["counts"], root / OUTPUT_FILES["counts_all"], runs=2)
    pruned_ids = [source.counts.feature_ids[index] for index, total in enumerate(rarefied_all.feature_totals) if total == 0]
    return {
        "samples": EXPECTED_SAMPLES,
        "all_features": EXPECTED_FEATURES,
        "analysis_features": len(retained_indices),
        "pruned_features": len(pruned_ids),
        "pruned_feature_ids": pruned_ids,
        "rarefied_total": sum(rarefied_all.sample_totals),
        "deterministic_recomputation_mismatches": 0 if run_oracle else None,
        "deterministic_oracle_recomputations": 2 if run_oracle else 0,
        "linkage": linkage,
    }


def build_readme() -> str:
    return textwrap.dedent(
        f"""\
        # 16S rarefied analysis package at 20,000 reads/sample (v1)

        This package is built from the verified 75-sample analysis tier in
        `{SOURCE_PACKAGE_ID}`. It uses `vegan::rrarefy(t(counts), {DEPTH})`,
        without replacement, after seed `{SEED}` and the explicit RNG contract
        Mersenne-Twister / Inversion / Rejection.

        ## Feature tiers

        - `{OUTPUT_FILES['counts_all']}` preserves all {EXPECTED_FEATURES:,}
          source ASVs, including ASVs that become zero after rarefaction.
        - `{OUTPUT_FILES['counts_analysis']}` removes only post-rarefaction
          zero-total ASVs. Taxonomy, representative FASTA, and the original
          headerless four-column SINTAX records are exact keyed subsets in the
          same retained-ASV order.

        Sample and feature provenance are recorded in the two QC ledgers.
        Source/output SHA-256 values are recorded in
        `{OUTPUT_FILES['manifest']}`. Protected source-manifest records are
        revalidated before package creation and each later check.

        ## Scope boundary

        No alpha/beta diversity, hypothesis test, taxonomy summary, manuscript
        table, or figure is part of this package.

        ## Recheck

        ```powershell
        python -B revision/re_analysis/scripts/rarefy_16s_20000_v1.py --check
        ```
        """
    )


def create_stage(output_root: Path) -> Path:
    parent = output_root.parent.resolve()
    require(parent.is_dir(), f"Output parent does not exist: {parent}")
    require(output_root.resolve().name == PACKAGE_ID, f"Output package must be named {PACKAGE_ID}")
    require(not output_root.exists(), f"Refusing to overwrite existing output package: {output_root}")
    stage = Path(tempfile.mkdtemp(prefix=f".{PACKAGE_ID}.staging.", dir=parent)).resolve()
    require(stage.parent == parent and stage.name.startswith(f".{PACKAGE_ID}.staging."), "Unsafe staging path")
    return stage


def cleanup_stage(stage: Path, output_parent: Path) -> None:
    resolved = stage.resolve()
    parent = output_parent.resolve()
    require(resolved.parent == parent and resolved.name.startswith(f".{PACKAGE_ID}.staging."), f"Refusing unsafe staging cleanup: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def build_package(
    output_root: Path,
    source: SourceData,
    rscript: Path,
    worker: Path,
    orchestrator: Path,
    software: dict[str, str],
    compatibility: dict[str, object],
) -> dict[str, object]:
    stage = create_stage(output_root)
    try:
        for subdir in ("results", "qc", "provenance", "logs"):
            (stage / subdir).mkdir()
        r_run = run_command(r_worker_command(rscript, worker, "--run", source.paths["counts"], stage / OUTPUT_FILES["counts_all"]))
        write_text_no_overwrite(stage / OUTPUT_FILES["r_stdout"], r_run.stdout)
        write_text_no_overwrite(stage / OUTPUT_FILES["r_stderr"], r_run.stderr)
        rarefied_all = read_count_table(stage / OUTPUT_FILES["counts_all"])
        retained_indices = assert_rarefaction_invariants(source.counts, rarefied_all)
        retained_ids = tuple(source.counts.feature_ids[index] for index in retained_indices)
        pruned_ids = [source.counts.feature_ids[index] for index, total in enumerate(rarefied_all.feature_totals) if total == 0]

        write_count_table(stage / OUTPUT_FILES["counts_analysis"], rarefied_all, retained_indices)
        write_keyed_subset(stage / OUTPUT_FILES["taxonomy"], source.taxonomy, retained_ids)
        write_fasta_subset(stage / OUTPUT_FILES["fasta"], source.fasta, retained_ids)
        write_headerless_subset(stage / OUTPUT_FILES["sintax"], source.sintax, retained_ids)
        write_tsv_no_overwrite(stage / OUTPUT_FILES["sample_ledger"], SAMPLE_LEDGER_HEADER, build_sample_ledger_rows(source, rarefied_all))
        write_tsv_no_overwrite(stage / OUTPUT_FILES["feature_ledger"], FEATURE_LEDGER_HEADER, build_feature_ledger_rows(source, rarefied_all))
        copy_no_overwrite(source.paths["metadata"], stage / OUTPUT_FILES["metadata"])

        source_records = {
            label: {"path": str(source.paths[label]), "bytes": source.sizes[label], "sha256": source.hashes[label]}
            for label in source.hashes
        }
        config = {
            "schema_version": "1.0",
            "package_id": PACKAGE_ID,
            "marker": "16S",
            "created_utc": utc_now(),
            "parameters": {
                "depth": DEPTH, "seed": SEED, "algorithm": "vegan::rrarefy", "replacement": False,
                "input_orientation": "ASV_x_sample", "rrarefy_orientation": "sample_x_ASV", "rng": RNG_CONTRACT,
            },
            "expected": {"samples": EXPECTED_SAMPLES, "source_features": EXPECTED_FEATURES, "rarefied_total": EXPECTED_TOTAL},
            "source_root": str(source.paths["counts"].parent),
            "output_root": str(output_root.resolve()),
            "rscript": str(rscript),
            "r_worker": str(worker.resolve()),
            "software": software,
            "software_compatibility": compatibility,
            "source_files": source_records,
            "source_manifest_validation": list(source.manifest_validation),
        }
        write_json_no_overwrite(stage / OUTPUT_FILES["config"], config)

        source_total = sum(source.counts.sample_totals)
        warnings: list[str] = []
        if compatibility["build_runtime_patch_mismatch"]:
            warnings.append(str(compatibility["message"]))
        summary = {
            "schema_version": "1.0", "package_id": PACKAGE_ID, "status": "validated",
            "rarefaction": {
                "performed": True, "depth": DEPTH, "seed": SEED, "algorithm": "vegan::rrarefy",
                "without_replacement": True, "rng": RNG_CONTRACT,
            },
            "dimensions": {
                "samples": EXPECTED_SAMPLES, "source_feature_tier": EXPECTED_FEATURES,
                "all_feature_rarefied_tier": EXPECTED_FEATURES,
                "analysis_feature_tier": len(retained_ids), "post_rarefaction_zero_pruned": len(pruned_ids),
            },
            "reads": {
                "source_total": source_total, "rarefied_total": EXPECTED_TOTAL,
                "removed_reads": source_total - EXPECTED_TOTAL,
                "retention_fraction": EXPECTED_TOTAL / source_total,
                "per_sample_expected": DEPTH,
                "per_sample_min": min(rarefied_all.sample_totals),
                "per_sample_max": max(rarefied_all.sample_totals),
            },
            "linkage": {
                "counts_taxonomy_fasta_sintax_exact_order": True,
                "metadata_sample_order_exact": True,
                "source_manifest_records_validated": True,
            },
            "software": software,
            "software_compatibility": compatibility,
            "outputs": OUTPUT_FILES,
            "warnings": warnings,
        }
        write_json_no_overwrite(stage / OUTPUT_FILES["summary"], summary)
        write_text_no_overwrite(stage / OUTPUT_FILES["readme"], build_readme())

        deterministic_oracle(rscript, worker, source.paths["counts"], stage / OUTPUT_FILES["counts_all"], runs=2)
        linkage = validate_linked_outputs(stage, source, rarefied_all, retained_indices)
        validate_ledgers(stage, source, rarefied_all, retained_indices)
        validation_payload = {
            "schema_version": "1.0", "validated_utc": utc_now(), "passed": True, "blockers": [],
            "metrics": {
                "samples": EXPECTED_SAMPLES, "source_features": EXPECTED_FEATURES,
                "all_rarefied_features": EXPECTED_FEATURES,
                "analysis_features": len(retained_ids), "pruned_features": len(pruned_ids),
                "pruned_feature_ids": pruned_ids, "rarefied_total": EXPECTED_TOTAL,
                "per_cell_bound_violations": 0, "deterministic_recomputation_mismatches": 0,
                "deterministic_oracle_recomputations": 2, "linked_id_mismatches": 0,
                "source_manifest_mismatches": 0,
            },
            "linkage": linkage,
            "warnings": warnings,
        }
        write_json_no_overwrite(stage / OUTPUT_FILES["validation"], validation_payload)
        write_manifest(stage, source, worker, orchestrator)
        # The two oracle calls above already close determinism for creation.
        # Revalidate every stored artifact here without silently running two
        # additional unreported recomputations.
        validate_package(stage, source, rscript, worker, orchestrator, run_oracle=False)

        for label, expected_hash in source.hashes.items():
            require(sha256_file(source.paths[label]) == expected_hash, f"Protected source changed during build: {label}")
        require(not output_root.exists(), f"Output appeared during staging: {output_root}")
        stage.replace(output_root)
        return validation_payload["metrics"]
    except BaseException:
        cleanup_stage(stage, output_root.parent)
        raise


def validate_scoped_paths(project_root: Path, source_root: Path, output_root: Path, worker: Path) -> None:
    reanalysis_root = (project_root / "revision" / "re_analysis").resolve()
    require(source_root.resolve().parent == reanalysis_root, f"Source must be directly under {reanalysis_root}")
    require(output_root.resolve().parent == reanalysis_root, f"Output must be directly under {reanalysis_root}")
    require(source_root.resolve().name == SOURCE_PACKAGE_ID, f"Unexpected source package: {source_root}")
    require(output_root.resolve().name == PACKAGE_ID, f"Unexpected output package: {output_root}")
    require(normalized_path(source_root) != normalized_path(output_root), "Source and output paths collide")
    require(worker.resolve().parent == (reanalysis_root / "scripts"), f"R worker must remain under {reanalysis_root / 'scripts'}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    script = Path(__file__).resolve()
    project_root = script.parents[3]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--run", action="store_true", help="Create the new package; refuses an existing output path")
    modes.add_argument("--dry-run", action="store_true", help="Validate sources, R, vegan, and the frozen plan without writing output")
    modes.add_argument("--check", action="store_true", help="Read-only full package and deterministic recomputation check")
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--source-root", type=Path, default=project_root / "revision" / "re_analysis" / SOURCE_PACKAGE_ID)
    parser.add_argument("--output-root", type=Path, default=project_root / "revision" / "re_analysis" / PACKAGE_ID)
    parser.add_argument("--r-worker", type=Path, default=script.with_suffix(".R"))
    parser.add_argument("--rscript", help="Explicit path to Rscript; otherwise Rscript is resolved from PATH")
    parser.add_argument("--depth", type=int, default=DEPTH, help="Frozen at 20000")
    parser.add_argument("--seed", type=int, default=SEED, help="Frozen at 1")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args(argv)
    require(args.depth == DEPTH, "This version is frozen at depth 20000")
    require(args.seed == SEED, "This version is frozen at seed 1")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        project_root = args.project_root.resolve()
        source_root = args.source_root.resolve()
        output_root = args.output_root.resolve()
        worker = args.r_worker.resolve()
        orchestrator = Path(__file__).resolve()
        validate_scoped_paths(project_root, source_root, output_root, worker)
        require(worker.is_file(), f"Missing R worker: {worker}")
        rscript = find_rscript(args.rscript)
        source = validate_source_data(source_root)
        software, compatibility = collect_software(rscript)

        if args.dry_run:
            require(not output_root.exists(), f"Dry-run found an existing output path: {output_root}")
            r_dry = run_command(r_worker_command(rscript, worker, "--dry-run", source.paths["counts"]))
            payload = {
                "status": "dry_run_validated", "package_id": PACKAGE_ID,
                "source_features": len(source.counts.feature_ids), "samples": len(source.counts.samples),
                "source_total": sum(source.counts.sample_totals),
                "source_depth_min": min(source.counts.sample_totals),
                "source_depth_max": max(source.counts.sample_totals),
                "planned_depth": DEPTH, "planned_total": EXPECTED_TOTAL,
                "planned_reads_removed": sum(source.counts.sample_totals) - EXPECTED_TOTAL,
                "planned_retention_fraction": EXPECTED_TOTAL / sum(source.counts.sample_totals),
                "post_rarefaction_zero_features": "determined_only_by_explicit_run",
                "source_manifest_records_validated": len(source.manifest_validation),
                "r_worker_status": r_dry.stdout.strip().splitlines()[-1] if r_dry.stdout.strip() else "unknown",
                "software": software, "software_compatibility": compatibility,
                "output_exists": False,
            }
            print(json.dumps(payload, indent=2, ensure_ascii=False))
            return 0
        if args.run:
            require(not output_root.exists(), f"Refusing to overwrite existing output package: {output_root}")
            metrics = build_package(output_root, source, rscript, worker, orchestrator, software, compatibility)
            print(json.dumps({"status": "run_complete", "output_root": str(output_root), "metrics": metrics}, indent=2, ensure_ascii=False))
            return 0

        require(output_root.is_dir(), f"Missing output package for --check: {output_root}")
        metrics = validate_package(output_root, source, rscript, worker, orchestrator, run_oracle=True)
        print(json.dumps({"status": "check_passed", "output_root": str(output_root), "metrics": metrics}, indent=2, ensure_ascii=False))
        return 0
    except (ValidationError, OSError, ValueError, KeyError, TypeError, json.JSONDecodeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
