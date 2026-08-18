#!/usr/bin/env python3
"""Build and verify the versioned 75-sample 16S reanalysis input package.

The processed source files supplied by the user are treated as read-only.
Exactly SS_BS_S1, SS_BS_S2, and SS_BS_S3 are excluded.  The package retains
both an all-feature table and a zero-sum-pruned analysis table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


VERSION = "v1"
EXCLUDED_SAMPLES = ("SS_BS_S1", "SS_BS_S2", "SS_BS_S3")
EXPECTED = {
    "input_samples": 78,
    "retained_samples": 75,
    "input_features": 18169,
    "analysis_features": 17562,
    "zero_sum_features": 607,
    "original_counts": 3655885,
    "excluded_counts": 180743,
    "retained_counts": 3475142,
}

OUTPUT_FILES = {
    "counts_all": "16s_asv_counts_all_75_v1.tsv",
    "counts_analysis": "16s_asv_counts_analysis_75_v1.tsv",
    "taxonomy_all": "16s_taxonomy_all_75_v1.tsv",
    "taxonomy_analysis": "16s_taxonomy_analysis_75_v1.tsv",
    "sintax_all": "16s_sintax_all_75_v1.tsv",
    "sintax_analysis": "16s_sintax_analysis_75_v1.tsv",
    "fasta_all": "16s_rep_seqs_all_75_v1.fasta",
    "fasta_analysis": "16s_rep_seqs_analysis_75_v1.fasta",
    "metadata": "metadata_16s_75_v1.tsv",
    "filtering": "read_filtering_16s_75_v1.tsv",
    "seqkit": "seqkit_16s_75_v1.tsv",
    "excluded_samples": "excluded_samples_16s_v1.tsv",
    "zero_sum": "zero_sum_asvs_after_SS_BS_exclusion_v1.tsv",
    "summary": "package_summary_16s_75_v1.json",
    "readme": "README_16s_75_v1.md",
    "manifest": "sha256_manifest_16s_75_v1.tsv",
}


class ValidationError(RuntimeError):
    """Raised when a source or generated artifact violates the contract."""


@dataclass(frozen=True)
class CountTable:
    samples: tuple[str, ...]
    feature_ids: tuple[str, ...]
    rows: dict[str, tuple[int, ...]]


@dataclass(frozen=True)
class PackageData:
    project_root: Path
    source_paths: dict[str, Path]
    count: CountTable
    retained_samples: tuple[str, ...]
    retained_indices: tuple[int, ...]
    excluded_indices: tuple[int, ...]
    retained_rows: dict[str, tuple[int, ...]]
    retained_sample_totals: tuple[int, ...]
    analysis_ids: tuple[str, ...]
    zero_sum_ids: tuple[str, ...]
    taxonomy_header: tuple[str, ...]
    taxonomy_rows: dict[str, tuple[str, ...]]
    sintax_rows: dict[str, tuple[str, ...]]
    fasta_rows: dict[str, str]
    metadata_header: tuple[str, ...]
    metadata_rows: dict[str, tuple[str, ...]]
    filtering_header: tuple[str, ...]
    filtering_rows: dict[str, tuple[str, ...]]
    seqkit_header: tuple[str, ...]
    seqkit_rows: dict[tuple[str, str], tuple[str, ...]]
    original_total: int
    excluded_total: int
    retained_total: int


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen:
            duplicates.append(value)
        seen.add(value)
    require(not duplicates, f"Duplicate {label}: {sorted(set(duplicates))[:10]}")


def read_tsv(path: Path, *, expected_columns: int | None = None) -> tuple[tuple[str, ...], list[tuple[str, ...]]]:
    require(path.is_file(), f"Missing required file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise ValidationError(f"Empty TSV file: {path}") from exc
        rows: list[tuple[str, ...]] = []
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            if expected_columns is not None:
                require(len(row) == expected_columns, f"{path}:{line_number} has {len(row)} columns; expected {expected_columns}")
            else:
                require(len(row) == len(header), f"{path}:{line_number} has {len(row)} columns; expected {len(header)}")
            rows.append(tuple(row))
    return header, rows


def read_count_table(path: Path) -> CountTable:
    require(path.is_file(), f"Missing count table: {path}")
    samples: tuple[str, ...]
    feature_ids: list[str] = []
    rows: dict[str, tuple[int, ...]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValidationError(f"Empty count table: {path}") from exc
        require(header and header[0] in {"#OTUID", "OTUID"}, f"Unexpected count-table ID header in {path}: {header[:1]}")
        samples = tuple(header[1:])
        ensure_unique(samples, f"sample IDs in {path}")
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            require(len(row) == len(header), f"{path}:{line_number} has {len(row)} columns; expected {len(header)}")
            feature_id = row[0]
            require(feature_id not in rows, f"Duplicate feature ID {feature_id} in {path}")
            try:
                counts = tuple(int(value) for value in row[1:])
            except ValueError as exc:
                raise ValidationError(f"Non-integer count at {path}:{line_number}") from exc
            require(all(value >= 0 for value in counts), f"Negative count at {path}:{line_number}")
            feature_ids.append(feature_id)
            rows[feature_id] = counts
    return CountTable(samples=samples, feature_ids=tuple(feature_ids), rows=rows)


def read_keyed_tsv(path: Path, expected_header: tuple[str, ...] | None = None) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]]]:
    header, rows = read_tsv(path)
    if expected_header is not None:
        require(header == expected_header, f"Unexpected header in {path}: {header}")
    mapping: dict[str, tuple[str, ...]] = {}
    for row in rows:
        require(row[0] not in mapping, f"Duplicate ID {row[0]} in {path}")
        mapping[row[0]] = row
    return header, mapping


def read_sintax(path: Path) -> dict[str, tuple[str, ...]]:
    require(path.is_file(), f"Missing SINTAX file: {path}")
    mapping: dict[str, tuple[str, ...]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or all(cell == "" for cell in row):
                continue
            require(len(row) == 4, f"{path}:{line_number} has {len(row)} fields; expected 4")
            require(row[0] not in mapping, f"Duplicate SINTAX ID {row[0]} in {path}")
            mapping[row[0]] = tuple(row)
    return mapping


def read_fasta(path: Path) -> dict[str, str]:
    require(path.is_file(), f"Missing FASTA file: {path}")
    records: dict[str, str] = {}
    current_id: str | None = None
    sequence_parts: list[str] = []

    def finish_record() -> None:
        nonlocal current_id, sequence_parts
        if current_id is None:
            return
        require(current_id not in records, f"Duplicate FASTA ID {current_id} in {path}")
        sequence = "".join(sequence_parts).replace(" ", "").upper()
        require(bool(sequence), f"Empty FASTA sequence for {current_id} in {path}")
        records[current_id] = sequence

    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish_record()
                current_id = line[1:].split()[0]
                require(bool(current_id), f"Empty FASTA ID at {path}:{line_number}")
                sequence_parts = []
            else:
                require(current_id is not None, f"Sequence before FASTA header at {path}:{line_number}")
                sequence_parts.append(line)
    finish_record()
    return records


SEQKIT_FILE_RE = re.compile(r"^(?P<sample>.+)_(?P<read>[12])\.fastq\.gz$")


def read_seqkit(path: Path) -> tuple[tuple[str, ...], dict[tuple[str, str], tuple[str, ...]]]:
    require(path.is_file(), f"Missing seqkit summary: {path}")
    with path.open("r", encoding="utf-8-sig") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    require(bool(lines), f"Empty seqkit summary: {path}")
    source_header = tuple(lines[0].split())
    require(source_header == ("file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"), f"Unexpected seqkit header in {path}: {source_header}")
    output_header = ("SampleID", "read", *source_header)
    mapping: dict[tuple[str, str], tuple[str, ...]] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        require(len(fields) == 8, f"{path}:{line_number} has {len(fields)} fields; expected 8")
        filename = fields[0].replace("\\", "/").rsplit("/", 1)[-1]
        match = SEQKIT_FILE_RE.fullmatch(filename)
        require(match is not None, f"Could not parse sample/read from {path}:{line_number}: {fields[0]}")
        sample_id = match.group("sample")
        read = f"R{match.group('read')}"
        normalized = list(fields)
        for index in (3, 4, 5, 6, 7):
            normalized[index] = normalized[index].replace(",", "")
        try:
            int(normalized[3])
            int(normalized[4])
            int(normalized[5])
            float(normalized[6])
            int(normalized[7])
        except ValueError as exc:
            raise ValidationError(f"Malformed seqkit numeric field at {path}:{line_number}") from exc
        key = (sample_id, read)
        require(key not in mapping, f"Duplicate seqkit row {key} in {path}")
        mapping[key] = (sample_id, read, *normalized)
    return output_header, mapping


def source_paths(project_root: Path) -> dict[str, Path]:
    base16 = project_root / "my_results" / "result_16s"
    return {
        "count_16s": base16 / "otutab.txt",
        "taxonomy_16s": base16 / "taxonomy.txt",
        "sintax_16s": base16 / "otus.sintax",
        "fasta_16s": base16 / "otus.fa",
        "metadata_16s": base16 / "metadata1.txt",
        "filtering_16s": base16 / "raw" / "otutab_nonBac.stat",
        "seqkit_16s": base16 / "seqkit.txt",
        "count_18s_comparison": project_root / "my_results" / "result_18s" / "otutab.txt",
    }


def load_and_validate(project_root: Path) -> PackageData:
    paths = source_paths(project_root)
    for path in paths.values():
        require(path.is_file(), f"Missing required source: {path}")

    count = read_count_table(paths["count_16s"])
    require(len(count.samples) == EXPECTED["input_samples"], f"Expected 78 input samples; found {len(count.samples)}")
    require(len(count.feature_ids) == EXPECTED["input_features"], f"Expected 18,169 input features; found {len(count.feature_ids)}")
    require(set(EXCLUDED_SAMPLES).issubset(count.samples), "One or more required SS_BS samples are absent")
    require(not (set(count.samples) & {name for name in count.samples if name.startswith("SS_BS") and name not in EXCLUDED_SAMPLES}), "Unexpected SS_BS sample ID")

    retained_samples = tuple(sample for sample in count.samples if sample not in EXCLUDED_SAMPLES)
    retained_indices = tuple(index for index, sample in enumerate(count.samples) if sample not in EXCLUDED_SAMPLES)
    excluded_indices = tuple(index for index, sample in enumerate(count.samples) if sample in EXCLUDED_SAMPLES)
    require(len(retained_samples) == EXPECTED["retained_samples"], f"Expected 75 retained samples; found {len(retained_samples)}")
    require(tuple(count.samples[index] for index in excluded_indices) == EXCLUDED_SAMPLES, "SS_BS source-column order differs from the frozen exclusion contract")

    count18 = read_count_table(paths["count_18s_comparison"])
    require(count18.samples == retained_samples, "Retained 16S sample IDs/order do not exactly match the 18S count-table header")

    taxonomy_header, taxonomy_rows = read_keyed_tsv(
        paths["taxonomy_16s"],
        ("OTUID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    )
    sintax_rows = read_sintax(paths["sintax_16s"])
    fasta_rows = read_fasta(paths["fasta_16s"])
    feature_set = set(count.feature_ids)
    require(set(taxonomy_rows) == feature_set, "16S taxonomy feature IDs do not match the count table")
    require(set(sintax_rows) == feature_set, "16S SINTAX feature IDs do not match the count table")
    require(set(fasta_rows) == feature_set, "16S FASTA feature IDs do not match the count table")

    metadata_header, metadata_rows = read_keyed_tsv(
        paths["metadata_16s"],
        ("SampleID", "site", "Group", "type", "subgroup", "subID"),
    )
    require(set(metadata_rows) == set(count.samples), "16S metadata sample set does not match the count table")

    filtering_header, filtering_rows = read_keyed_tsv(
        paths["filtering_16s"],
        ("SampleID", "total_reads", "nonspecific_reads", "chloroplast_reads", "mitochondria_reads", "filtered_reads"),
    )
    require(set(filtering_rows) == set(count.samples), "Filtering-stat sample set does not match the count table")

    seqkit_header, seqkit_rows = read_seqkit(paths["seqkit_16s"])
    seqkit_samples = {sample for sample, _read in seqkit_rows}
    require(seqkit_samples == set(count.samples), "Seqkit sample set does not match the count table")
    require(len(seqkit_rows) == 2 * len(count.samples), "Seqkit summary must contain exactly R1 and R2 for every sample")
    for sample in count.samples:
        require((sample, "R1") in seqkit_rows and (sample, "R2") in seqkit_rows, f"Missing paired seqkit rows for {sample}")
        require(seqkit_rows[(sample, "R1")][5] == seqkit_rows[(sample, "R2")][5], f"R1/R2 num_seqs mismatch for {sample}")

    retained_rows: dict[str, tuple[int, ...]] = {}
    original_total = 0
    excluded_total = 0
    retained_total = 0
    retained_sample_totals = [0] * len(retained_samples)
    zero_sum_ids: list[str] = []
    analysis_ids: list[str] = []
    for feature_id in count.feature_ids:
        values = count.rows[feature_id]
        original_feature_total = sum(values)
        excluded_feature_total = sum(values[index] for index in excluded_indices)
        retained = tuple(values[index] for index in retained_indices)
        retained_feature_total = sum(retained)
        retained_rows[feature_id] = retained
        original_total += original_feature_total
        excluded_total += excluded_feature_total
        retained_total += retained_feature_total
        for index, value in enumerate(retained):
            retained_sample_totals[index] += value
        if retained_feature_total == 0:
            zero_sum_ids.append(feature_id)
        else:
            analysis_ids.append(feature_id)

    require(original_total == EXPECTED["original_counts"], f"Original count total drift: {original_total}")
    require(excluded_total == EXPECTED["excluded_counts"], f"Excluded count total drift: {excluded_total}")
    require(retained_total == EXPECTED["retained_counts"], f"Retained count total drift: {retained_total}")
    require(original_total == excluded_total + retained_total, "Original counts do not reconcile to excluded + retained")
    require(len(zero_sum_ids) == EXPECTED["zero_sum_features"], f"Expected 607 zero-sum features; found {len(zero_sum_ids)}")
    require(len(analysis_ids) == EXPECTED["analysis_features"], f"Expected 17,562 analysis features; found {len(analysis_ids)}")
    require(all(value > 0 for value in retained_sample_totals), "At least one retained sample has zero counts")

    for sample, count_sum in zip(retained_samples, retained_sample_totals):
        try:
            filtered_reads = int(filtering_rows[sample][5])
        except ValueError as exc:
            raise ValidationError(f"Non-integer filtered_reads for {sample}") from exc
        require(filtered_reads == count_sum, f"Count-table total ({count_sum}) != filtered_reads ({filtered_reads}) for {sample}")

    return PackageData(
        project_root=project_root,
        source_paths=paths,
        count=count,
        retained_samples=retained_samples,
        retained_indices=retained_indices,
        excluded_indices=excluded_indices,
        retained_rows=retained_rows,
        retained_sample_totals=tuple(retained_sample_totals),
        analysis_ids=tuple(analysis_ids),
        zero_sum_ids=tuple(zero_sum_ids),
        taxonomy_header=taxonomy_header,
        taxonomy_rows=taxonomy_rows,
        sintax_rows=sintax_rows,
        fasta_rows=fasta_rows,
        metadata_header=metadata_header,
        metadata_rows=metadata_rows,
        filtering_header=filtering_header,
        filtering_rows=filtering_rows,
        seqkit_header=seqkit_header,
        seqkit_rows=seqkit_rows,
        original_total=original_total,
        excluded_total=excluded_total,
        retained_total=retained_total,
    )


def write_tsv(path: Path, header: Sequence[str] | None, rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header is not None:
            writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def write_fasta(path: Path, feature_ids: Sequence[str], sequences: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for feature_id in feature_ids:
            handle.write(f">{feature_id}\n")
            for line in textwrap.wrap(sequences[feature_id], width=80):
                handle.write(f"{line}\n")


def summary_dict(data: PackageData) -> dict[str, object]:
    group_counts: dict[str, int] = {}
    type_counts: dict[str, int] = {}
    for sample in data.retained_samples:
        metadata = data.metadata_rows[sample]
        group_counts[metadata[2]] = group_counts.get(metadata[2], 0) + 1
        type_counts[metadata[3]] = type_counts.get(metadata[3], 0) + 1
    original_global_zero = sum(1 for feature_id in data.zero_sum_ids if sum(data.count.rows[feature_id]) == 0)
    return {
        "package_version": VERSION,
        "marker": "16S rRNA gene",
        "transformation": "exclude SS_BS samples; preserve all features; create a retained-sample nonzero analysis tier",
        "excluded_samples": list(EXCLUDED_SAMPLES),
        "sample_counts": {"input": len(data.count.samples), "retained": len(data.retained_samples), "excluded": len(EXCLUDED_SAMPLES)},
        "feature_counts": {
            "all_feature_tier": len(data.count.feature_ids),
            "analysis_nonzero_tier": len(data.analysis_ids),
            "zero_sum_after_exclusion": len(data.zero_sum_ids),
            "SS_BS_only": len(data.zero_sum_ids) - original_global_zero,
            "original_global_zero": original_global_zero,
        },
        "read_counts": {"original": data.original_total, "excluded": data.excluded_total, "retained": data.retained_total},
        "retained_library_size": {
            "minimum": min(data.retained_sample_totals),
            "maximum": max(data.retained_sample_totals),
            "total": sum(data.retained_sample_totals),
            "zero_count_samples": sum(value == 0 for value in data.retained_sample_totals),
        },
        "design": {"Group": dict(sorted(group_counts.items())), "type": dict(sorted(type_counts.items()))},
        "sample_order_matches_existing_18S_count_header": True,
        "metadata_reordered_by_SampleID_to_count_header": True,
        "existing_18S_feature_provenance_validated": False,
        "notes": [
            "The all-feature tier preserves 607 zero-sum rows for provenance.",
            "Use the analysis_nonzero tier for downstream methods that require positive feature totals.",
            "The existing 18S table is used only to freeze the 75-sample ID/order contract; its feature provenance is outside this package.",
        ],
    }


def write_package(build_dir: Path, data: PackageData) -> None:
    all_ids = data.count.feature_ids
    analysis_ids = data.analysis_ids

    write_tsv(
        build_dir / OUTPUT_FILES["counts_all"],
        ("#OTUID", *data.retained_samples),
        ((feature_id, *data.retained_rows[feature_id]) for feature_id in all_ids),
    )
    write_tsv(
        build_dir / OUTPUT_FILES["counts_analysis"],
        ("#OTUID", *data.retained_samples),
        ((feature_id, *data.retained_rows[feature_id]) for feature_id in analysis_ids),
    )
    write_tsv(
        build_dir / OUTPUT_FILES["taxonomy_all"],
        data.taxonomy_header,
        (data.taxonomy_rows[feature_id] for feature_id in all_ids),
    )
    write_tsv(
        build_dir / OUTPUT_FILES["taxonomy_analysis"],
        data.taxonomy_header,
        (data.taxonomy_rows[feature_id] for feature_id in analysis_ids),
    )
    write_tsv(build_dir / OUTPUT_FILES["sintax_all"], None, (data.sintax_rows[feature_id] for feature_id in all_ids))
    write_tsv(build_dir / OUTPUT_FILES["sintax_analysis"], None, (data.sintax_rows[feature_id] for feature_id in analysis_ids))
    write_fasta(build_dir / OUTPUT_FILES["fasta_all"], all_ids, data.fasta_rows)
    write_fasta(build_dir / OUTPUT_FILES["fasta_analysis"], analysis_ids, data.fasta_rows)

    write_tsv(
        build_dir / OUTPUT_FILES["metadata"],
        data.metadata_header,
        (data.metadata_rows[sample] for sample in data.retained_samples),
    )
    write_tsv(
        build_dir / OUTPUT_FILES["filtering"],
        data.filtering_header,
        (data.filtering_rows[sample] for sample in data.retained_samples),
    )
    write_tsv(
        build_dir / OUTPUT_FILES["seqkit"],
        data.seqkit_header,
        (data.seqkit_rows[(sample, read)] for sample in data.retained_samples for read in ("R1", "R2")),
    )

    excluded_rows = []
    for sample in EXCLUDED_SAMPLES:
        sample_index = data.count.samples.index(sample)
        counts = sum(data.count.rows[feature_id][sample_index] for feature_id in data.count.feature_ids)
        filtered_reads = int(data.filtering_rows[sample][5])
        excluded_rows.append((*data.metadata_rows[sample], counts, filtered_reads, "excluded_exact_SS_BS_sample"))
    write_tsv(
        build_dir / OUTPUT_FILES["excluded_samples"],
        (*data.metadata_header, "count_table_total", "filter_stat_filtered_reads", "reason"),
        excluded_rows,
    )

    zero_rows = []
    for feature_id in data.zero_sum_ids:
        values = data.count.rows[feature_id]
        original_total = sum(values)
        excluded_total = sum(values[index] for index in data.excluded_indices)
        classification = "original_global_zero" if original_total == 0 else "SS_BS_only"
        zero_rows.append((feature_id, classification, original_total, excluded_total, 0))
    write_tsv(
        build_dir / OUTPUT_FILES["zero_sum"],
        ("OTUID", "classification", "original_total", "excluded_SS_BS_total", "retained_total"),
        zero_rows,
    )

    summary = summary_dict(data)
    with (build_dir / OUTPUT_FILES["summary"]).open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    readme = f"""# 16S 75-sample reanalysis input package ({VERSION})

This package was derived from the canonical final 16S ASV outputs under
`my_results/result_16s` without modifying those sources.

## Frozen transformation

- Excluded samples: `{EXCLUDED_SAMPLES[0]}`, `{EXCLUDED_SAMPLES[1]}`, `{EXCLUDED_SAMPLES[2]}`.
- Retained samples: {len(data.retained_samples)}; their IDs and order exactly match the existing 18S count-table header.
- Original counts: {data.original_total:,}; excluded: {data.excluded_total:,}; retained: {data.retained_total:,}.
- All-feature tier: {len(data.count.feature_ids):,} ASVs (including zero-sum rows for provenance).
- Analysis tier: {len(data.analysis_ids):,} ASVs after removing {len(data.zero_sum_ids):,} rows with a zero total across the retained samples.
- Metadata and read summaries were joined and reordered by explicit `SampleID`.

## Intended use

Use `16s_asv_counts_analysis_75_v1.tsv` and its matching taxonomy/SINTAX/FASTA files
for downstream analysis. The all-feature files preserve complete lineage from the
prior final 16S table. Do not use old rarefied, mean-aggregated, ordination,
differential-abundance, machine-learning, neutral-model, or network outputs as
75-sample results; regenerate them from this package.

## 18S caution

The existing 18S count table was consulted only to confirm the common 75-sample
ID/order contract. This package does not validate or copy existing 18S features.

## Reproduction and validation

```powershell
python -B revision\\re_analysis\\scripts\\prepare_16s_75_inputs_v1.py --dry-run
python -B revision\\re_analysis\\scripts\\prepare_16s_75_inputs_v1.py --check --output-dir revision\\re_analysis\\01_input_16s_75_v1
```

The SHA-256 manifest records every protected source and every generated artifact
except the manifest itself. `--check` recomputes hashes and validates linked IDs,
sample order, dimensions, count totals, ledgers, and metadata/read-stat alignment.
"""
    (build_dir / OUTPUT_FILES["readme"]).write_text(readme, encoding="utf-8", newline="\n")

    manifest_rows: list[tuple[object, ...]] = []
    for label, path in sorted(data.source_paths.items()):
        manifest_rows.append(("source", label, path.relative_to(data.project_root).as_posix(), path.stat().st_size, sha256_file(path)))
    for key, filename in sorted(OUTPUT_FILES.items()):
        if key == "manifest":
            continue
        path = build_dir / filename
        manifest_rows.append(("output", key, filename, path.stat().st_size, sha256_file(path)))
    write_tsv(
        build_dir / OUTPUT_FILES["manifest"],
        ("category", "label", "path", "bytes", "sha256"),
        manifest_rows,
    )


def validate_manifest(output_dir: Path, data: PackageData) -> None:
    manifest_path = output_dir / OUTPUT_FILES["manifest"]
    header, rows = read_tsv(manifest_path)
    require(header == ("category", "label", "path", "bytes", "sha256"), f"Unexpected manifest header: {header}")
    seen: set[tuple[str, str]] = set()
    expected_keys = {("source", key) for key in data.source_paths} | {("output", key) for key in OUTPUT_FILES if key != "manifest"}
    for row in rows:
        category, label, relative_path, size_text, expected_hash = row
        key = (category, label)
        require(key not in seen, f"Duplicate manifest row: {key}")
        seen.add(key)
        if category == "source":
            require(label in data.source_paths, f"Unknown source manifest label: {label}")
            path = data.project_root / Path(relative_path)
            require(path.resolve() == data.source_paths[label].resolve(), f"Source manifest path mismatch for {label}")
        elif category == "output":
            require(label in OUTPUT_FILES and label != "manifest", f"Unknown output manifest label: {label}")
            require(relative_path == OUTPUT_FILES[label], f"Output manifest path mismatch for {label}")
            path = output_dir / relative_path
        else:
            raise ValidationError(f"Unknown manifest category: {category}")
        require(path.is_file(), f"Manifest target missing: {path}")
        require(path.stat().st_size == int(size_text), f"Size drift for {path}")
        require(sha256_file(path) == expected_hash, f"SHA-256 drift for {path}")
    require(seen == expected_keys, f"Manifest row-set mismatch; missing={sorted(expected_keys - seen)}, extra={sorted(seen - expected_keys)}")


def keyed_output_ids(path: Path, *, headered: bool) -> tuple[str, ...]:
    ids: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if headered:
            try:
                next(reader)
            except StopIteration as exc:
                raise ValidationError(f"Empty output file: {path}") from exc
        for line_number, row in enumerate(reader, start=2 if headered else 1):
            if not row or all(cell == "" for cell in row):
                continue
            require(bool(row[0]), f"Empty ID in {path}:{line_number}")
            ids.append(row[0])
    ensure_unique(ids, f"feature IDs in {path}")
    return tuple(ids)


def check_package(output_dir: Path, data: PackageData) -> dict[str, object]:
    require(output_dir.is_dir(), f"Output directory does not exist: {output_dir}")
    expected_names = set(OUTPUT_FILES.values())
    actual_names = {path.name for path in output_dir.iterdir() if path.is_file()}
    require(actual_names == expected_names, f"Output file set mismatch; missing={sorted(expected_names - actual_names)}, extra={sorted(actual_names - expected_names)}")
    require(not any(path.is_dir() for path in output_dir.iterdir()), "Unexpected subdirectory in output package")
    validate_manifest(output_dir, data)

    all_counts = read_count_table(output_dir / OUTPUT_FILES["counts_all"])
    analysis_counts = read_count_table(output_dir / OUTPUT_FILES["counts_analysis"])
    require(all_counts.samples == data.retained_samples, "All-feature count sample order mismatch")
    require(analysis_counts.samples == data.retained_samples, "Analysis count sample order mismatch")
    require(all_counts.feature_ids == data.count.feature_ids, "All-feature count ID/order mismatch")
    require(analysis_counts.feature_ids == data.analysis_ids, "Analysis count ID/order mismatch")
    require(sum(sum(values) for values in all_counts.rows.values()) == data.retained_total, "All-feature count total mismatch")
    require(sum(sum(values) for values in analysis_counts.rows.values()) == data.retained_total, "Analysis count total mismatch")
    require(all(sum(values) > 0 for values in analysis_counts.rows.values()), "Analysis count table contains a zero-sum feature")

    expected_tiers = (("all", data.count.feature_ids), ("analysis", data.analysis_ids))
    for tier, expected_ids in expected_tiers:
        taxonomy_ids = keyed_output_ids(output_dir / OUTPUT_FILES[f"taxonomy_{tier}"], headered=True)
        sintax_ids = keyed_output_ids(output_dir / OUTPUT_FILES[f"sintax_{tier}"], headered=False)
        fasta_ids = tuple(read_fasta(output_dir / OUTPUT_FILES[f"fasta_{tier}"]).keys())
        require(taxonomy_ids == expected_ids, f"{tier} taxonomy ID/order mismatch")
        require(sintax_ids == expected_ids, f"{tier} SINTAX ID/order mismatch")
        require(fasta_ids == expected_ids, f"{tier} FASTA ID/order mismatch")

    metadata_header, metadata_rows = read_tsv(output_dir / OUTPUT_FILES["metadata"])
    require(metadata_header == data.metadata_header, "Generated metadata header mismatch")
    require(tuple(row[0] for row in metadata_rows) == data.retained_samples, "Generated metadata sample order mismatch")
    filtering_header, filtering_rows = read_tsv(output_dir / OUTPUT_FILES["filtering"])
    require(filtering_header == data.filtering_header, "Generated filtering header mismatch")
    require(tuple(row[0] for row in filtering_rows) == data.retained_samples, "Generated filtering sample order mismatch")
    for row, expected_total in zip(filtering_rows, data.retained_sample_totals):
        require(int(row[5]) == expected_total, f"Generated filtering/count mismatch for {row[0]}")

    seqkit_header, seqkit_rows = read_tsv(output_dir / OUTPUT_FILES["seqkit"])
    require(seqkit_header == data.seqkit_header, "Generated seqkit header mismatch")
    expected_seqkit_keys = tuple((sample, read) for sample in data.retained_samples for read in ("R1", "R2"))
    require(tuple((row[0], row[1]) for row in seqkit_rows) == expected_seqkit_keys, "Generated seqkit sample/read order mismatch")

    excluded_header, excluded_rows = read_tsv(output_dir / OUTPUT_FILES["excluded_samples"])
    require(excluded_header[-3:] == ("count_table_total", "filter_stat_filtered_reads", "reason"), "Excluded-sample ledger header mismatch")
    require(tuple(row[0] for row in excluded_rows) == EXCLUDED_SAMPLES, "Excluded-sample ledger IDs/order mismatch")
    require(sum(int(row[-3]) for row in excluded_rows) == data.excluded_total, "Excluded-sample ledger count total mismatch")
    require(all(row[-3] == row[-2] for row in excluded_rows), "Excluded count total and filter stat differ")

    zero_header, zero_rows = read_tsv(output_dir / OUTPUT_FILES["zero_sum"])
    require(zero_header == ("OTUID", "classification", "original_total", "excluded_SS_BS_total", "retained_total"), "Zero-sum ledger header mismatch")
    require(tuple(row[0] for row in zero_rows) == data.zero_sum_ids, "Zero-sum ledger ID/order mismatch")
    require(all(int(row[4]) == 0 for row in zero_rows), "Zero-sum ledger contains a nonzero retained total")
    require(sum(row[1] == "SS_BS_only" for row in zero_rows) == 606, "Expected 606 SS_BS-only ASVs")
    require(sum(row[1] == "original_global_zero" for row in zero_rows) == 1, "Expected one original global-zero ASV")

    with (output_dir / OUTPUT_FILES["summary"]).open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    require(summary == summary_dict(data), "JSON package summary differs from recomputed source summary")
    return summary


def resolve_project_root(script_path: Path, supplied: str | None) -> Path:
    if supplied is not None:
        root = Path(supplied).expanduser().resolve()
    else:
        root = script_path.resolve().parents[3]
    require((root / "my_results" / "result_16s").is_dir(), f"Could not locate 16S source tree under project root: {root}")
    return root


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare or verify the versioned 75-sample 16S reanalysis input package.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--project-root", help="2024-amplicon project root; inferred from this script when omitted")
    parser.add_argument("--output-dir", help="Output package directory; defaults under revision/re_analysis")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="Validate sources and print the planned summary without writing")
    mode.add_argument("--check", action="store_true", help="Read-only validation of an existing generated package")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    project_root = resolve_project_root(Path(__file__), args.project_root)
    output_dir = Path(args.output_dir).expanduser() if args.output_dir else project_root / "revision" / "re_analysis" / "01_input_16s_75_v1"
    if not output_dir.is_absolute():
        output_dir = (Path.cwd() / output_dir).resolve()
    else:
        output_dir = output_dir.resolve()
    data = load_and_validate(project_root)

    if args.dry_run:
        print(json.dumps(summary_dict(data), ensure_ascii=False, indent=2, sort_keys=True))
        print(f"DRY RUN PASSED: no files written; planned output={output_dir}")
        return 0

    if args.check:
        summary = check_package(output_dir, data)
        print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
        print(f"CHECK PASSED: {output_dir}")
        return 0

    require(not output_dir.exists(), f"Refusing to overwrite existing output path: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    build_path = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}_build_", dir=output_dir.parent))
    try:
        write_package(build_path, data)
        check_package(build_path, data)
        build_path.rename(output_dir)
    except Exception:
        if build_path.exists() and build_path.parent.resolve() == output_dir.parent.resolve() and build_path.name.startswith(f".{output_dir.name}_build_"):
            shutil.rmtree(build_path)
        raise
    print(f"PACKAGE CREATED AND VERIFIED: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
