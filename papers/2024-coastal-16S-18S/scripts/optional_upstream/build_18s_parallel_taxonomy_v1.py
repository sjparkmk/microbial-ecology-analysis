#!/usr/bin/env python3
"""Build a PR2-primary, SILVA-v123-comparator 18S taxonomy package.

This helper starts from the completed upstream 0.97 mapping outputs.  It does
not run SINTAX, rarefaction, diversity analysis, or plotting.  PR2 and SILVA
calls are kept as independent evidence surfaces and are joined only by the
explicit feature identifier.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import re
import shutil
import sys
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


VERSION = "v1"
UNASSIGNED = "Unassigned"
CANONICAL_RANKS = (
    "Domain",
    "Kingdom",
    "Division_or_Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
)
PR2_TAGS = ("k", "d", "p", "c", "o", "f", "g", "s")
SILVA_TAGS = ("d", "k", "p", "c", "o", "f", "g", "s")
PR2_SOURCE_NAMES = {
    "k": "domain",
    "d": "supergroup",
    "p": "division_subdivision_composite",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
SILVA_SOURCE_NAMES = {
    "d": "domain",
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
PR2_COMPARTMENT_CODES = frozenset({"plas", "apic", "chro", "mito", "nucl"})
PLASTID_CODES = frozenset({"plas", "apic", "chro"})
OUTPUT_FILES = {
    "taxonomy_all": "results/18s_taxonomy_pr2_primary_all_v1.tsv",
    "taxonomy_analysis": "results/18s_taxonomy_pr2_primary_analysis_v1.tsv",
    "parallel": "results/18s_taxonomy_parallel_pr2_silva123_all_v1.tsv",
    "counts_analysis": "results/18s_asv_counts_pr2_analysis_v1.tsv",
    "fasta_analysis": "results/18s_rep_seqs_pr2_analysis_v1.fasta",
    "ledger": "qc/18s_asv_flags_ledger_v1.tsv",
    "sample_depths": "qc/18s_sample_depths_pr2_analysis_v1.tsv",
    "rarefaction_candidates": "qc/18s_rarefaction_candidates_v1.tsv",
    "summary": "provenance/package_summary_18s_parallel_taxonomy_v1.json",
    "readme": "provenance/README_18s_parallel_taxonomy_v1.md",
    "manifest": "provenance/sha256_manifest_18s_parallel_taxonomy_v1.tsv",
}


class ValidationError(RuntimeError):
    """Raised when an input or generated artifact violates the contract."""


@dataclass(frozen=True)
class CountTable:
    id_header: str
    samples: tuple[str, ...]
    feature_ids: tuple[str, ...]
    rows: dict[str, tuple[int, ...]]


@dataclass(frozen=True)
class RankCall:
    tag: str
    raw_values: tuple[str, ...]
    confidences: tuple[float, ...]
    accepted_values: tuple[str, ...]
    status: str

    @property
    def accepted_value(self) -> str | None:
        if self.status == "accepted_unique":
            return self.accepted_values[0]
        return None


@dataclass(frozen=True)
class SintaxAssignment:
    feature_id: str
    raw_assignment: str
    strand: str
    cutoff_assignment: str
    calls: dict[str, RankCall]
    unexpected_tags: tuple[str, ...]
    cutoff_matches_raw_confidence: bool


@dataclass(frozen=True)
class Pr2Crosswalk:
    composite_to_ranks: dict[str, tuple[str, str]]
    record_count: int
    unique_lineage_count: int
    unique_composite_count: int


@dataclass(frozen=True)
class FeatureResult:
    feature_id: str
    pr2_canonical: tuple[str, ...]
    silva_canonical: tuple[str, ...]
    pr2_subdivision: str
    pr2_compartments: tuple[str, ...]
    flags: dict[str, bool]
    analysis_keep: bool
    exclusion_reason: str
    agreement_level: str
    agreement_status: str
    agreement_taxon: str
    pr2_agreement_rank: str
    silva_agreement_rank: str


@dataclass(frozen=True)
class PackageData:
    source_paths: dict[str, Path]
    source_hashes: dict[str, str]
    count: CountTable
    fasta: dict[str, str]
    metadata_header: tuple[str, ...]
    metadata_rows: dict[str, tuple[str, ...]]
    frozen_order: tuple[str, ...]
    pr2: dict[str, SintaxAssignment]
    silva: dict[str, SintaxAssignment]
    crosswalk: Pr2Crosswalk
    results: dict[str, FeatureResult]
    analysis_ids: tuple[str, ...]
    cutoff: float
    feature_prefix: str


RAW_CALL_RE = re.compile(
    r"^(?P<tag>[A-Za-z]):(?P<value>.*)\((?P<confidence>(?:0(?:\.\d+)?)|(?:1(?:\.0+)?))\)$"
)
PR2_SUFFIX_RE = re.compile(r":(?P<code>plas|apic|chro|mito|nucl)$", re.IGNORECASE)
TRAILING_ANNOTATION_RE = re.compile(r"_\([^()]*\)$")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def format_confidence(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def ensure_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    require(not duplicates, f"Duplicate {label}: {sorted(duplicates)[:10]}")


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[tuple[str, ...]]]:
    require(path.is_file(), f"Missing TSV: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise ValidationError(f"Empty TSV: {path}") from exc
        require(bool(header), f"Empty TSV header: {path}")
        rows: list[tuple[str, ...]] = []
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            require(
                len(row) == len(header),
                f"{path}:{line_number} has {len(row)} columns; expected {len(header)}",
            )
            rows.append(tuple(row))
    return header, rows


def read_count_table(path: Path) -> CountTable:
    require(path.is_file(), f"Missing 0.97 count table: {path}")
    feature_ids: list[str] = []
    rows: dict[str, tuple[int, ...]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValidationError(f"Empty count table: {path}") from exc
        require(len(header) >= 2 and bool(header[0]), f"Malformed count-table header: {path}")
        samples = tuple(header[1:])
        require(all(samples), f"Empty sample ID in count-table header: {path}")
        ensure_unique(samples, f"sample IDs in {path}")
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            require(
                len(row) == len(header),
                f"{path}:{line_number} has {len(row)} columns; expected {len(header)}",
            )
            feature_id = row[0]
            require(bool(feature_id), f"Empty feature ID at {path}:{line_number}")
            require(feature_id not in rows, f"Duplicate feature ID {feature_id} in {path}")
            try:
                counts = tuple(int(value) for value in row[1:])
            except ValueError as exc:
                raise ValidationError(f"Non-integer count at {path}:{line_number}") from exc
            require(all(value >= 0 for value in counts), f"Negative count at {path}:{line_number}")
            feature_ids.append(feature_id)
            rows[feature_id] = counts
    require(bool(feature_ids), f"Count table has no feature rows: {path}")
    return CountTable(
        id_header=header[0],
        samples=samples,
        feature_ids=tuple(feature_ids),
        rows=rows,
    )


def read_fasta(path: Path) -> dict[str, str]:
    require(path.is_file(), f"Missing representative FASTA: {path}")
    records: dict[str, str] = {}
    current_id: str | None = None
    parts: list[str] = []

    def finish() -> None:
        nonlocal current_id, parts
        if current_id is None:
            return
        require(current_id not in records, f"Duplicate FASTA ID {current_id} in {path}")
        sequence = "".join(parts).replace(" ", "").upper()
        require(bool(sequence), f"Empty FASTA sequence for {current_id} in {path}")
        require(
            re.fullmatch(r"[ACGTUNRYKMSWBDHV.-]+", sequence) is not None,
            f"Invalid nucleotide symbol in FASTA sequence {current_id} in {path}",
        )
        records[current_id] = sequence

    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                current_id = line[1:].split()[0]
                require(bool(current_id), f"Empty FASTA ID at {path}:{line_number}")
                parts = []
            else:
                require(current_id is not None, f"Sequence before FASTA header at {path}:{line_number}")
                parts.append(line)
    finish()
    require(bool(records), f"Representative FASTA has no records: {path}")
    return records


def read_metadata(path: Path) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]], tuple[str, ...]]:
    header, rows = read_tsv(path)
    require(
        header[0].lstrip("#").lower() in {"sampleid", "sample_id", "sample"},
        f"Metadata first column must identify SampleID: {path} ({header[0]!r})",
    )
    mapping: dict[str, tuple[str, ...]] = {}
    order: list[str] = []
    for row in rows:
        sample = row[0]
        require(bool(sample), f"Empty SampleID in metadata: {path}")
        require(sample not in mapping, f"Duplicate metadata SampleID {sample} in {path}")
        mapping[sample] = row
        order.append(sample)
    require(bool(order), f"Metadata contains no samples: {path}")
    return header, mapping, tuple(order)


def read_frozen_order(path: Path) -> tuple[str, ...]:
    """Read a count header, metadata-like TSV, or one-ID-per-line order file."""

    require(path.is_file(), f"Missing frozen-order file: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    require(bool(lines), f"Empty frozen-order file: {path}")
    first = lines[0].split("\t")
    first_key = first[0].lstrip("#").lower()
    if len(first) > 1 and first_key in {"otuid", "otu id", "featureid", "feature_id"}:
        order = tuple(first[1:])
    elif first_key in {"sampleid", "sample_id", "sample"}:
        order = tuple(line.split("\t", 1)[0] for line in lines[1:])
    elif len(first) == 1:
        order = tuple(line.split("\t", 1)[0] for line in lines)
    else:
        raise ValidationError(
            f"Could not identify frozen-order format in {path}; use a count header, metadata TSV, or one ID per line"
        )
    require(bool(order) and all(order), f"Frozen order has an empty or missing SampleID: {path}")
    ensure_unique(order, f"sample IDs in frozen order {path}")
    return order


def _parse_raw_calls(
    text: str,
    *,
    path: Path,
    line_number: int,
) -> list[tuple[str, str, float]]:
    calls: list[tuple[str, str, float]] = []
    if not text:
        return calls
    for token in text.split(","):
        match = RAW_CALL_RE.fullmatch(token)
        require(match is not None, f"Malformed raw SINTAX token at {path}:{line_number}: {token!r}")
        assert match is not None
        tag = match.group("tag").lower()
        value = match.group("value")
        confidence = float(match.group("confidence"))
        require(bool(value), f"Empty SINTAX taxon at {path}:{line_number}: {token!r}")
        require(0.0 <= confidence <= 1.0, f"SINTAX confidence outside [0,1] at {path}:{line_number}")
        calls.append((tag, value, confidence))
    return calls


def _parse_cutoff_calls(
    text: str,
    *,
    path: Path,
    line_number: int,
) -> list[tuple[str, str]]:
    calls: list[tuple[str, str]] = []
    if not text:
        return calls
    for token in text.split(","):
        require(":" in token, f"Malformed cutoff SINTAX token at {path}:{line_number}: {token!r}")
        tag, value = token.split(":", 1)
        require(len(tag) == 1 and tag.isalpha(), f"Malformed SINTAX rank tag at {path}:{line_number}: {tag!r}")
        require(bool(value), f"Empty cutoff SINTAX taxon at {path}:{line_number}: {token!r}")
        calls.append((tag.lower(), value))
    return calls


def _rank_call(
    tag: str,
    raw_calls: Sequence[tuple[str, str, float]],
    cutoff: float,
) -> RankCall:
    selected = [(value, confidence) for call_tag, value, confidence in raw_calls if call_tag == tag]
    raw_values = tuple(value for value, _ in selected)
    confidences = tuple(confidence for _, confidence in selected)
    accepted_ordered = [value for value, confidence in selected if confidence + 1e-12 >= cutoff]
    accepted_values = tuple(dict.fromkeys(accepted_ordered))
    distinct_raw = tuple(dict.fromkeys(raw_values))
    if not raw_values:
        status = "missing"
    elif len(distinct_raw) > 1:
        status = "ambiguous" if accepted_values else "ambiguous_below_cutoff"
    elif not accepted_values:
        status = "below_cutoff"
    else:
        status = "accepted_unique"
    return RankCall(
        tag=tag,
        raw_values=raw_values,
        confidences=confidences,
        accepted_values=accepted_values,
        status=status,
    )


def read_sintax(
    path: Path,
    *,
    source: str,
    expected_tags: Sequence[str],
    cutoff: float,
) -> dict[str, SintaxAssignment]:
    require(path.is_file(), f"Missing {source} SINTAX table: {path}")
    expected = set(expected_tags)
    assignments: dict[str, SintaxAssignment] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or all(cell == "" for cell in row):
                continue
            if line_number == 1 and row[0].lstrip("#").lower() in {
                "asv_id",
                "feature_id",
                "query",
                "otuid",
            }:
                require(len(row) == 4, f"Unexpected {source} SINTAX header width in {path}")
                continue
            require(
                len(row) == 4,
                f"{path}:{line_number} has {len(row)} fields; expected raw 4-column VSEARCH SINTAX output",
            )
            feature_id, raw_text, strand, cutoff_text = row
            require(bool(feature_id), f"Empty SINTAX feature ID at {path}:{line_number}")
            require(feature_id not in assignments, f"Duplicate SINTAX ID {feature_id} in {path}")
            raw_calls = _parse_raw_calls(raw_text, path=path, line_number=line_number)
            cutoff_calls = _parse_cutoff_calls(cutoff_text, path=path, line_number=line_number)
            unexpected = tuple(
                dict.fromkeys(tag for tag, _, _ in raw_calls if tag not in expected)
            )
            unexpected_cutoff = tuple(
                dict.fromkeys(tag for tag, _ in cutoff_calls if tag not in expected)
            )
            require(
                not unexpected and not unexpected_cutoff,
                f"Unexpected {source} rank tag(s) at {path}:{line_number}: "
                f"{sorted(set(unexpected) | set(unexpected_cutoff))}",
            )
            expected_cutoff_calls = [
                (tag, value)
                for tag, value, confidence in raw_calls
                if confidence + 1e-12 >= cutoff
            ]
            cutoff_matches = expected_cutoff_calls == cutoff_calls
            require(
                cutoff_matches,
                f"{source} cutoff column does not match raw confidence >= {cutoff:g} "
                f"at {path}:{line_number} ({feature_id})",
            )
            calls = {
                tag: _rank_call(tag, raw_calls, cutoff)
                for tag in expected_tags
            }
            assignments[feature_id] = SintaxAssignment(
                feature_id=feature_id,
                raw_assignment=raw_text,
                strand=strand,
                cutoff_assignment=cutoff_text,
                calls=calls,
                unexpected_tags=unexpected,
                cutoff_matches_raw_confidence=cutoff_matches,
            )
    require(bool(assignments), f"{source} SINTAX table contains no feature rows: {path}")
    return assignments


def read_pr2_dada2_crosswalk(path: Path) -> Pr2Crosswalk:
    """Build the official PR2 p-tag crosswalk from all 9-rank DADA2 headers."""

    require(path.is_file(), f"Missing PR2 DADA2 FASTA.gz: {path}")
    require(path.name.lower().endswith((".fasta.gz", ".fa.gz", ".fas.gz")), f"Expected gzipped FASTA: {path}")
    composite_to_ranks: dict[str, tuple[str, str]] = {}
    # Store fixed-size digests instead of all nine taxonomy strings so the
    # exhaustive schema audit remains memory-bounded for the 240k-record DB.
    unique_lineages: set[bytes] = set()
    record_count = 0
    current_has_sequence = False
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_count:
                    require(current_has_sequence, f"Empty PR2 DADA2 sequence before {path}:{line_number}")
                header = line[1:]
                fields = header.split(";")
                if fields and fields[-1] == "":
                    fields.pop()
                require(
                    len(fields) == 9,
                    f"PR2 DADA2 header at {path}:{line_number} has {len(fields)} ranks; expected 9",
                )
                require(all(fields), f"Empty PR2 DADA2 rank at {path}:{line_number}")
                lineage = tuple(fields)
                unique_lineages.add(hashlib.sha256("\x1f".join(lineage).encode("utf-8")).digest())
                division, subdivision = fields[2], fields[3]
                composite = f"{division}-{subdivision}"
                previous = composite_to_ranks.get(composite)
                require(
                    previous is None or previous == (division, subdivision),
                    f"PR2 composite crosswalk collision for {composite!r}: {previous} vs {(division, subdivision)}",
                )
                composite_to_ranks[composite] = (division, subdivision)
                record_count += 1
                current_has_sequence = False
            else:
                require(record_count > 0, f"Sequence before PR2 DADA2 header at {path}:{line_number}")
                require(
                    re.fullmatch(r"[ACGTUNRYKMSWBDHV.-]+", line.upper()) is not None,
                    f"Invalid nucleotide symbol in PR2 DADA2 FASTA at {path}:{line_number}",
                )
                current_has_sequence = True
    require(record_count > 0 and current_has_sequence, f"Empty or truncated PR2 DADA2 FASTA: {path}")
    require(bool(composite_to_ranks), f"No PR2 division/subdivision crosswalk entries in {path}")
    return Pr2Crosswalk(
        composite_to_ranks=composite_to_ranks,
        record_count=record_count,
        unique_lineage_count=len(unique_lineages),
        unique_composite_count=len(composite_to_ranks),
    )


def strip_pr2_compartment(value: str) -> tuple[str, str | None]:
    match = PR2_SUFFIX_RE.search(value)
    if match is None:
        return value, None
    return value[: match.start()], match.group("code").lower()


def accepted_value(assignment: SintaxAssignment, tag: str) -> str | None:
    return assignment.calls[tag].accepted_value


def _accepted_pr2_values(assignment: SintaxAssignment) -> tuple[str, ...]:
    values: list[str] = []
    for tag in PR2_TAGS:
        value = accepted_value(assignment, tag)
        if value is not None:
            values.append(value)
    return tuple(values)


def _accepted_silva_values(assignment: SintaxAssignment) -> tuple[str, ...]:
    return tuple(
        value
        for tag in SILVA_TAGS
        if (value := accepted_value(assignment, tag)) is not None
    )


def _canonical_or_unassigned(value: str | None, *, strip_pr2: bool = False) -> str:
    if value is None:
        return UNASSIGNED
    if strip_pr2:
        value, _ = strip_pr2_compartment(value)
    return value or UNASSIGNED


def canonicalize_pr2(
    assignment: SintaxAssignment,
    crosswalk: Pr2Crosswalk,
) -> tuple[tuple[str, ...], str, bool, tuple[str, ...]]:
    p_value = accepted_value(assignment, "p")
    crosswalk_missing = False
    division: str | None = None
    subdivision: str | None = None
    if p_value is not None:
        mapped = crosswalk.composite_to_ranks.get(p_value)
        if mapped is None:
            crosswalk_missing = True
        else:
            division, subdivision = mapped
    source_values = _accepted_pr2_values(assignment)
    compartments = tuple(
        sorted(
            {
                code
                for value in source_values
                for _, code in (strip_pr2_compartment(value),)
                if code is not None
            }
        )
    )
    canonical = (
        _canonical_or_unassigned(accepted_value(assignment, "k"), strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "d"), strip_pr2=True),
        _canonical_or_unassigned(division, strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "c"), strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "o"), strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "f"), strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "g"), strip_pr2=True),
        _canonical_or_unassigned(accepted_value(assignment, "s"), strip_pr2=True),
    )
    subdivision_text = _canonical_or_unassigned(subdivision, strip_pr2=True)
    return canonical, subdivision_text, crosswalk_missing, compartments


def canonicalize_silva(assignment: SintaxAssignment) -> tuple[str, ...]:
    return tuple(
        _canonical_or_unassigned(accepted_value(assignment, tag))
        for tag in SILVA_TAGS
    )


def is_placeholder(value: str) -> bool:
    normalized = value.strip().lower().replace(" ", "_")
    if not normalized or normalized == UNASSIGNED.lower():
        return True
    exact = {
        "unknown",
        "unidentified",
        "uncultured",
        "environmental_sample",
        "incertae_sedis",
        "incertae_sedis_x",
    }
    return (
        normalized in exact
        or normalized.endswith("_x")
        or normalized.endswith("_sp")
        or normalized.endswith("_sp.")
        or "environmental_sample" in normalized
        or normalized.startswith("uncultured_")
        or normalized.startswith("unidentified_")
    )


def comparison_key(value: str) -> str | None:
    if is_placeholder(value):
        return None
    stripped, _ = strip_pr2_compartment(value)
    stripped = TRAILING_ANNOTATION_RE.sub("", stripped)
    normalized = re.sub(r"[_\s-]+", "_", stripped.strip()).strip("_").lower()
    return normalized or None


def _agreement(
    pr2_canonical: tuple[str, ...],
    pr2_subdivision: str,
    silva_canonical: tuple[str, ...],
) -> tuple[str, str, str, str, str]:
    pr2_domain_key = comparison_key(pr2_canonical[0])
    silva_domain_key = comparison_key(silva_canonical[0])
    if pr2_domain_key is None and silva_domain_key is None:
        return "not_assessable", "both_unassigned", "", "", ""
    if pr2_domain_key is None:
        return "not_assessable", "pr2_unassigned", "", "", ""
    if silva_domain_key is None:
        return "not_assessable", "silva_unassigned", "", "", ""
    if pr2_domain_key != silva_domain_key:
        return "none", "domain_discordant", "", "", ""

    pr2_lineage = [
        ("Domain", pr2_canonical[0]),
        ("Supergroup", pr2_canonical[1]),
        ("Division", pr2_canonical[2]),
        ("Subdivision", pr2_subdivision),
        ("Class", pr2_canonical[3]),
        ("Order", pr2_canonical[4]),
        ("Family", pr2_canonical[5]),
        ("Genus", pr2_canonical[6]),
        ("Species", pr2_canonical[7]),
    ]
    silva_lineage = list(zip(CANONICAL_RANKS, silva_canonical))
    pr2_items = [
        (rank, value, comparison_key(value), index)
        for index, (rank, value) in enumerate(pr2_lineage)
        if comparison_key(value) is not None
    ]
    silva_items = [
        (rank, value, comparison_key(value), index)
        for index, (rank, value) in enumerate(silva_lineage)
        if comparison_key(value) is not None
    ]
    matches: list[tuple[int, int, str, str, str, str]] = []
    for pr2_rank, pr2_value, pr2_key, pr2_index in pr2_items:
        for silva_rank, silva_value, silva_key, silva_index in silva_items:
            if pr2_key == silva_key:
                matches.append(
                    (pr2_index, silva_index, pr2_rank, silva_rank, pr2_value, pr2_key or "")
                )
    if not matches:
        return "none", "discordant", "", "", ""
    deepest = max(matches, key=lambda row: (row[0] + row[1], min(row[0], row[1])))
    _, _, pr2_rank, silva_rank, agreement_taxon, _ = deepest
    if pr2_rank == "Domain" and silva_rank == "Domain":
        has_pr2_below = any(rank != "Domain" for rank, *_ in pr2_items)
        has_silva_below = any(rank != "Domain" for rank, *_ in silva_items)
        status = "discordant_below_domain" if has_pr2_below and has_silva_below else "domain_only"
        return "domain", status, agreement_taxon, pr2_rank, silva_rank
    same_rank = (
        pr2_rank == silva_rank
        or (pr2_rank == "Division" and silva_rank == "Division_or_Phylum")
    )
    level = pr2_rank.lower() if same_rank else "shared_taxon"
    status = "same_rank_concordant" if same_rank else "rank_shift_shared_taxon"
    return level, status, agreement_taxon, pr2_rank, silva_rank


FLAG_NAMES = (
    "flag_pr2_unassigned_root",
    "flag_pr2_ambiguous_assignment",
    "flag_pr2_crosswalk_missing",
    "flag_pr2_low_confidence",
    "flag_pr2_non_eukaryota",
    "flag_pr2_plastid",
    "flag_pr2_mitochondria",
    "flag_pr2_nucleomorph",
    "flag_pr2_exact_chloroplast_token",
    "flag_pr2_exact_mitochondria_token",
    "flag_pr2_metazoa",
    "flag_pr2_chloroplastida",
    "flag_pr2_placeholder",
    "flag_silva_unassigned_root",
    "flag_silva_ambiguous_assignment",
    "flag_silva_low_confidence",
    "flag_silva_metazoa",
    "flag_silva_chloroplastida",
    "flag_silva_placeholder",
    "flag_silva_non_eukaryota_not_assessable",
    "flag_crossdb_discordant",
    "flag_zero_total_count",
)


def _casefold_values(values: Iterable[str], *, strip_pr2: bool = False) -> set[str]:
    normalized: set[str] = set()
    for value in values:
        if strip_pr2:
            value, _ = strip_pr2_compartment(value)
        normalized.add(value.casefold())
    return normalized


def build_feature_result(
    feature_id: str,
    *,
    pr2: SintaxAssignment,
    silva: SintaxAssignment,
    crosswalk: Pr2Crosswalk,
    total_reads: int,
) -> FeatureResult:
    pr2_canonical, subdivision, crosswalk_missing, compartments = canonicalize_pr2(pr2, crosswalk)
    silva_canonical = canonicalize_silva(silva)
    pr2_values = _accepted_pr2_values(pr2)
    silva_values = _accepted_silva_values(silva)
    pr2_plain = _casefold_values(pr2_values, strip_pr2=True)
    silva_plain = _casefold_values(silva_values)
    root_call = pr2.calls["k"]
    root = accepted_value(pr2, "k")
    root_plain = strip_pr2_compartment(root)[0] if root is not None else None
    pr2_ambiguous = any(call.status.startswith("ambiguous") for call in pr2.calls.values())
    silva_ambiguous = any(call.status.startswith("ambiguous") for call in silva.calls.values())
    level, status, taxon, pr2_rank, silva_rank = _agreement(
        pr2_canonical,
        subdivision,
        silva_canonical,
    )
    flags = {
        "flag_pr2_unassigned_root": root_call.status in {"missing", "below_cutoff"},
        "flag_pr2_ambiguous_assignment": pr2_ambiguous,
        "flag_pr2_crosswalk_missing": crosswalk_missing,
        "flag_pr2_low_confidence": any(call.status.endswith("below_cutoff") or call.status == "below_cutoff" for call in pr2.calls.values()),
        "flag_pr2_non_eukaryota": root_plain is not None and root_plain.casefold() != "eukaryota",
        "flag_pr2_plastid": bool(set(compartments) & PLASTID_CODES) or "chloroplast" in pr2_plain,
        "flag_pr2_mitochondria": "mito" in compartments or "mitochondria" in pr2_plain,
        "flag_pr2_nucleomorph": "nucl" in compartments,
        "flag_pr2_exact_chloroplast_token": "chloroplast" in pr2_plain,
        "flag_pr2_exact_mitochondria_token": "mitochondria" in pr2_plain,
        "flag_pr2_metazoa": "metazoa" in pr2_plain,
        "flag_pr2_chloroplastida": "chloroplastida" in pr2_plain,
        "flag_pr2_placeholder": any(
            value != UNASSIGNED and is_placeholder(value)
            for value in pr2_canonical + (subdivision,)
        ),
        "flag_silva_unassigned_root": accepted_value(silva, "d") is None,
        "flag_silva_ambiguous_assignment": silva_ambiguous,
        "flag_silva_low_confidence": any(call.status.endswith("below_cutoff") or call.status == "below_cutoff" for call in silva.calls.values()),
        "flag_silva_metazoa": any(value in {"metazoa", "metazoa_(animalia)"} for value in silva_plain),
        "flag_silva_chloroplastida": "chloroplastida" in silva_plain,
        "flag_silva_placeholder": any(
            value != UNASSIGNED and is_placeholder(value)
            for value in silva_canonical
        ),
        # The audited SILVA v123 training FASTA is d:Eukaryota-only, so this source
        # cannot independently assess non-eukaryote status.
        "flag_silva_non_eukaryota_not_assessable": True,
        "flag_crossdb_discordant": status in {"domain_discordant", "discordant", "discordant_below_domain"},
        "flag_zero_total_count": total_reads == 0,
    }
    require(set(flags) == set(FLAG_NAMES), "Internal flag schema mismatch")

    exclusion_reasons: list[str] = []
    if root_call.status.startswith("ambiguous"):
        exclusion_reasons.append("pr2_ambiguous_root")
    elif flags["flag_pr2_unassigned_root"]:
        exclusion_reasons.append("pr2_unassigned_root")
    elif flags["flag_pr2_non_eukaryota"]:
        exclusion_reasons.append("pr2_non_eukaryota_root")
    if bool(set(compartments) & PLASTID_CODES):
        exclusion_reasons.append("pr2_plastid_compartment")
    if "mito" in compartments:
        exclusion_reasons.append("pr2_mitochondrial_compartment")
    if flags["flag_pr2_exact_chloroplast_token"]:
        exclusion_reasons.append("pr2_exact_chloroplast_token")
    if flags["flag_pr2_exact_mitochondria_token"]:
        exclusion_reasons.append("pr2_exact_mitochondria_token")
    if flags["flag_zero_total_count"]:
        exclusion_reasons.append("zero_total_count")
    analysis_keep = not exclusion_reasons
    return FeatureResult(
        feature_id=feature_id,
        pr2_canonical=pr2_canonical,
        silva_canonical=silva_canonical,
        pr2_subdivision=subdivision,
        pr2_compartments=compartments,
        flags=flags,
        analysis_keep=analysis_keep,
        exclusion_reason="retained_pr2_eukaryota" if analysis_keep else ";".join(exclusion_reasons),
        agreement_level=level,
        agreement_status=status,
        agreement_taxon=taxon,
        pr2_agreement_rank=pr2_rank,
        silva_agreement_rank=silva_rank,
    )


def assert_same_id_set(
    expected_ids: Sequence[str],
    observed_ids: Iterable[str],
    *,
    label: str,
) -> None:
    expected = set(expected_ids)
    observed = set(observed_ids)
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    require(
        not missing and not extra,
        f"Feature-ID set mismatch for {label}; missing={missing[:10]}, extra={extra[:10]}",
    )


def load_and_validate(
    *,
    count_path: Path,
    fasta_path: Path,
    pr2_sintax_path: Path,
    silva_sintax_path: Path,
    pr2_dada2_path: Path,
    metadata_path: Path,
    frozen_order_path: Path | None,
    cutoff: float,
    expected_samples: int,
    feature_prefix: str,
) -> PackageData:
    require(0.0 <= cutoff <= 1.0, f"SINTAX cutoff must be within [0,1], got {cutoff}")
    require(expected_samples > 0, f"Expected sample count must be positive, got {expected_samples}")
    require(bool(feature_prefix), "Feature prefix must not be empty")
    resolved_pr2 = pr2_sintax_path.resolve()
    resolved_silva = silva_sintax_path.resolve()
    require(resolved_pr2 != resolved_silva, "PR2 and SILVA SINTAX inputs must be different files")

    count = read_count_table(count_path)
    fasta = read_fasta(fasta_path)
    metadata_header, metadata_rows, metadata_order = read_metadata(metadata_path)
    frozen_order = read_frozen_order(frozen_order_path) if frozen_order_path else metadata_order
    pr2 = read_sintax(
        pr2_sintax_path,
        source="PR2 v5.1.1",
        expected_tags=PR2_TAGS,
        cutoff=cutoff,
    )
    silva = read_sintax(
        silva_sintax_path,
        source="SILVA 18S v123",
        expected_tags=SILVA_TAGS,
        cutoff=cutoff,
    )
    crosswalk = read_pr2_dada2_crosswalk(pr2_dada2_path)

    require(
        len(count.samples) == expected_samples,
        f"0.97 count table has {len(count.samples)} samples; expected {expected_samples}",
    )
    require(
        len(frozen_order) == expected_samples,
        f"Frozen order has {len(frozen_order)} samples; expected {expected_samples}",
    )
    require(
        count.samples == frozen_order,
        "0.97 count-table sample order does not exactly match the frozen order",
    )
    require(
        set(metadata_rows) == set(frozen_order),
        f"Metadata/frozen sample-set mismatch; missing={sorted(set(frozen_order) - set(metadata_rows))[:10]}, "
        f"extra={sorted(set(metadata_rows) - set(frozen_order))[:10]}",
    )
    ss_bs = [sample for sample in frozen_order if "SS_BS" in sample.upper()]
    require(not ss_bs, f"Frozen 18S design contains forbidden SS_BS sample(s): {ss_bs}")
    bad_prefix = [feature_id for feature_id in count.feature_ids if not feature_id.startswith(feature_prefix)]
    require(
        not bad_prefix,
        f"Feature IDs do not use required prefix {feature_prefix!r}: {bad_prefix[:10]}",
    )
    assert_same_id_set(count.feature_ids, fasta, label="representative FASTA")
    assert_same_id_set(count.feature_ids, pr2, label="PR2 SINTAX")
    assert_same_id_set(count.feature_ids, silva, label="SILVA v123 SINTAX")

    source_paths = {
        "count_0.97": count_path.resolve(),
        "representative_fasta": fasta_path.resolve(),
        "pr2_sintax": pr2_sintax_path.resolve(),
        "silva_v123_sintax": silva_sintax_path.resolve(),
        "pr2_dada2_9rank_crosswalk": pr2_dada2_path.resolve(),
        "metadata": metadata_path.resolve(),
    }
    if frozen_order_path is not None:
        source_paths["frozen_order"] = frozen_order_path.resolve()
    source_hashes = {label: sha256_file(path) for label, path in source_paths.items()}
    require(
        source_hashes["pr2_sintax"] != source_hashes["silva_v123_sintax"],
        "PR2 and SILVA SINTAX inputs have the same SHA-256; refusing possible comparator mix-up",
    )

    results: dict[str, FeatureResult] = {}
    crosswalk_missing_ids: list[str] = []
    for feature_id in count.feature_ids:
        total_reads = sum(count.rows[feature_id])
        result = build_feature_result(
            feature_id,
            pr2=pr2[feature_id],
            silva=silva[feature_id],
            crosswalk=crosswalk,
            total_reads=total_reads,
        )
        results[feature_id] = result
        if result.flags["flag_pr2_crosswalk_missing"]:
            crosswalk_missing_ids.append(feature_id)
    require(
        not crosswalk_missing_ids,
        "Accepted PR2 p-tag values are absent from the supplied official DADA2 9-rank crosswalk; "
        f"possible PR2 version mismatch for {crosswalk_missing_ids[:10]}",
    )
    analysis_ids = tuple(
        feature_id
        for feature_id in count.feature_ids
        if results[feature_id].analysis_keep
    )
    require(bool(analysis_ids), "PR2 conservative target filter retained no positive-count ASVs")
    for sample_index, sample in enumerate(count.samples):
        all_depth = sum(row[sample_index] for row in count.rows.values())
        analysis_depth = sum(count.rows[feature_id][sample_index] for feature_id in analysis_ids)
        require(
            analysis_depth <= all_depth,
            f"Internal sample-depth additivity failure for {sample}",
        )

    return PackageData(
        source_paths=source_paths,
        source_hashes=source_hashes,
        count=count,
        fasta=fasta,
        metadata_header=metadata_header,
        metadata_rows=metadata_rows,
        frozen_order=frozen_order,
        pr2=pr2,
        silva=silva,
        crosswalk=crosswalk,
        results=results,
        analysis_ids=analysis_ids,
        cutoff=cutoff,
        feature_prefix=feature_prefix,
    )


def write_tsv(
    path: Path,
    header: Sequence[object] | None,
    rows: Iterable[Sequence[object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header is not None:
            writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def write_fasta(path: Path, feature_ids: Sequence[str], sequences: Mapping[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for feature_id in feature_ids:
            handle.write(f">{feature_id}\n")
            for line in textwrap.wrap(sequences[feature_id], width=80):
                handle.write(f"{line}\n")


def canonical_header() -> tuple[str, ...]:
    return ("ASV_ID", *CANONICAL_RANKS)


def canonical_rows(
    data: PackageData,
    feature_ids: Sequence[str],
) -> Iterator[tuple[object, ...]]:
    for feature_id in feature_ids:
        yield (feature_id, *data.results[feature_id].pr2_canonical)


def _source_rank_columns(prefix: str, names: Mapping[str, str], tags: Sequence[str]) -> list[str]:
    columns: list[str] = []
    for tag in tags:
        base = f"{prefix}_{names[tag]}"
        columns.extend(
            (
                f"{base}_raw_values",
                f"{base}_confidences",
                f"{base}_accepted_values",
                f"{base}_status",
            )
        )
    return columns


def parallel_header() -> tuple[str, ...]:
    columns = [
        "ASV_ID",
        "pr2_source_database",
        "pr2_raw_assignment",
        "pr2_strand",
        "pr2_cutoff_assignment",
        "pr2_cutoff",
    ]
    columns.extend(_source_rank_columns("pr2", PR2_SOURCE_NAMES, PR2_TAGS))
    columns.extend(f"pr2_canonical_{rank.lower()}" for rank in CANONICAL_RANKS)
    columns.extend(("pr2_canonical_subdivision", "pr2_compartment_codes"))
    columns.extend(
        (
            "silva_source_database",
            "silva_raw_assignment",
            "silva_strand",
            "silva_cutoff_assignment",
            "silva_cutoff",
        )
    )
    columns.extend(_source_rank_columns("silva", SILVA_SOURCE_NAMES, SILVA_TAGS))
    columns.extend(f"silva_canonical_{rank.lower()}" for rank in CANONICAL_RANKS)
    columns.extend(
        (
            "agreement_level",
            "agreement_status",
            "agreement_taxon",
            "pr2_agreement_source_rank",
            "silva_agreement_source_rank",
            "agreement_method",
        )
    )
    return tuple(columns)


def _rank_call_cells(call: RankCall) -> tuple[str, str, str, str]:
    return (
        "|".join(call.raw_values),
        "|".join(format_confidence(value) for value in call.confidences),
        "|".join(call.accepted_values),
        call.status,
    )


def parallel_rows(data: PackageData) -> Iterator[tuple[object, ...]]:
    for feature_id in data.count.feature_ids:
        pr2 = data.pr2[feature_id]
        silva = data.silva[feature_id]
        result = data.results[feature_id]
        cells: list[object] = [
            feature_id,
            "PR2_v5.1.1_SSU_UTAX",
            pr2.raw_assignment,
            pr2.strand,
            pr2.cutoff_assignment,
            format_confidence(data.cutoff),
        ]
        for tag in PR2_TAGS:
            cells.extend(_rank_call_cells(pr2.calls[tag]))
        cells.extend(result.pr2_canonical)
        cells.extend((result.pr2_subdivision, ";".join(result.pr2_compartments)))
        cells.extend(
            (
                "SILVA_18S_v123_SINTAX_comparator",
                silva.raw_assignment,
                silva.strand,
                silva.cutoff_assignment,
                format_confidence(data.cutoff),
            )
        )
        for tag in SILVA_TAGS:
            cells.extend(_rank_call_cells(silva.calls[tag]))
        cells.extend(result.silva_canonical)
        cells.extend(
            (
                result.agreement_level,
                result.agreement_status,
                result.agreement_taxon,
                result.pr2_agreement_rank,
                result.silva_agreement_rank,
                "exact normalized taxon labels only; controlled PR2 compartment suffix and terminal parenthetical annotation removed; no fuzzy matching or cross-source backfill",
            )
        )
        require(len(cells) == len(parallel_header()), f"Internal parallel-row width error for {feature_id}")
        yield tuple(cells)


def ledger_header() -> tuple[str, ...]:
    return (
        "ASV_ID",
        "total_reads_0.97",
        "positive_sample_count_0.97",
        "analysis_keep",
        "exclusion_reason",
        *FLAG_NAMES,
        "pr2_accepted_root_raw",
        "pr2_canonical_domain",
        "pr2_compartment_codes",
        "silva_accepted_root_raw",
        "agreement_level",
        "agreement_status",
    )


def ledger_rows(data: PackageData) -> Iterator[tuple[object, ...]]:
    for feature_id in data.count.feature_ids:
        counts = data.count.rows[feature_id]
        result = data.results[feature_id]
        pr2_root = accepted_value(data.pr2[feature_id], "k") or ""
        silva_root = accepted_value(data.silva[feature_id], "d") or ""
        yield (
            feature_id,
            sum(counts),
            sum(value > 0 for value in counts),
            bool_text(result.analysis_keep),
            result.exclusion_reason,
            *(bool_text(result.flags[name]) for name in FLAG_NAMES),
            pr2_root,
            result.pr2_canonical[0],
            ";".join(result.pr2_compartments),
            silva_root,
            result.agreement_level,
            result.agreement_status,
        )


def sample_depth_header(data: PackageData) -> tuple[str, ...]:
    metadata_columns = tuple(f"metadata_{name.lstrip('#')}" for name in data.metadata_header[1:])
    return (
        "SampleID",
        "all_feature_reads_0.97",
        "pr2_analysis_reads",
        "excluded_reads",
        "pr2_analysis_retained_fraction",
        "zero_depth_after_pr2_filter",
        *metadata_columns,
    )


def sample_depth_values(data: PackageData) -> tuple[tuple[int, ...], tuple[int, ...]]:
    all_depths: list[int] = []
    analysis_depths: list[int] = []
    for index, _sample in enumerate(data.frozen_order):
        all_depths.append(sum(row[index] for row in data.count.rows.values()))
        analysis_depths.append(
            sum(data.count.rows[feature_id][index] for feature_id in data.analysis_ids)
        )
    return tuple(all_depths), tuple(analysis_depths)


def sample_depth_rows(data: PackageData) -> Iterator[tuple[object, ...]]:
    all_depths, analysis_depths = sample_depth_values(data)
    for sample, all_depth, analysis_depth in zip(data.frozen_order, all_depths, analysis_depths):
        excluded = all_depth - analysis_depth
        require(excluded >= 0 and analysis_depth + excluded == all_depth, f"Depth additivity failure for {sample}")
        retained_fraction = analysis_depth / all_depth if all_depth else 0.0
        yield (
            sample,
            all_depth,
            analysis_depth,
            excluded,
            f"{retained_fraction:.8f}",
            bool_text(analysis_depth == 0),
            *data.metadata_rows[sample][1:],
        )


def _nearest_rank(values: Sequence[int], fraction: float) -> int:
    require(bool(values), "Cannot calculate a quantile from an empty depth vector")
    ordered = sorted(values)
    index = max(0, math.ceil(fraction * len(ordered)) - 1)
    return ordered[index]


def rarefaction_candidate_header() -> tuple[str, ...]:
    return (
        "candidate_depth",
        "candidate_sources",
        "samples_retained",
        "samples_dropped",
        "sample_retention_fraction",
        "projected_rarefied_reads",
        "observed_analysis_reads_in_retained_samples",
        "projected_fraction_of_retained_observed_reads",
        "all_frozen_samples_retained",
        "dropped_sample_ids",
    )


def rarefaction_candidate_rows(data: PackageData) -> Iterator[tuple[object, ...]]:
    _, analysis_depths = sample_depth_values(data)
    positive_depths = [value for value in analysis_depths if value > 0]
    require(bool(positive_depths), "No positive PR2-analysis sample depth is available")
    candidates: dict[int, set[str]] = {}

    def add(depth: int, label: str) -> None:
        if depth > 0:
            candidates.setdefault(depth, set()).add(label)

    for sample, depth in zip(data.frozen_order, analysis_depths):
        add(depth, f"observed_depth:{sample}")
    for fraction, label in (
        (0.05, "nearest_rank_q05"),
        (0.10, "nearest_rank_q10"),
        (0.25, "nearest_rank_q25"),
        (0.50, "nearest_rank_median"),
        (0.75, "nearest_rank_q75"),
    ):
        add(_nearest_rank(positive_depths, fraction), label)
    add(min(positive_depths), "minimum_positive_depth")
    add(max(positive_depths), "maximum_depth")
    maximum = max(positive_depths)
    for standard in (1000, 2500, 5000, 7500, 10000, 15000, 20000, 25000, 30000, 50000, 75000, 100000):
        if standard <= maximum:
            add(standard, "standard_decision_depth")

    sample_count = len(data.frozen_order)
    for depth in sorted(candidates):
        retained = [
            sample
            for sample, observed in zip(data.frozen_order, analysis_depths)
            if observed >= depth
        ]
        dropped = [
            sample
            for sample, observed in zip(data.frozen_order, analysis_depths)
            if observed < depth
        ]
        retained_observed_reads = sum(
            observed for observed in analysis_depths if observed >= depth
        )
        projected = depth * len(retained)
        projected_fraction = projected / retained_observed_reads if retained_observed_reads else 0.0
        yield (
            depth,
            ";".join(sorted(candidates[depth])),
            len(retained),
            len(dropped),
            f"{len(retained) / sample_count:.8f}",
            projected,
            retained_observed_reads,
            f"{projected_fraction:.8f}",
            bool_text(len(retained) == sample_count),
            ";".join(dropped),
        )


def summary_dict(data: PackageData) -> dict[str, object]:
    all_depths, analysis_depths = sample_depth_values(data)
    flag_counts = {
        name: sum(data.results[feature_id].flags[name] for feature_id in data.count.feature_ids)
        for name in FLAG_NAMES
    }
    exclusion_counts: dict[str, int] = {}
    agreement_counts: dict[str, int] = {}
    for feature_id in data.count.feature_ids:
        result = data.results[feature_id]
        exclusion_counts[result.exclusion_reason] = exclusion_counts.get(result.exclusion_reason, 0) + 1
        agreement_counts[result.agreement_status] = agreement_counts.get(result.agreement_status, 0) + 1
    total_all = sum(all_depths)
    total_analysis = sum(analysis_depths)
    return {
        "package_version": VERSION,
        "marker": "18S rRNA gene V4",
        "primary_taxonomy": "PR2 v5.1.1 SSU UTAX SINTAX",
        "comparator_taxonomy": "legacy SILVA 18S v123 SINTAX",
        "sintax_cutoff": data.cutoff,
        "sample_contract": {
            "count": len(data.frozen_order),
            "count_header_exactly_matches_frozen_order": True,
            "metadata_joined_by_sample_id": True,
            "contains_SS_BS": False,
        },
        "feature_contract": {
            "all_features": len(data.count.feature_ids),
            "pr2_analysis_features": len(data.analysis_ids),
            "excluded_features": len(data.count.feature_ids) - len(data.analysis_ids),
            "count_fasta_pr2_sintax_silva_sintax_id_sets_equal": True,
            "feature_prefix": data.feature_prefix,
        },
        "read_counts": {
            "all_features_0.97": total_all,
            "pr2_analysis": total_analysis,
            "excluded": total_all - total_analysis,
            "retained_fraction": total_analysis / total_all if total_all else 0.0,
        },
        "sample_depths": {
            "all_minimum": min(all_depths),
            "all_maximum": max(all_depths),
            "analysis_minimum": min(analysis_depths),
            "analysis_maximum": max(analysis_depths),
            "analysis_zero_depth_samples": sum(value == 0 for value in analysis_depths),
        },
        "pr2_dada2_crosswalk": {
            "records_scanned": data.crosswalk.record_count,
            "unique_9rank_lineages": data.crosswalk.unique_lineage_count,
            "unique_division_subdivision_composites": data.crosswalk.unique_composite_count,
            "accepted_pr2_p_values_missing_from_crosswalk": 0,
        },
        "flag_counts": dict(sorted(flag_counts.items())),
        "exclusion_reason_counts": dict(sorted(exclusion_counts.items())),
        "agreement_status_counts": dict(sorted(agreement_counts.items())),
        "filter_contract": {
            "auto_excluded": [
                "PR2 root unassigned at cutoff",
                "PR2 root ambiguous at cutoff",
                "PR2 root other than Eukaryota after controlled compartment-suffix normalization",
                "exactly parsed PR2 plastid/apicoplast/chromatophore suffix (:plas, :apic, :chro)",
                "exactly parsed PR2 mitochondrial suffix (:mito)",
                "exact accepted PR2 token Chloroplast",
                "exact accepted PR2 token Mitochondria",
                "zero total count across frozen samples",
            ],
            "flag_only": [
                "Metazoa",
                "Chloroplastida",
                "PR2 nucleomorph compartment suffix",
            ],
        },
        "normalization_policy": {
            "pr2_rank_semantics": {
                "k": "domain",
                "d": "supergroup (represented in canonical Kingdom column with source semantics preserved in wide table)",
                "p": "division-subdivision composite resolved against official PR2 v5.1.1 DADA2 9-rank headers",
            },
            "silva_rank_semantics": "d/domain, k/kingdom, p/phylum, c/class, o/order, f/family, g/genus, s/species",
            "cross_source_backfill": False,
            "agreement": "exact normalized labels only; no fuzzy or synonym inference",
        },
        "rarefaction_executed": False,
        "plots_generated": False,
    }


def readme_text(data: PackageData) -> str:
    summary = summary_dict(data)
    return f"""# 18S parallel taxonomy post-processing package ({VERSION})

This helper package was built from the completed 0.97 read-recruitment count
table and representative ASVs. It does not run SINTAX, rarefaction, diversity
analysis, or plotting.

## Taxonomy contract

- PR2 v5.1.1 is the primary taxonomy. SILVA 18S v123 is an independent legacy
  comparator and never fills a missing PR2 rank.
- Every raw assignment, rank confidence, cutoff-accepted value, and parsing
  status is retained in `{OUTPUT_FILES['parallel']}`.
- PR2 rank letters are interpreted from the audited database schema: `k` is
  domain, `d` is supergroup, and `p` is the composite
  `division-subdivision`. The `p` value is resolved through all
  {data.crosswalk.record_count:,} headers in the supplied official PR2 DADA2
  9-rank FASTA; it is never blindly split on punctuation.
- The primary canonical table uses Domain, Kingdom, Division_or_Phylum, Class,
  Order, Family, Genus, and Species. Its Kingdom column contains the PR2
  supergroup by explicit mapping, while the wide table preserves that source
  meaning and the PR2 subdivision separately.
- Duplicate distinct values for a source rank are `ambiguous`; no later value
  silently overwrites an earlier one. Agreement uses exact normalized labels
  only. Controlled PR2 compartment suffixes and terminal parenthetical common
  names are removed for comparison, with no fuzzy matching or synonym inference.

## Analysis filter and flags

At SINTAX cutoff {data.cutoff:g}, {len(data.analysis_ids):,} of
{len(data.count.feature_ids):,} ASVs are retained in the PR2 analysis tier.
The conservative core contract excludes an unassigned/ambiguous/non-Eukaryota
PR2 root, exactly parsed PR2-native compartment suffixes `:plas`, `:apic`,
`:chro`, or `:mito`, exact accepted tokens `Chloroplast` or `Mitochondria`, and
zero-total features. `Chloroplastida` is a nuclear plant lineage and remains,
as does Metazoa; both are flags only. The PR2 `:nucl` nucleomorph suffix is
also a flag only. All controlled suffixes remain visible in the ledger.

The SILVA v123 comparator database is Eukaryota-rooted and therefore cannot
independently assess non-eukaryote status. This limitation is explicit in the
ASV ledger.

## Depth-decision files

`{OUTPUT_FILES['sample_depths']}` reports the all-feature and PR2-analysis
library size for every frozen sample with an additive excluded-read column.
`{OUTPUT_FILES['rarefaction_candidates']}` is a deterministic decision table
covering observed depths, nearest-rank quantiles, and common round thresholds.
It estimates sample retention only; no count table is subsampled.

## Recheck

Run the helper again with the same explicit inputs and `--check`. Source and
output SHA-256 values, ID sets, sample order, all table cells, FASTA sequences,
flag ledgers, and depth calculations are recomputed.

Summary: {summary['feature_contract']['all_features']:,} all ASVs;
{summary['feature_contract']['pr2_analysis_features']:,} analysis ASVs;
{summary['read_counts']['pr2_analysis']:,} retained reads across
{summary['sample_contract']['count']} samples.
"""


def write_package(build_root: Path, data: PackageData) -> None:
    write_tsv(
        build_root / OUTPUT_FILES["taxonomy_all"],
        canonical_header(),
        canonical_rows(data, data.count.feature_ids),
    )
    write_tsv(
        build_root / OUTPUT_FILES["taxonomy_analysis"],
        canonical_header(),
        canonical_rows(data, data.analysis_ids),
    )
    write_tsv(
        build_root / OUTPUT_FILES["parallel"],
        parallel_header(),
        parallel_rows(data),
    )
    write_tsv(
        build_root / OUTPUT_FILES["counts_analysis"],
        (data.count.id_header, *data.frozen_order),
        (
            (feature_id, *data.count.rows[feature_id])
            for feature_id in data.analysis_ids
        ),
    )
    write_fasta(
        build_root / OUTPUT_FILES["fasta_analysis"],
        data.analysis_ids,
        data.fasta,
    )
    write_tsv(
        build_root / OUTPUT_FILES["ledger"],
        ledger_header(),
        ledger_rows(data),
    )
    write_tsv(
        build_root / OUTPUT_FILES["sample_depths"],
        sample_depth_header(data),
        sample_depth_rows(data),
    )
    write_tsv(
        build_root / OUTPUT_FILES["rarefaction_candidates"],
        rarefaction_candidate_header(),
        rarefaction_candidate_rows(data),
    )

    summary_path = build_root / OUTPUT_FILES["summary"]
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(summary_dict(data), handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    readme_path = build_root / OUTPUT_FILES["readme"]
    readme_path.parent.mkdir(parents=True, exist_ok=True)
    readme_path.write_text(readme_text(data), encoding="utf-8", newline="\n")

    manifest_rows: list[tuple[object, ...]] = []
    for label, path in sorted(data.source_paths.items()):
        manifest_rows.append(
            ("source", label, str(path), path.stat().st_size, data.source_hashes[label])
        )
    for label, relative_path in sorted(OUTPUT_FILES.items()):
        if label == "manifest":
            continue
        path = build_root / relative_path
        manifest_rows.append(
            ("output", label, relative_path, path.stat().st_size, sha256_file(path))
        )
    write_tsv(
        build_root / OUTPUT_FILES["manifest"],
        ("category", "label", "path", "bytes", "sha256"),
        manifest_rows,
    )


def assert_tsv_exact(
    path: Path,
    expected_header: Sequence[object],
    expected_rows: Iterable[Sequence[object]],
    *,
    label: str,
) -> None:
    require(path.is_file(), f"Missing generated {label}: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            observed_header = tuple(next(reader))
        except StopIteration as exc:
            raise ValidationError(f"Empty generated {label}: {path}") from exc
        expected_header_text = tuple(str(cell) for cell in expected_header)
        require(
            observed_header == expected_header_text,
            f"Generated {label} header mismatch: {observed_header} != {expected_header_text}",
        )
        observed_rows = (
            tuple(row)
            for row in reader
            if row and not all(cell == "" for cell in row)
        )
        expected_text_rows = (tuple(str(cell) for cell in row) for row in expected_rows)
        sentinel = object()
        for row_number, (observed, expected) in enumerate(
            itertools.zip_longest(observed_rows, expected_text_rows, fillvalue=sentinel),
            start=2,
        ):
            require(observed is not sentinel, f"Generated {label} ended early before row {row_number}")
            require(expected is not sentinel, f"Generated {label} has an extra row at {row_number}")
            require(
                observed == expected,
                f"Generated {label} differs at row {row_number}: observed={observed!r}, expected={expected!r}",
            )


def validate_manifest(output_root: Path, data: PackageData) -> None:
    manifest_path = output_root / OUTPUT_FILES["manifest"]
    header, rows = read_tsv(manifest_path)
    require(
        header == ("category", "label", "path", "bytes", "sha256"),
        f"Unexpected helper manifest header: {header}",
    )
    expected_keys = {
        *(('source', label) for label in data.source_paths),
        *(('output', label) for label in OUTPUT_FILES if label != 'manifest'),
    }
    seen: set[tuple[str, str]] = set()
    for row in rows:
        category, label, recorded_path, bytes_text, recorded_hash = row
        key = (category, label)
        require(key not in seen, f"Duplicate helper manifest row: {key}")
        seen.add(key)
        if category == "source":
            require(label in data.source_paths, f"Unknown source manifest label: {label}")
            path = data.source_paths[label]
            require(recorded_path == str(path), f"Source manifest path mismatch for {label}")
            require(recorded_hash == data.source_hashes[label], f"Source SHA-256 drift for {label}")
        elif category == "output":
            require(label in OUTPUT_FILES and label != "manifest", f"Unknown output manifest label: {label}")
            require(recorded_path == OUTPUT_FILES[label], f"Output manifest path mismatch for {label}")
            path = output_root / recorded_path
        else:
            raise ValidationError(f"Unknown manifest category: {category}")
        require(path.is_file(), f"Manifest target missing: {path}")
        require(path.stat().st_size == int(bytes_text), f"Manifest byte-size drift for {path}")
        require(sha256_file(path) == recorded_hash, f"Manifest SHA-256 drift for {path}")
    require(
        seen == expected_keys,
        f"Helper manifest row-set mismatch; missing={sorted(expected_keys - seen)}, extra={sorted(seen - expected_keys)}",
    )


def check_package(output_root: Path, data: PackageData) -> dict[str, object]:
    require(output_root.is_dir(), f"Core output package does not exist: {output_root}")
    for relative_path in OUTPUT_FILES.values():
        require((output_root / relative_path).is_file(), f"Missing helper output: {output_root / relative_path}")
    validate_manifest(output_root, data)

    assert_tsv_exact(
        output_root / OUTPUT_FILES["taxonomy_all"],
        canonical_header(),
        canonical_rows(data, data.count.feature_ids),
        label="PR2 primary all-feature taxonomy",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["taxonomy_analysis"],
        canonical_header(),
        canonical_rows(data, data.analysis_ids),
        label="PR2 primary analysis taxonomy",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["parallel"],
        parallel_header(),
        parallel_rows(data),
        label="PR2/SILVA parallel taxonomy",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["counts_analysis"],
        (data.count.id_header, *data.frozen_order),
        ((feature_id, *data.count.rows[feature_id]) for feature_id in data.analysis_ids),
        label="PR2 analysis count table",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["ledger"],
        ledger_header(),
        ledger_rows(data),
        label="ASV flag ledger",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["sample_depths"],
        sample_depth_header(data),
        sample_depth_rows(data),
        label="sample-depth table",
    )
    assert_tsv_exact(
        output_root / OUTPUT_FILES["rarefaction_candidates"],
        rarefaction_candidate_header(),
        rarefaction_candidate_rows(data),
        label="rarefaction candidate table",
    )

    observed_fasta = read_fasta(output_root / OUTPUT_FILES["fasta_analysis"])
    require(
        tuple(observed_fasta) == data.analysis_ids,
        "Generated PR2 analysis FASTA ID/order differs from the filtered count table",
    )
    for feature_id in data.analysis_ids:
        require(
            observed_fasta[feature_id] == data.fasta[feature_id],
            f"Generated PR2 analysis FASTA sequence differs for {feature_id}",
        )

    with (output_root / OUTPUT_FILES["summary"]).open("r", encoding="utf-8") as handle:
        observed_summary = json.load(handle)
    expected_summary = summary_dict(data)
    require(observed_summary == expected_summary, "Generated helper summary differs from recomputed inputs")
    observed_readme = (output_root / OUTPUT_FILES["readme"]).read_text(encoding="utf-8")
    require(observed_readme == readme_text(data), "Generated helper README differs from recomputed inputs")
    return expected_summary


def resolve_path(value: str, *, project_root: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = project_root / path
    return path.resolve()


def resolve_project_root(script_path: Path, supplied: str | None) -> Path:
    root = Path(supplied).expanduser().resolve() if supplied else script_path.resolve().parents[3]
    require((root / "revision" / "re_analysis").is_dir(), f"Invalid 2024-amplicon project root: {root}")
    return root


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build or verify PR2-primary canonical taxonomy, a source-separated PR2/SILVA-v123 "
            "parallel table, conservative PR2 analysis artifacts, and depth-decision tables."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--project-root", help="2024-amplicon project root; inferred from this script when omitted")
    parser.add_argument(
        "--count-table",
        default="revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_asv_counts_97_all_v1.tsv",
        help="Completed upstream 0.97 ASV count table",
    )
    parser.add_argument(
        "--rep-fasta",
        default="revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_rep_seqs_all_v1.fasta",
        help="Representative ASV FASTA with the exact count-table feature set",
    )
    parser.add_argument(
        "--pr2-sintax",
        default="revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_sintax_pr2_v5.1.1_all_v1.tsv",
        help="Raw 4-column PR2 v5.1.1 VSEARCH SINTAX output",
    )
    parser.add_argument(
        "--silva-sintax",
        default="revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_sintax_silva_v123_all_v1.tsv",
        help="Raw 4-column legacy SILVA 18S v123 VSEARCH SINTAX output",
    )
    parser.add_argument(
        "--pr2-dada2-fasta-gz",
        default="revision/re_analysis/new_db/pr2_version_5.1.1_SSU_dada2.fasta.gz",
        help="Official PR2 v5.1.1 DADA2 FASTA.gz with 9-rank taxonomy headers",
    )
    parser.add_argument(
        "--metadata",
        default="revision/re_analysis/02_18s_unoise_pr2_v1/provenance/metadata_18s_75_v1.tsv",
        help="75-sample metadata TSV keyed by SampleID",
    )
    parser.add_argument(
        "--frozen-order",
        default="revision/re_analysis/01_input_16s_75_v1/16s_asv_counts_all_75_v1.tsv",
        help="Frozen sample order as a count header, metadata TSV, or one SampleID per line",
    )
    parser.add_argument(
        "--output-dir",
        default="revision/re_analysis/02_18s_unoise_pr2_v1",
        help="Existing core package root; helper files are added under results/qc/provenance",
    )
    parser.add_argument("--cutoff", type=float, default=0.1, help="SINTAX confidence cutoff used by both input tables")
    parser.add_argument("--expected-samples", type=int, default=75, help="Exact frozen sample count")
    parser.add_argument("--feature-prefix", default="18S_ASV_", help="Required prefix for every feature ID")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true", help="Create helper outputs after full validation; refuses every existing target")
    mode.add_argument("--dry-run", action="store_true", help="Read and validate all inputs, then print the plan without writing")
    mode.add_argument("--check", action="store_true", help="Read-only validation of existing helper outputs against current inputs")
    return parser.parse_args(argv)


def _load_from_args(args: argparse.Namespace, project_root: Path) -> PackageData:
    return load_and_validate(
        count_path=resolve_path(args.count_table, project_root=project_root),
        fasta_path=resolve_path(args.rep_fasta, project_root=project_root),
        pr2_sintax_path=resolve_path(args.pr2_sintax, project_root=project_root),
        silva_sintax_path=resolve_path(args.silva_sintax, project_root=project_root),
        pr2_dada2_path=resolve_path(args.pr2_dada2_fasta_gz, project_root=project_root),
        metadata_path=resolve_path(args.metadata, project_root=project_root),
        frozen_order_path=resolve_path(args.frozen_order, project_root=project_root) if args.frozen_order else None,
        cutoff=args.cutoff,
        expected_samples=args.expected_samples,
        feature_prefix=args.feature_prefix,
    )


def _refuse_existing_targets(output_root: Path) -> None:
    existing = [relative for relative in OUTPUT_FILES.values() if (output_root / relative).exists()]
    require(not existing, f"Refusing to overwrite existing helper output(s): {existing}")


def install_staged_outputs(stage_root: Path, output_root: Path) -> None:
    _refuse_existing_targets(output_root)
    moved: list[Path] = []
    try:
        for relative_path in OUTPUT_FILES.values():
            source = stage_root / relative_path
            target = output_root / relative_path
            require(source.is_file(), f"Staged helper output is missing: {source}")
            require(not target.exists(), f"Refusing late overwrite of helper output: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            source.replace(target)
            moved.append(target)
    except Exception:
        for target in reversed(moved):
            if target.is_file():
                target.unlink()
        raise


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    project_root = resolve_project_root(Path(__file__), args.project_root)
    output_root = resolve_path(args.output_dir, project_root=project_root)
    data = _load_from_args(args, project_root)

    if args.dry_run:
        payload = summary_dict(data)
        payload["planned_helper_outputs"] = list(OUTPUT_FILES.values())
        payload["output_root"] = str(output_root)
        print(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True))
        print("DRY RUN PASSED: all inputs validated; no files written")
        return 0

    if args.check:
        summary = check_package(output_root, data)
        print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
        print(f"CHECK PASSED: {output_root}")
        return 0

    require(args.run, "Internal mode guard: helper outputs may be created only with explicit --run")
    require(output_root.is_dir(), f"Core output package does not exist: {output_root}")
    for required_subdir in ("results", "qc", "provenance"):
        require(
            (output_root / required_subdir).is_dir(),
            f"Core package is missing required directory {required_subdir!r}: {output_root}",
        )
    _refuse_existing_targets(output_root)
    stage_root = Path(tempfile.mkdtemp(prefix=".parallel_taxonomy_v1_build_", dir=output_root))
    try:
        write_package(stage_root, data)
        check_package(stage_root, data)
        install_staged_outputs(stage_root, output_root)
    finally:
        if (
            stage_root.exists()
            and stage_root.parent.resolve() == output_root.resolve()
            and stage_root.name.startswith(".parallel_taxonomy_v1_build_")
        ):
            shutil.rmtree(stage_root)
    print(f"HELPER OUTPUTS CREATED AND VERIFIED: {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValidationError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
