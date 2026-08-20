#!/usr/bin/env python3
"""Run the versioned 75-sample 18S V4 UNOISE3 reanalysis.

The confirmed original analysis notebook is reconstructed as checkpointed,
Windows-native commands.  Two confirmed defects are corrected: the input is
strictly pinned to raw_data/seq_18s and merged reads are stripped by the
observed 20-nt forward and 18-nt reverse primers.  Rarefaction is never run by
this program.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence


VERSION = "1.1.0"
RUN_SLUG = "02_18s_unoise_pr2_v1"
EXPECTED_SAMPLES = 75
EXPECTED_FASTQ_FILES = 150
FORWARD_PRIMER = "CCAGCASCYGCGGTAATTCC"
REVERSE_PRIMER = "ACTTTCGTTCTTGATYRA"
FORWARD_REGEX = re.compile(r"^CCAGCA[GC]C[CT]GCGGTAATTCC", re.IGNORECASE)
REVERSE_REGEX = re.compile(r"^ACTTTCGTTCTTGAT[CT][AG]A", re.IGNORECASE)
WRONG_16S_REGEX = re.compile(r"^CCTACGGG", re.IGNORECASE)
STAGE_NAMES = (
    "01_seqkit",
    "02_merge",
    "03_filter",
    "04_derep",
    "05_unoise",
    "06_chimera",
    "07_map97",
    "08_map100",
    "09_tax_pr2",
    "10_tax_silva",
    "11_core_qc",
)


class ContractError(RuntimeError):
    """Raised when an input, command, or output violates the run contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_fingerprint(path: Path, *, include_hash: bool = True) -> dict[str, object]:
    require(path.is_file(), f"Missing file: {path}")
    stat = path.stat()
    record: dict[str, object] = {
        "path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_hash:
        record["sha256"] = sha256_file(path)
    return record


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def write_tsv(path: Path, header: Sequence[object], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def unique(values: Sequence[str], label: str) -> None:
    counts = Counter(values)
    duplicates = sorted(value for value, count in counts.items() if count > 1)
    require(not duplicates, f"Duplicate {label}: {duplicates[:10]}")


def project_root_from_script() -> Path:
    return Path("/user/link/path")  # insert your path


@dataclass(frozen=True)
class SampleInput:
    sample_id: str
    r1: Path
    r2: Path
    metadata: tuple[str, ...]


@dataclass(frozen=True)
class ToolPaths:
    seqkit: Path
    vsearch: Path
    usearch12: Path
    python: Path


@dataclass(frozen=True)
class RunPaths:
    project_root: Path
    output_root: Path
    frozen_count: Path
    metadata: Path
    fastq_dir: Path
    pr2_utax: Path
    pr2_dada2: Path
    silva_v123: Path
    original_script: Path
    helper_script: Path


def resolve_executable(name: str) -> Path:
    candidate = shutil.which(name)
    require(candidate is not None, f"Required executable is not on PATH: {name}")
    return Path(candidate).resolve()


def read_frozen_order(path: Path) -> tuple[str, ...]:
    require(path.is_file(), f"Missing frozen count table: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ContractError(f"Empty frozen count table: {path}") from exc
    require(header and header[0] in {"#OTUID", "#OTU ID", "OTUID"}, f"Unexpected count header: {header[:1]}")
    samples = tuple(header[1:])
    unique(samples, "frozen sample IDs")
    require(len(samples) == EXPECTED_SAMPLES, f"Expected 75 frozen samples; found {len(samples)}")
    require(not any(sample.startswith("SS_BS") for sample in samples), "SS_BS leaked into frozen sample order")
    return samples


def read_metadata(path: Path) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]]]:
    require(path.is_file(), f"Missing metadata: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise ContractError(f"Empty metadata: {path}") from exc
        require(header == ("SampleID", "site", "Group", "type", "subgroup", "subID"), f"Unexpected metadata header: {header}")
        rows: dict[str, tuple[str, ...]] = {}
        for line_number, row in enumerate(reader, start=2):
            if not row or all(cell == "" for cell in row):
                continue
            require(len(row) == len(header), f"{path}:{line_number} has {len(row)} columns; expected {len(header)}")
            require(row[0] not in rows, f"Duplicate metadata SampleID {row[0]}")
            rows[row[0]] = tuple(row)
    return header, rows


def build_sample_inputs(paths: RunPaths) -> tuple[tuple[str, ...], tuple[str, ...], list[SampleInput]]:
    sample_order = read_frozen_order(paths.frozen_count)
    metadata_header, metadata_rows = read_metadata(paths.metadata)
    require(set(metadata_rows) == set(sample_order), "Metadata sample set differs from frozen 75-sample order")
    require(paths.fastq_dir.is_dir(), f"Missing FASTQ directory: {paths.fastq_dir}")
    expected_names: set[str] = set()
    samples: list[SampleInput] = []
    for sample_id in sample_order:
        r1 = paths.fastq_dir / f"18S_{sample_id}_1.fastq.gz"
        r2 = paths.fastq_dir / f"18S_{sample_id}_2.fastq.gz"
        require(r1.is_file(), f"Missing 18S R1: {r1}")
        require(r2.is_file(), f"Missing 18S R2: {r2}")
        expected_names.update((r1.name, r2.name))
        samples.append(SampleInput(sample_id, r1.resolve(), r2.resolve(), metadata_rows[sample_id]))
    actual_names = {path.name for path in paths.fastq_dir.glob("*.fastq.gz")}
    require(actual_names == expected_names, f"FASTQ set drift: missing={sorted(expected_names-actual_names)[:5]}, extra={sorted(actual_names-expected_names)[:5]}")
    require(len(actual_names) == EXPECTED_FASTQ_FILES, f"Expected 150 FASTQs; found {len(actual_names)}")
    return sample_order, metadata_header, samples


def normalize_fastq_read_id(header: str) -> str:
    identifier = header[1:].split()[0]
    return re.sub(r"/[12]$", "", identifier)


def read_fastq_prefixes(path: Path, limit: int) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            for index in range(limit):
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().strip()
                plus = handle.readline()
                quality = handle.readline().strip()
                require(header.startswith("@"), f"Malformed FASTQ header in {path} at sampled record {index + 1}")
                require(plus.startswith("+"), f"Malformed FASTQ plus line in {path} at sampled record {index + 1}")
                require(len(sequence) == len(quality), f"Sequence/quality length mismatch in {path} at sampled record {index + 1}")
                records.append((normalize_fastq_read_id(header), sequence))
    except (OSError, EOFError) as exc:
        raise ContractError(f"Cannot read gzip FASTQ {path}: {exc}") from exc
    require(records, f"No FASTQ records sampled from {path}")
    return records


def primer_audit(samples: Sequence[SampleInput], records_per_read: int) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_r1 = total_r2 = total_forward = total_reverse = total_wrong16 = 0
    for sample in samples:
        r1_records = read_fastq_prefixes(sample.r1, records_per_read)
        r2_records = read_fastq_prefixes(sample.r2, records_per_read)
        require(len(r1_records) == len(r2_records), f"Sampled R1/R2 record-count mismatch for {sample.sample_id}")
        paired_ids_match = all(r1[0] == r2[0] for r1, r2 in zip(r1_records, r2_records))
        require(paired_ids_match, f"Sampled R1/R2 read IDs do not pair for {sample.sample_id}")
        r1_sequences = [record[1] for record in r1_records]
        r2_sequences = [record[1] for record in r2_records]
        forward_hits = sum(bool(FORWARD_REGEX.match(seq)) for seq in r1_sequences)
        reverse_hits = sum(bool(REVERSE_REGEX.match(seq)) for seq in r2_sequences)
        wrong16_hits = sum(bool(WRONG_16S_REGEX.match(seq)) for seq in r1_sequences)
        total_r1 += len(r1_sequences)
        total_r2 += len(r2_sequences)
        total_forward += forward_hits
        total_reverse += reverse_hits
        total_wrong16 += wrong16_hits
        forward_rate = forward_hits / len(r1_sequences)
        reverse_rate = reverse_hits / len(r2_sequences)
        wrong16_rate = wrong16_hits / len(r1_sequences)
        require(forward_rate >= 0.70, f"Per-sample 18S R1 primer rate too low for {sample.sample_id}: {forward_rate:.3f}")
        require(reverse_rate >= 0.70, f"Per-sample 18S R2 primer rate too low for {sample.sample_id}: {reverse_rate:.3f}")
        require(wrong16_rate <= 0.05, f"Per-sample 16S R1 motif guard failed for {sample.sample_id}: {wrong16_rate:.3f}")
        rows.append({
            "SampleID": sample.sample_id,
            "R1_sampled": len(r1_sequences),
            "R1_18S_primer_hits": forward_hits,
            "R1_18S_primer_rate": forward_rate,
            "R1_16S_motif_hits": wrong16_hits,
            "R2_sampled": len(r2_sequences),
            "R2_18S_primer_hits": reverse_hits,
            "R2_18S_primer_rate": reverse_rate,
            "sampled_paired_read_ids_match": paired_ids_match,
        })
    summary = {
        "records_per_read_requested": records_per_read,
        "r1_sampled": total_r1,
        "r2_sampled": total_r2,
        "r1_18s_primer_hits": total_forward,
        "r1_18s_primer_rate": total_forward / total_r1,
        "r2_18s_primer_hits": total_reverse,
        "r2_18s_primer_rate": total_reverse / total_r2,
        "r1_16s_motif_hits": total_wrong16,
        "r1_16s_motif_rate": total_wrong16 / total_r1,
        "minimum_per_sample_r1_18s_primer_rate": min(float(row["R1_18S_primer_rate"]) for row in rows),
        "minimum_per_sample_r2_18s_primer_rate": min(float(row["R2_18S_primer_rate"]) for row in rows),
        "maximum_per_sample_r1_16s_motif_rate": max(float(row["R1_16S_motif_hits"]) / int(row["R1_sampled"]) for row in rows),
        "all_sampled_paired_read_ids_match": all(bool(row["sampled_paired_read_ids_match"]) for row in rows),
    }
    require(summary["r1_18s_primer_rate"] >= 0.80, f"18S R1 primer rate too low: {summary['r1_18s_primer_rate']:.3f}")
    require(summary["r2_18s_primer_rate"] >= 0.80, f"18S R2 primer rate too low: {summary['r2_18s_primer_rate']:.3f}")
    require(summary["r1_16s_motif_rate"] <= 0.01, f"16S R1 motif contamination guard failed: {summary['r1_16s_motif_rate']:.3f}")
    return rows, summary


def capture_version(command: Sequence[str]) -> str:
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    output = completed.stdout.strip()
    require(output, f"No version output from: {' '.join(command)}")
    return output


def resolve_tools() -> tuple[ToolPaths, dict[str, str]]:
    tools = ToolPaths(
        seqkit=resolve_executable("seqkit"),
        vsearch=resolve_executable("vsearch"),
        usearch12=resolve_executable("usearch_win_12.0-beta"),
        python=Path(sys.executable).resolve(),
    )
    versions = {
        "seqkit": capture_version([str(tools.seqkit), "version"]),
        "vsearch": capture_version([str(tools.vsearch), "--version"]),
        "usearch12": capture_version([str(tools.usearch12)]),
        "python": sys.version,
    }
    require("v12.0" in versions["usearch12"], f"UNOISE3 runner requires USEARCH 12; got {versions['usearch12']!r}")
    return tools, versions


def run_paths(project_root: Path, output_root: Path | None) -> RunPaths:
    selected_output = output_root if output_root is not None else project_root / "revision" / "re_analysis" / RUN_SLUG
    return RunPaths(
        project_root=project_root,
        output_root=selected_output.resolve(),
        frozen_count=(project_root / "revision" / "re_analysis" / "01_input_16s_75_v1" / "16s_asv_counts_all_75_v1.tsv").resolve(),
        metadata=(project_root / "revision" / "re_analysis" / "01_input_16s_75_v1" / "metadata_16s_75_v1.tsv").resolve(),
        fastq_dir=(project_root / "raw_data" / "seq_18s").resolve(),
        pr2_utax=(project_root / "revision" / "re_analysis" / "new_db" / "pr2_version_5.1.1_SSU_UTAX.fasta.gz").resolve(),
        pr2_dada2=(project_root / "revision" / "re_analysis" / "new_db" / "pr2_version_5.1.1_SSU_dada2.fasta.gz").resolve(),
        silva_v123=(project_root / "my_results" / "db" / "silva_18s_v123.fa").resolve(),
        original_script=(project_root / "revision" / "re_analysis" / "ori_analy_script" / "Illumina_pipeline-2024nibr_20260211_18s.sh").resolve(),
        helper_script=(project_root / "revision" / "re_analysis" / "scripts" / "build_18s_parallel_taxonomy_v1.py").resolve(),
    )


def static_input_audit(paths: RunPaths, primer_records: int) -> tuple[tuple[str, ...], tuple[str, ...], list[SampleInput], ToolPaths, dict[str, str], list[dict[str, object]], dict[str, object]]:
    for path in (paths.pr2_utax, paths.pr2_dada2, paths.silva_v123, paths.original_script, paths.helper_script):
        require(path.is_file(), f"Missing required source: {path}")
    sample_order, metadata_header, samples = build_sample_inputs(paths)
    tools, versions = resolve_tools()
    primer_rows, primer_summary = primer_audit(samples, primer_records)
    return sample_order, metadata_header, samples, tools, versions, primer_rows, primer_summary


def stage_marker(output_root: Path, stage_name: str) -> Path:
    return output_root / "_state" / f"{stage_name}.json"


def output_records(paths: Sequence[Path]) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in paths:
        require(path.is_file(), f"Stage output missing: {path}")
        records.append(file_fingerprint(path, include_hash=True))
    return records


def validate_stage_marker(output_root: Path, stage_name: str) -> dict[str, object]:
    marker = stage_marker(output_root, stage_name)
    require(marker.is_file(), f"Missing stage marker: {marker}")
    with marker.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    require(data.get("stage") == stage_name, f"Stage marker name mismatch: {marker}")
    for record in data.get("outputs", []):
        path = Path(record["path"])
        require(path.is_file(), f"Marked output missing: {path}")
        require(path.stat().st_size == record["size_bytes"], f"Marked output size drift: {path}")
        require(sha256_file(path) == record["sha256"], f"Marked output hash drift: {path}")
    return data


def complete_stage(output_root: Path, stage_name: str, outputs: Sequence[Path], started_at: str) -> None:
    record = {
        "stage": stage_name,
        "started_at_utc": started_at,
        "completed_at_utc": utc_now(),
        "outputs": output_records(outputs),
    }
    write_json(stage_marker(output_root, stage_name), record)


def prepare_stage(output_root: Path, stage_name: str, resume: bool) -> bool:
    marker = stage_marker(output_root, stage_name)
    stage_dir = output_root / "stages" / stage_name
    if marker.exists():
        require(resume, f"Completed stage exists; rerun requires --resume: {stage_name}")
        validate_stage_marker(output_root, stage_name)
        print(f"[SKIP] {stage_name} (verified marker)")
        return False
    require(not stage_dir.exists(), f"Unmarked stage directory exists and will not be overwritten: {stage_dir}")
    stage_dir.mkdir(parents=True, exist_ok=False)
    return True


def run_command(command: Sequence[str], log_path: Path, *, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        log.write("COMMAND\t" + json.dumps(list(command), ensure_ascii=False) + "\n")
        log.write("START_UTC\t" + utc_now() + "\n")
        log.flush()
        completed = subprocess.run(command, cwd=cwd, text=True, stdout=log, stderr=subprocess.STDOUT, check=False)
        log.write(f"EXIT_CODE\t{completed.returncode}\n")
        log.write(f"ELAPSED_SECONDS\t{time.monotonic() - started:.3f}\n")
    require(completed.returncode == 0, f"Command failed with exit {completed.returncode}; see {log_path}")


def create_run_root(
    paths: RunPaths,
    sample_order: Sequence[str],
    metadata_header: Sequence[str],
    samples: Sequence[SampleInput],
    tools: ToolPaths,
    versions: dict[str, str],
    primer_rows: Sequence[dict[str, object]],
    primer_summary: dict[str, object],
    args: argparse.Namespace,
) -> None:
    require(not paths.output_root.exists(), f"Output already exists and will not be overwritten: {paths.output_root}")
    paths.output_root.mkdir(parents=True, exist_ok=False)
    for name in ("_state", "logs", "provenance", "qc", "results", "stages"):
        (paths.output_root / name).mkdir()

    config = {
        "runner_version": VERSION,
        "created_at_utc": utc_now(),
        "project_root": str(paths.project_root),
        "output_root": str(paths.output_root),
        "sample_order": list(sample_order),
        "parameters": {
            "forward_primer": FORWARD_PRIMER,
            "reverse_primer": REVERSE_PRIMER,
            "strip_left": 20,
            "strip_right": 18,
            "fastq_maxee_rate": args.maxee_rate,
            "derep_minuniquesize": args.min_size,
            "unoise_minsize": args.min_size,
            "primary_mapping_identity": 0.97,
            "exact_mapping_identity": 1.0,
            "sintax_cutoff": args.sintax_cutoff,
            "minimum_reference_chimera_retention": args.min_chimera_retention,
            "reference_chimera_mode": "diagnostic_only",
            "threads": args.threads,
            "merge_jobs": args.merge_jobs,
            "rarefaction_performed": False,
        },
        "paths": {
            "frozen_count": str(paths.frozen_count),
            "metadata": str(paths.metadata),
            "fastq_dir": str(paths.fastq_dir),
            "pr2_utax": str(paths.pr2_utax),
            "pr2_dada2": str(paths.pr2_dada2),
            "silva_v123": str(paths.silva_v123),
            "original_script": str(paths.original_script),
            "helper_script": str(paths.helper_script),
        },
        "tools": {
            "seqkit": str(tools.seqkit),
            "vsearch": str(tools.vsearch),
            "usearch12": str(tools.usearch12),
            "python": str(tools.python),
        },
        "tool_versions": versions,
        "intentional_corrections": [
            "Pin merge input to raw_data/seq_18s rather than legacy seq/",
            "Strip merged reads left 20 and right 18 based on observed primers rather than legacy left 18/right 20",
            "Run PR2 v5.1.1 UCHIME reference classification as diagnostic-only QC; do not hard-filter UNOISE3 ASVs because PR2 documents the UTAX file for taxonomy assignment, not as a validated UCHIME reference",
            "Retain Chloroplastida and Metazoa as flag-only eukaryotic lineages",
            "Do not rarefy before reviewing new sample depths",
        ],
    }
    write_json(paths.output_root / "provenance" / "run_config_v1.json", config)

    metadata_rows = [sample.metadata for sample in samples]
    write_tsv(paths.output_root / "provenance" / "metadata_18s_75_v1.tsv", metadata_header, metadata_rows)
    write_tsv(
        paths.output_root / "qc" / "primer_audit_18s_v1.tsv",
        ("SampleID", "R1_sampled", "R1_18S_primer_hits", "R1_18S_primer_rate", "R1_16S_motif_hits", "R2_sampled", "R2_18S_primer_hits", "R2_18S_primer_rate", "sampled_paired_read_ids_match"),
        (
            (
                row["SampleID"], row["R1_sampled"], row["R1_18S_primer_hits"], f"{row['R1_18S_primer_rate']:.6f}",
                row["R1_16S_motif_hits"], row["R2_sampled"], row["R2_18S_primer_hits"], f"{row['R2_18S_primer_rate']:.6f}",
                row["sampled_paired_read_ids_match"],
            )
            for row in primer_rows
        ),
    )
    write_json(paths.output_root / "qc" / "primer_audit_summary_18s_v1.json", primer_summary)

    source_files: list[tuple[str, Path]] = [
        ("frozen_count", paths.frozen_count),
        ("metadata", paths.metadata),
        ("pr2_utax", paths.pr2_utax),
        ("pr2_dada2", paths.pr2_dada2),
        ("silva_v123", paths.silva_v123),
        ("original_script", paths.original_script),
        ("runner", Path(__file__).resolve()),
        ("taxonomy_helper", paths.helper_script),
    ]
    for sample in samples:
        source_files.extend(((f"{sample.sample_id}_R1", sample.r1), (f"{sample.sample_id}_R2", sample.r2)))
    manifest_rows = []
    for role, path in source_files:
        fingerprint = file_fingerprint(path, include_hash=True)
        manifest_rows.append((role, fingerprint["path"], fingerprint["size_bytes"], fingerprint["mtime_ns"], fingerprint["sha256"]))
    write_tsv(paths.output_root / "provenance" / "source_sha256_manifest_v1.tsv", ("role", "absolute_path", "size_bytes", "mtime_ns", "sha256"), manifest_rows)


def verify_source_manifest(path: Path) -> None:
    require(path.is_file(), f"Missing source manifest: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(
            reader.fieldnames == ["role", "absolute_path", "size_bytes", "mtime_ns", "sha256"],
            f"Unexpected source-manifest header: {reader.fieldnames}",
        )
        rows = list(reader)
    require(len(rows) == 158, f"Expected 158 manifested sources; found {len(rows)}")
    for row in rows:
        source = Path(row["absolute_path"])
        require(source.is_file(), f"Manifested source missing: {source}")
        require(source.stat().st_size == int(row["size_bytes"]), f"Manifested source size drift: {source}")
        require(sha256_file(source) == row["sha256"], f"Manifested source hash drift: {source}")


def verify_resume_config(
    paths: RunPaths,
    sample_order: Sequence[str],
    tools: ToolPaths,
    versions: dict[str, str],
    args: argparse.Namespace,
) -> None:
    config_path = paths.output_root / "provenance" / "run_config_v1.json"
    require(config_path.is_file(), f"Resume output is not runner-owned; missing {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    require(config.get("runner_version") == VERSION, "Runner version differs from existing output")
    require(config.get("sample_order") == list(sample_order), "Frozen sample order differs from existing run")
    params = config.get("parameters", {})
    expected = {
        "fastq_maxee_rate": args.maxee_rate,
        "derep_minuniquesize": args.min_size,
        "unoise_minsize": args.min_size,
        "sintax_cutoff": args.sintax_cutoff,
        "minimum_reference_chimera_retention": args.min_chimera_retention,
        "threads": args.threads,
        "merge_jobs": args.merge_jobs,
    }
    for key, value in expected.items():
        require(params.get(key) == value, f"Resume parameter drift for {key}: existing={params.get(key)!r}, requested={value!r}")
    configured_paths = config.get("paths", {})
    for key in ("frozen_count", "metadata", "fastq_dir", "pr2_utax", "pr2_dada2", "silva_v123", "original_script", "helper_script"):
        require(configured_paths.get(key) == str(getattr(paths, key)), f"Resume path drift for {key}")
    configured_tools = config.get("tools", {})
    for key in ("seqkit", "vsearch", "usearch12", "python"):
        require(configured_tools.get(key) == str(getattr(tools, key)), f"Resume tool-path drift for {key}")
    require(config.get("tool_versions") == versions, "Resume tool-version drift")
    verify_source_manifest(paths.output_root / "provenance" / "source_sha256_manifest_v1.tsv")


def stage_seqkit(paths: RunPaths, samples: Sequence[SampleInput], tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "01_seqkit"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    output = stage_dir / "seqkit_stats_18s_v1.tsv"
    command = [str(tools.seqkit), "stats", "-T", "-a", "-b", "-j", str(min(args.threads, 12)), "-o", str(output)]
    for sample in samples:
        command.extend((str(sample.r1), str(sample.r2)))
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    complete_stage(paths.output_root, name, (output,), started)


def merge_one(sample: SampleInput, stage_dir: Path, vsearch: Path) -> Path:
    output = stage_dir / f"{sample.sample_id}.merged.fastq"
    log_path = stage_dir / "logs" / f"{sample.sample_id}.log"
    command = [
        str(vsearch), "--fastq_mergepairs", str(sample.r1), "--reverse", str(sample.r2),
        "--fastqout", str(output), "--relabel", f"{sample.sample_id}.",
    ]
    run_command(command, log_path)
    return output


def stage_merge(paths: RunPaths, samples: Sequence[SampleInput], tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "02_merge"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    (stage_dir / "logs").mkdir()
    futures = {}
    with ThreadPoolExecutor(max_workers=args.merge_jobs) as pool:
        for sample in samples:
            futures[pool.submit(merge_one, sample, stage_dir, tools.vsearch)] = sample.sample_id
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                future.result()
            except Exception as exc:
                raise ContractError(f"Merge failed for {sample_id}: {exc}") from exc
            print(f"[MERGED] {sample_id}")
    merged_paths = [stage_dir / f"{sample.sample_id}.merged.fastq" for sample in samples]
    combined = stage_dir / "all_18s_75.merged.fastq"
    with combined.open("wb") as destination:
        for merged in merged_paths:
            require(merged.is_file() and merged.stat().st_size > 0, f"Missing/empty merged FASTQ: {merged}")
            with merged.open("rb") as source:
                shutil.copyfileobj(source, destination, length=16 * 1024 * 1024)
    merge_logs = [stage_dir / "logs" / f"{sample.sample_id}.log" for sample in samples]
    complete_stage(paths.output_root, name, (*merged_paths, combined, *merge_logs), started)


def stage_filter(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "03_filter"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    source = paths.output_root / "stages" / "02_merge" / "all_18s_75.merged.fastq"
    output = stage_dir / "filtered_18s_v1.fasta"
    command = [
        str(tools.vsearch), "--fastx_filter", str(source),
        "--fastq_stripleft", "20", "--fastq_stripright", "18",
        "--fastq_maxee_rate", str(args.maxee_rate), "--fastaout", str(output),
    ]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    complete_stage(paths.output_root, name, (output,), started)


def stage_derep(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "04_derep"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    source = paths.output_root / "stages" / "03_filter" / "filtered_18s_v1.fasta"
    output = stage_dir / "uniques_18s_v1.fasta"
    command = [
        str(tools.vsearch), "--derep_fulllength", str(source), "--minuniquesize", str(args.min_size),
        "--sizeout", "--relabel", "Uni_", "--output", str(output),
    ]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    complete_stage(paths.output_root, name, (output,), started)


def relabel_zotus(source: Path, output: Path) -> int:
    count = 0
    with source.open("r", encoding="ascii") as reader, output.open("w", encoding="ascii", newline="\n") as writer:
        for line in reader:
            if line.startswith(">"):
                match = re.match(r">Zotu_?(\d+)(.*)$", line.rstrip("\r\n"), re.IGNORECASE)
                require(match is not None, f"Unexpected UNOISE3 header: {line.strip()}")
                count += 1
                writer.write(f">18S_ASV_{match.group(1)}{match.group(2)}\n")
            else:
                writer.write(line.replace("\r\n", "\n"))
    require(count > 0, "UNOISE3 produced zero ZOTUs")
    return count


def stage_unoise(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "05_unoise"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    source = paths.output_root / "stages" / "04_derep" / "uniques_18s_v1.fasta"
    zotus = stage_dir / "zotus_18s_v1.fasta"
    asvs = stage_dir / "asvs_18s_v1.fasta"
    command = [str(tools.usearch12), "-unoise3", str(source), "-minsize", str(args.min_size), "-zotus", str(zotus)]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    count = relabel_zotus(zotus, asvs)
    write_json(stage_dir / "unoise_summary_v1.json", {"asv_count": count, "minsize": args.min_size})
    complete_stage(paths.output_root, name, (zotus, asvs, stage_dir / "unoise_summary_v1.json"), started)


def stage_chimera(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace) -> None:
    name = "06_chimera"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    source = paths.output_root / "stages" / "05_unoise" / "asvs_18s_v1.fasta"
    nonchimeras = stage_dir / "18s_pr2_uchime_nonchimeras_qc_v1.fasta"
    chimeras = stage_dir / "18s_pr2_uchime_chimeras_qc_v1.fasta"
    borderline = stage_dir / "18s_pr2_uchime_borderline_qc_v1.fasta"
    uchimeout = stage_dir / "18s_uchime_ref_v1.tsv"
    command = [
        str(tools.vsearch), "--uchime_ref", str(source), "--db", str(paths.pr2_utax),
        "--nonchimeras", str(nonchimeras), "--chimeras", str(chimeras),
        "--borderline", str(borderline),
        "--uchimeout", str(uchimeout), "--threads", str(args.threads),
    ]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    for path in (nonchimeras, chimeras, borderline, uchimeout):
        if not path.exists():
            path.touch()
    require(nonchimeras.stat().st_size > 0, "PR2 UCHIME diagnostic classified zero ASVs as nonchimeric")
    input_order = fasta_ids(source)
    input_ids = set(input_order)
    retained_ids = set(fasta_ids(nonchimeras))
    chimera_ids = set(fasta_ids(chimeras, allow_empty=True))
    borderline_ids = set(fasta_ids(borderline, allow_empty=True))
    require(not (retained_ids & chimera_ids or retained_ids & borderline_ids or chimera_ids & borderline_ids), "Reference-chimera output ID partitions overlap")
    require(retained_ids | chimera_ids | borderline_ids == input_ids, "Reference-chimera ID partitions do not reconcile to UNOISE3 input")
    retention_rate = len(retained_ids) / len(input_ids)
    review_gate_passed = retention_rate >= args.min_chimera_retention
    status_by_id = {feature_id: "nonchimera" for feature_id in retained_ids}
    status_by_id.update({feature_id: "chimera" for feature_id in chimera_ids})
    status_by_id.update({feature_id: "borderline" for feature_id in borderline_ids})
    ledger = stage_dir / "18s_pr2_uchime_diagnostic_ledger_v1.tsv"
    write_tsv(
        ledger,
        ("ASV_ID", "pr2_uchime_ref_status", "primary_unoise3_keep"),
        ((feature_id, status_by_id[feature_id], "TRUE") for feature_id in input_order),
    )
    summary = stage_dir / "reference_chimera_summary_v1.json"
    write_json(summary, {
        "input_unoise_asvs": len(input_ids),
        "retained_nonchimeras": len(retained_ids),
        "reference_chimeras": len(chimera_ids),
        "borderline_excluded": len(borderline_ids),
        "retention_rate": retention_rate,
        "minimum_retention_gate": args.min_chimera_retention,
        "review_gate_passed": review_gate_passed,
        "reference_database": str(paths.pr2_utax),
        "mode": "diagnostic_only",
        "primary_filter_applied": False,
        "primary_asvs_after_unoise3": len(input_ids),
    })
    complete_stage(paths.output_root, name, (nonchimeras, chimeras, borderline, uchimeout, ledger, summary), started)


def stage_map(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace, *, identity: float, name: str, filename: str) -> None:
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    reads = paths.output_root / "stages" / "03_filter" / "filtered_18s_v1.fasta"
    database = paths.output_root / "stages" / "05_unoise" / "asvs_18s_v1.fasta"
    output = stage_dir / filename
    command = [
        str(tools.vsearch), "--usearch_global", str(reads), "--db", str(database),
        "--id", f"{identity:.2f}", "--threads", str(args.threads), "--otutabout", str(output),
    ]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    complete_stage(paths.output_root, name, (output,), started)


def stage_taxonomy(paths: RunPaths, tools: ToolPaths, args: argparse.Namespace, *, name: str, database: Path, filename: str) -> None:
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    query = paths.output_root / "stages" / "05_unoise" / "asvs_18s_v1.fasta"
    output = stage_dir / filename
    command = [
        str(tools.vsearch), "--sintax", str(query), "--db", str(database),
        "--sintax_cutoff", str(args.sintax_cutoff), "--tabbedout", str(output),
        "--threads", str(args.threads),
    ]
    run_command(command, paths.output_root / "logs" / f"{name}.log")
    complete_stage(paths.output_root, name, (output,), started)


def parse_seqkit_counts(path: Path) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None and "file" in reader.fieldnames and "num_seqs" in reader.fieldnames, f"Unexpected SeqKit table: {path}")
        for row in reader:
            filename = Path(row["file"]).name
            match = re.fullmatch(r"18S_(.+)_([12])\.fastq\.gz", filename)
            require(match is not None, f"Cannot parse SeqKit FASTQ name: {filename}")
            counts[(match.group(1), f"R{match.group(2)}")] = int(row["num_seqs"].replace(",", ""))
    return counts


MERGED_RE = re.compile(r"^\s*(\d+)\s+Merged\s+\(", re.MULTILINE)
PAIRS_RE = re.compile(r"^\s*(\d+)\s+Pairs\s*$", re.MULTILINE)


def parse_merge_log(path: Path) -> tuple[int, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    pairs_match = PAIRS_RE.search(text)
    merged_match = MERGED_RE.search(text)
    require(pairs_match is not None and merged_match is not None, f"Cannot parse VSEARCH merge log: {path}")
    return int(pairs_match.group(1)), int(merged_match.group(1))


def fasta_sample_counts(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with path.open("r", encoding="ascii", errors="strict") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            identifier = line[1:].split()[0]
            sample_id = identifier.split(".", 1)[0]
            counts[sample_id] += 1
    return counts


def count_table_sample_sums(path: Path) -> tuple[tuple[str, ...], dict[str, int], tuple[str, ...]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ContractError(f"Empty count table: {path}") from exc
        samples = tuple(header[1:])
        sums = [0] * len(samples)
        feature_ids: list[str] = []
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            if not row:
                continue
            require(len(row) == len(header), f"{path}:{line_number} has {len(row)} columns; expected {len(header)}")
            require(row[0] not in seen, f"Duplicate feature ID in {path}: {row[0]}")
            seen.add(row[0])
            feature_ids.append(row[0])
            for index, value in enumerate(row[1:]):
                integer = int(value)
                require(integer >= 0, f"Negative count in {path}:{line_number}")
                sums[index] += integer
    return samples, dict(zip(samples, sums)), tuple(feature_ids)


def fasta_ids(path: Path, *, allow_empty: bool = False) -> tuple[str, ...]:
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="ascii") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            identifier = line[1:].split()[0]
            require(identifier and identifier not in seen, f"Duplicate/empty FASTA ID at {path}:{line_number}")
            seen.add(identifier)
            identifiers.append(identifier)
    require(allow_empty or identifiers, f"No FASTA records in {path}")
    return tuple(identifiers)


def sintax_ids(path: Path) -> tuple[str, ...]:
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            require(len(row) == 4, f"{path}:{line_number} has {len(row)} SINTAX fields; expected 4")
            require(row[0] and row[0] not in seen, f"Duplicate/empty SINTAX ID at {path}:{line_number}")
            seen.add(row[0])
            identifiers.append(row[0])
    return tuple(identifiers)


def reorder_count_table(source: Path, destination: Path, frozen_samples: Sequence[str]) -> tuple[int, int]:
    with source.open("r", encoding="utf-8-sig", newline="") as reader_handle:
        reader = csv.reader(reader_handle, delimiter="\t")
        header = next(reader)
        source_samples = tuple(header[1:])
        require(set(source_samples) == set(frozen_samples), f"Count sample set differs from frozen order in {source}")
        indices = [source_samples.index(sample) + 1 for sample in frozen_samples]
        feature_count = 0
        total = 0
        with destination.open("w", encoding="utf-8", newline="") as writer_handle:
            writer = csv.writer(writer_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#OTUID", *frozen_samples])
            for line_number, row in enumerate(reader, start=2):
                if not row:
                    continue
                require(len(row) == len(header), f"{source}:{line_number} malformed")
                values = [int(row[index]) for index in indices]
                require(all(value >= 0 for value in values), f"Negative count in {source}:{line_number}")
                writer.writerow([row[0], *values])
                feature_count += 1
                total += sum(values)
    return feature_count, total


def copy_core_results(paths: RunPaths, sample_order: Sequence[str]) -> tuple[Path, ...]:
    results = paths.output_root / "results"
    counts97 = results / "18s_asv_counts_97_all_v1.tsv"
    counts100 = results / "18s_asv_counts_exact_all_v1.tsv"
    reorder_count_table(paths.output_root / "stages" / "07_map97" / "18s_asv_counts_97_raw_v1.tsv", counts97, sample_order)
    reorder_count_table(paths.output_root / "stages" / "08_map100" / "18s_asv_counts_exact_raw_v1.tsv", counts100, sample_order)
    rep = results / "18s_rep_seqs_all_v1.fasta"
    pr2 = results / "18s_sintax_pr2_v5.1.1_all_v1.tsv"
    silva = results / "18s_sintax_silva_v123_all_v1.tsv"
    shutil.copy2(paths.output_root / "stages" / "05_unoise" / "asvs_18s_v1.fasta", rep)
    shutil.copy2(paths.output_root / "stages" / "09_tax_pr2" / "18s_sintax_pr2_v5.1.1_raw_v1.tsv", pr2)
    shutil.copy2(paths.output_root / "stages" / "10_tax_silva" / "18s_sintax_silva_v123_raw_v1.tsv", silva)
    return counts97, counts100, rep, pr2, silva


def stage_core_qc(paths: RunPaths, sample_order: Sequence[str], samples: Sequence[SampleInput], args: argparse.Namespace) -> None:
    name = "11_core_qc"
    if not prepare_stage(paths.output_root, name, args.resume):
        return
    started = utc_now()
    stage_dir = paths.output_root / "stages" / name
    core_results = copy_core_results(paths, sample_order)
    chimera_ledger_qc = paths.output_root / "qc" / "18s_pr2_uchime_diagnostic_ledger_v1.tsv"
    chimera_summary_qc = paths.output_root / "qc" / "18s_pr2_uchime_diagnostic_summary_v1.json"
    shutil.copy2(
        paths.output_root / "stages" / "06_chimera" / "18s_pr2_uchime_diagnostic_ledger_v1.tsv",
        chimera_ledger_qc,
    )
    shutil.copy2(
        paths.output_root / "stages" / "06_chimera" / "reference_chimera_summary_v1.json",
        chimera_summary_qc,
    )
    seqkit_counts = parse_seqkit_counts(paths.output_root / "stages" / "01_seqkit" / "seqkit_stats_18s_v1.tsv")
    filtered_counts = fasta_sample_counts(paths.output_root / "stages" / "03_filter" / "filtered_18s_v1.fasta")
    map97_samples, map97_sums, ids97 = count_table_sample_sums(core_results[0])
    map100_samples, map100_sums, ids100 = count_table_sample_sums(core_results[1])
    require(map97_samples == tuple(sample_order), "97% count table is not in frozen sample order")
    require(map100_samples == tuple(sample_order), "Exact count table is not in frozen sample order")
    require(set(ids97) == set(ids100), "97% and exact mapping feature sets differ")
    rep_ids = fasta_ids(core_results[2])
    pr2_ids = sintax_ids(core_results[3])
    silva_ids = sintax_ids(core_results[4])
    require(set(ids97) == set(rep_ids), "Count/representative-sequence feature sets differ")
    require(set(ids97) == set(pr2_ids), "Count/PR2-SINTAX feature sets differ")
    require(set(ids97) == set(silva_ids), "Count/SILVA-SINTAX feature sets differ")
    rows = []
    for sample in samples:
        raw_r1 = seqkit_counts[(sample.sample_id, "R1")]
        raw_r2 = seqkit_counts[(sample.sample_id, "R2")]
        require(raw_r1 == raw_r2, f"Raw R1/R2 count mismatch for {sample.sample_id}")
        pairs, merged = parse_merge_log(paths.output_root / "stages" / "02_merge" / "logs" / f"{sample.sample_id}.log")
        require(pairs == raw_r1, f"Merge input count differs from SeqKit for {sample.sample_id}")
        filtered = filtered_counts.get(sample.sample_id, 0)
        mapped97 = map97_sums[sample.sample_id]
        mapped100 = map100_sums[sample.sample_id]
        require(raw_r1 >= merged >= filtered >= mapped97 >= 0, f"Non-monotonic primary retention for {sample.sample_id}")
        require(filtered >= mapped100 >= 0, f"Non-monotonic exact retention for {sample.sample_id}")
        require(mapped97 >= mapped100, f"Exact mapped reads exceed 97% mapped reads for {sample.sample_id}")
        require(mapped97 > 0, f"Sample has zero primary mapped reads: {sample.sample_id}")
        rows.append((
            sample.sample_id, raw_r1, merged, f"{merged/raw_r1:.6f}", filtered, f"{filtered/raw_r1:.6f}",
            mapped97, f"{mapped97/raw_r1:.6f}", mapped100, f"{mapped100/raw_r1:.6f}",
            sample.metadata[1], sample.metadata[2], sample.metadata[3], sample.metadata[4], sample.metadata[5],
        ))
    retention = paths.output_root / "qc" / "read_retention_by_sample_18s_v1.tsv"
    write_tsv(
        retention,
        ("SampleID", "raw_pairs", "merged_reads", "merged_rate", "ee_filtered_reads", "ee_filtered_rate",
         "mapped_97_reads", "mapped_97_rate", "mapped_exact_reads", "mapped_exact_rate",
         "site", "Group", "type", "subgroup", "subID"),
        rows,
    )
    summary = {
        "samples": len(sample_order),
        "unoise3_de_novo_asvs": len(ids97),
        "pr2_reference_chimera_mode": "diagnostic_only",
        "raw_pairs": sum(row[1] for row in rows),
        "merged_reads": sum(row[2] for row in rows),
        "ee_filtered_reads": sum(row[4] for row in rows),
        "mapped_97_reads": sum(row[6] for row in rows),
        "mapped_exact_reads": sum(row[8] for row in rows),
        "rarefaction_performed": False,
    }
    summary_path = paths.output_root / "qc" / "core_pipeline_summary_18s_v1.json"
    write_json(summary_path, summary)
    marker_note = stage_dir / "core_qc_complete_v1.json"
    write_json(marker_note, summary)
    complete_stage(
        paths.output_root,
        name,
        (*core_results, chimera_ledger_qc, chimera_summary_qc, retention, summary_path, marker_note),
        started,
    )


def verify_final_manifest(path: Path, output_root: Path) -> None:
    require(path.is_file(), f"Missing final manifest: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames == ["relative_path", "size_bytes", "sha256"], f"Unexpected final-manifest header: {reader.fieldnames}")
        rows = list(reader)
    require(rows, f"Empty final manifest: {path}")
    for row in rows:
        target = output_root / Path(row["relative_path"])
        require(target.is_file(), f"Manifest target missing: {target}")
        require(target.stat().st_size == int(row["size_bytes"]), f"Manifest target size drift: {target}")
        require(sha256_file(target) == row["sha256"], f"Manifest target hash drift: {target}")


def final_manifest(paths: RunPaths) -> Path:
    output = paths.output_root / "provenance" / "final_sha256_manifest_core_v1.tsv"
    if output.exists():
        verify_final_manifest(output, paths.output_root)
        return output
    rows = []
    for folder in ("results", "qc", "provenance", "_state", "logs"):
        for path in sorted((paths.output_root / folder).rglob("*")):
            if path.is_file() and path != output:
                rows.append((path.relative_to(paths.output_root).as_posix(), path.stat().st_size, sha256_file(path)))
    write_tsv(output, ("relative_path", "size_bytes", "sha256"), rows)
    return output


def verify_core(paths: RunPaths) -> dict[str, object]:
    require(paths.output_root.is_dir(), f"Output package does not exist: {paths.output_root}")
    for stage_name in STAGE_NAMES:
        validate_stage_marker(paths.output_root, stage_name)
    config_path = paths.output_root / "provenance" / "run_config_v1.json"
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    sample_order = tuple(config["sample_order"])
    require(len(sample_order) == EXPECTED_SAMPLES and len(set(sample_order)) == EXPECTED_SAMPLES, "Config sample order invalid")
    results = paths.output_root / "results"
    count97 = results / "18s_asv_counts_97_all_v1.tsv"
    count100 = results / "18s_asv_counts_exact_all_v1.tsv"
    samples97, sums97, ids97 = count_table_sample_sums(count97)
    samples100, sums100, ids100 = count_table_sample_sums(count100)
    require(samples97 == sample_order == samples100, "Final count-table sample order drift")
    require(set(ids97) == set(ids100), "Final mapping feature-set drift")
    require(all(feature.startswith("18S_ASV_") for feature in ids97), "Non-18S feature ID in final counts")
    rep_ids = fasta_ids(results / "18s_rep_seqs_all_v1.fasta")
    pr2_ids = sintax_ids(results / "18s_sintax_pr2_v5.1.1_all_v1.tsv")
    silva_ids = sintax_ids(results / "18s_sintax_silva_v123_all_v1.tsv")
    require(set(ids97) == set(rep_ids), "Final count/representative-sequence feature-set drift")
    require(set(ids97) == set(pr2_ids), "Final count/PR2-SINTAX feature-set drift")
    require(set(ids97) == set(silva_ids), "Final count/SILVA-SINTAX feature-set drift")
    require(sum(sums97.values()) >= sum(sums100.values()), "Exact mapped count exceeds 97% mapped count")
    require(not any(path.name.startswith("otutab_rare") for path in paths.output_root.rglob("*")), "Unexpected rarefied artifact")
    manifest = paths.output_root / "provenance" / "final_sha256_manifest_core_v1.tsv"
    require(manifest.is_file(), f"Missing completed core manifest: {manifest}")
    verify_final_manifest(manifest, paths.output_root)
    return {
        "samples": len(sample_order),
        "asvs": len(ids97),
        "mapped_97_reads": sum(sums97.values()),
        "mapped_exact_reads": sum(sums100.values()),
        "rarefaction_performed": False,
    }


def print_dry_run(paths: RunPaths, sample_order: Sequence[str], versions: dict[str, str], primer_summary: dict[str, object], args: argparse.Namespace) -> None:
    plan = {
        "status": "DRY_RUN_OK",
        "runner_version": VERSION,
        "output_would_be": str(paths.output_root),
        "samples": len(sample_order),
        "fastq_files": EXPECTED_FASTQ_FILES,
        "first_sample": sample_order[0],
        "last_sample": sample_order[-1],
        "primer_summary": primer_summary,
        "parameters": {
            "strip_left": 20,
            "strip_right": 18,
            "maxee_rate": args.maxee_rate,
            "min_size": args.min_size,
            "map_identities": [0.97, 1.0],
            "sintax_cutoff": args.sintax_cutoff,
            "minimum_reference_chimera_retention": args.min_chimera_retention,
            "reference_chimera_db": str(paths.pr2_utax),
            "reference_chimera_mode": "diagnostic_only_no_primary_filter",
            "primary_taxonomy_db": str(paths.pr2_utax),
            "comparison_taxonomy_db": str(paths.silva_v123),
            "rarefaction": "NOT_RUN",
        },
        "tool_versions": versions,
        "stages": list(STAGE_NAMES),
    }
    print(json.dumps(plan, ensure_ascii=False, indent=2))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run or verify the 75-sample 18S V4 UNOISE3 + PR2 v5.1.1 reanalysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true", help="Validate inputs, tools, sample order, and primers without creating output.")
    mode.add_argument("--run", action="store_true", help="Execute the checkpointed core pipeline.")
    mode.add_argument("--check", action="store_true", help="Verify an existing completed core output package.")
    parser.add_argument("--project-root", type=Path, default=project_root_from_script(), help="Project root.")
    parser.add_argument("--output", type=Path, default=None, help="Versioned output root.")
    parser.add_argument("--resume", action="store_true", help="Resume only marker-verified completed stages under an identical contract.")
    parser.add_argument("--threads", type=int, default=20, help="VSEARCH worker threads.")
    parser.add_argument("--merge-jobs", type=int, default=2, help="Concurrent per-sample merge jobs.")
    parser.add_argument("--maxee-rate", type=float, default=0.01, help="VSEARCH expected-error rate threshold.")
    parser.add_argument("--min-size", type=int, default=10, help="Dereplication and UNOISE3 minimum abundance.")
    parser.add_argument("--sintax-cutoff", type=float, default=0.1, help="SINTAX reporting cutoff.")
    parser.add_argument("--min-chimera-retention", type=float, default=0.8, help="Diagnostic PR2 UCHIME retention review threshold; recorded but not used to filter or fail the primary UNOISE3 ASV set.")
    parser.add_argument("--primer-records", type=int, default=100, help="Records sampled per R1/R2 file for assay guard.")
    args = parser.parse_args(argv)
    require(args.threads >= 1, "--threads must be >=1")
    require(1 <= args.merge_jobs <= 8, "--merge-jobs must be between 1 and 8")
    require(0 < args.maxee_rate <= 1, "--maxee-rate must be in (0,1]")
    require(args.min_size >= 2, "--min-size must be >=2")
    require(0 <= args.sintax_cutoff <= 1, "--sintax-cutoff must be in [0,1]")
    require(0 < args.min_chimera_retention <= 1, "--min-chimera-retention must be in (0,1]")
    require(args.primer_records >= 10, "--primer-records must be >=10")
    require(not args.resume or args.run, "--resume is valid only with --run")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    project_root = args.project_root.resolve()
    paths = run_paths(project_root, args.output.resolve() if args.output is not None else None)

    if args.check:
        summary = verify_core(paths)
        print(json.dumps({"status": "CHECK_OK", **summary}, ensure_ascii=False, indent=2))
        return 0

    audit = static_input_audit(paths, args.primer_records)
    sample_order, metadata_header, samples, tools, versions, primer_rows, primer_summary = audit
    if args.dry_run:
        require(not paths.output_root.exists(), f"Dry-run target already exists; choose a new version or use --check: {paths.output_root}")
        print_dry_run(paths, sample_order, versions, primer_summary, args)
        return 0

    if args.resume:
        require(paths.output_root.is_dir(), f"Cannot resume missing output: {paths.output_root}")
        verify_resume_config(paths, sample_order, tools, versions, args)
    else:
        create_run_root(paths, sample_order, metadata_header, samples, tools, versions, primer_rows, primer_summary, args)

    print(f"[RUN] output={paths.output_root}")
    stage_seqkit(paths, samples, tools, args)
    stage_merge(paths, samples, tools, args)
    stage_filter(paths, tools, args)
    stage_derep(paths, tools, args)
    stage_unoise(paths, tools, args)
    stage_chimera(paths, tools, args)
    stage_map(paths, tools, args, identity=0.97, name="07_map97", filename="18s_asv_counts_97_raw_v1.tsv")
    stage_map(paths, tools, args, identity=1.0, name="08_map100", filename="18s_asv_counts_exact_raw_v1.tsv")
    stage_taxonomy(paths, tools, args, name="09_tax_pr2", database=paths.pr2_utax, filename="18s_sintax_pr2_v5.1.1_raw_v1.tsv")
    stage_taxonomy(paths, tools, args, name="10_tax_silva", database=paths.silva_v123, filename="18s_sintax_silva_v123_raw_v1.tsv")
    stage_core_qc(paths, sample_order, samples, args)
    manifest = final_manifest(paths)
    summary = verify_core(paths)
    print(json.dumps({"status": "CORE_RUN_OK", "manifest": str(manifest), **summary}, ensure_ascii=False, indent=2))
    print("[NEXT] Run build_18s_parallel_taxonomy_v1.py after reviewing its dry-run.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ContractError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
