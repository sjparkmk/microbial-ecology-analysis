#!/usr/bin/env python3
"""Read-only public-release audit for the paper-specific code package."""

from __future__ import annotations

import argparse
import ast
import hashlib
import re
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PACKAGE_ROOT.parents[1]
MANIFEST = PACKAGE_ROOT / "provenance" / "public_file_manifest_v1.tsv"

TEXT_EXTENSIONS = {
    ".R", ".py", ".tsv", ".json", ".yml", ".yaml", ".md", ".txt", ".cff"
}
ALLOWED_NO_EXTENSION = {"LICENSE"}
EXCLUDED_DIR_NAMES = {
    "data", "raw_data", "inputs", "references", "output", "outputs", "results",
    "generated", "cache", "caches", "__pycache__", ".git"
}
EXCLUDED_SUFFIXES = {
    ".fastq", ".fq", ".fasta", ".fa", ".fna", ".bam", ".sam", ".cram",
    ".zip", ".7z", ".tar", ".gz", ".exe", ".dll", ".so", ".dylib", ".pkl",
    ".rds", ".RData", ".sqlite", ".db"
}
PATH_PATTERNS = [
    re.compile(r"(?i)C:\\Users\\|[A-Z]:/Users/|/Users/|/home/|/mnt/"),
    re.compile(r"(?i)/user/link/path|Program Files|AppData|OneDrive"),
]
SECRET_PATTERNS = [
    re.compile(r"(?i)-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(r"(?i)\b(?:api[_-]?key|secret|password|passwd|token)\s*[:=]\s*['\"]?[A-Za-z0-9_./+\-]{16,}"),
    re.compile(r"\b(?:ghp|github_pat|sk|AKIA)[A-Za-z0-9_\-]{12,}\b"),
    re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return path.relative_to(PACKAGE_ROOT).as_posix()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig", errors="replace")


def package_files() -> list[Path]:
    return sorted(
        path for path in PACKAGE_ROOT.rglob("*")
        if path.is_file() and path.name != MANIFEST.name
    )


def check_manifest(files: list[Path], errors: list[str]) -> None:
    if not MANIFEST.is_file():
        errors.append(f"missing manifest: {rel(MANIFEST)}")
        return
    rows = MANIFEST.read_text(encoding="utf-8-sig").splitlines()
    if not rows or rows[0] != "Path\tBytes\tSHA256":
        errors.append("manifest header must be Path\\tBytes\\tSHA256")
        return
    observed: dict[str, tuple[int, str]] = {}
    for line_number, line in enumerate(rows[1:], start=2):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            errors.append(f"manifest line {line_number}: expected 3 tab-separated fields")
            continue
        path_text, byte_text, hash_text = fields
        relative_path = Path(path_text)
        path = PACKAGE_ROOT / relative_path
        if relative_path.is_absolute() or ".." in relative_path.parts:
            errors.append(f"manifest line {line_number}: non-relative path {path_text}")
            continue
        if not path.is_file():
            errors.append(f"manifest line {line_number}: missing file {path_text}")
            continue
        try:
            expected_bytes = int(byte_text)
        except ValueError:
            errors.append(f"manifest line {line_number}: invalid byte count {byte_text}")
            continue
        observed[path_text] = (expected_bytes, hash_text.lower())
        actual_bytes = path.stat().st_size
        actual_hash = sha256(path)
        if actual_bytes != expected_bytes or actual_hash != hash_text.lower():
            errors.append(f"manifest drift: {path_text}")
    expected_paths = {rel(path) for path in files}
    if set(observed) != expected_paths:
        missing = sorted(expected_paths - set(observed))
        extra = sorted(set(observed) - expected_paths)
        if missing:
            errors.append("manifest missing: " + ", ".join(missing))
        if extra:
            errors.append("manifest has unexpected entries: " + ", ".join(extra))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="run the read-only audit")
    args = parser.parse_args()
    if not args.check:
        parser.error("use --check")

    errors: list[str] = []
    warnings: list[str] = []
    files = package_files()

    for path in PACKAGE_ROOT.rglob("*"):
        if path.is_dir() and path.name in EXCLUDED_DIR_NAMES:
            errors.append(f"excluded directory present: {rel(path)}")
    for path in files:
        if path.suffix in EXCLUDED_SUFFIXES:
            errors.append(f"excluded file type present: {rel(path)}")
        if path.suffix not in TEXT_EXTENSIONS and path.name not in ALLOWED_NO_EXTENSION:
            errors.append(f"unexpected file type: {rel(path)}")
        if (path.suffix in TEXT_EXTENSIONS or path.name in ALLOWED_NO_EXTENSION) and path.resolve() != Path(__file__).resolve():
            text = read_text(path)
            for pattern in PATH_PATTERNS:
                if pattern.search(text):
                    errors.append(f"machine-specific path pattern in {rel(path)}")
                    break
            for pattern in SECRET_PATTERNS:
                if pattern.search(text):
                    errors.append(f"credential/identifier pattern in {rel(path)}")
                    break
        if path.suffix == ".py" and path.name != Path(__file__).name:
            try:
                ast.parse(read_text(path), filename=str(path))
            except SyntaxError as exc:
                errors.append(f"Python syntax error in {rel(path)}: {exc.msg} line {exc.lineno}")

    required = [
        PACKAGE_ROOT / "README.md",
        PACKAGE_ROOT / "config" / "paths.example.yml",
        PACKAGE_ROOT / "config" / "input_layout_v1.tsv",
        PACKAGE_ROOT / "config" / "advanced_ecology_status_gate_v2.tsv",
        PACKAGE_ROOT / "scripts" / "advanced_current" / "run_rf_option_c_v1.R",
        PACKAGE_ROOT / "scripts" / "advanced_current" / "run_residual_network_interaction_v1.R",
        PACKAGE_ROOT / "scripts" / "runnable" / "advanced" / "advanced_ecology_v1.R",
    ]
    for path in required:
        if not path.is_file():
            errors.append(f"required release file missing: {rel(path)}")

    readme = read_text(PACKAGE_ROOT / "README.md")
    if "/your/data/" not in readme or "/outfile/" not in readme:
        errors.append("paper README must document /your/data/ and /outfile/")
    status_text = read_text(PACKAGE_ROOT / "config" / "advanced_ecology_status_gate_v2.tsv")
    if not all(term in status_text for term in ("NCM", "tNST", "not_informative")):
        errors.append("NCM/tNST status gate is incomplete")
    path_template = read_text(PACKAGE_ROOT / "config" / "paths.example.yml")
    if "pr2_dada2_database" in path_template.lower():
        errors.append("DADA2 database path must not be advertised")

    if not (REPOSITORY_ROOT / "LICENSE").is_file():
        errors.append("repository MIT LICENSE is missing")
    if not (REPOSITORY_ROOT / "README.md").is_file():
        errors.append("repository README is missing")
    check_manifest(files, errors)

    try:
        import shutil
        if shutil.which("Rscript"):
            warnings.append("R syntax is available for an optional external parse check")
        else:
            warnings.append("Rscript not found; R syntax was not parsed in this environment")
    except Exception as exc:  # pragma: no cover
        warnings.append(f"R syntax check unavailable: {exc}")

    print("PUBLIC_RELEASE_AUDIT=" + ("PASS" if not errors else "FAIL"))
    print(f"FILES={len(files)}")
    for warning in warnings:
        print("WARNING=" + warning)
    for error in errors:
        print("ERROR=" + error)
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
