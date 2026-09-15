from __future__ import annotations

import argparse
import hashlib
import shutil
import zipfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-docx", type=Path, default=Path("/input/youre/data/figures.docx"))
    parser.add_argument("--output-png", type=Path, default=Path("/output/your/data/Figure_1_original_v2.png"))
    parser.add_argument("--log", type=Path, default=Path("/output/your/data/Figure_1_extraction.log"))
    args = parser.parse_args()
    if args.output_png.exists() or args.log.exists():
        raise FileExistsError("Refusing to overwrite an existing Figure 1 output or log")
    with zipfile.ZipFile(args.source_docx) as archive:
        payload = archive.read("word/media/image1.png")
    args.output_png.parent.mkdir(parents=True, exist_ok=True)
    args.output_png.write_bytes(payload)
    source_extract = args.output_png.parent / "source_data" / f"{args.output_png.stem}_embedded_image1.png"
    if source_extract.exists():
        raise FileExistsError(source_extract)
    source_extract.parent.mkdir(parents=True, exist_ok=True)
    source_extract.write_bytes(payload)
    if sha256(args.output_png) != sha256(source_extract):
        raise RuntimeError("Extracted Figure 1 copies do not match")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text(
        "Figure 1 original-image extraction\n"
        f"source_docx={args.source_docx}\n"
        "source_member=word/media/image1.png\n"
        f"output_png={args.output_png}\n"
        f"source_extract={source_extract}\n"
        f"sha256={sha256(args.output_png)}\n"
        f"bytes={args.output_png.stat().st_size}\n",
        encoding="utf-8",
    )
    print(f"WROTE {args.output_png}")


if __name__ == "__main__":
    main()
