# Public-release validation

Run from the paper directory:

```text
python validation/check_public_script_release_v1.py --check
```

The check is read-only. It verifies the absence of excluded data/archive
classes, local absolute paths and common credential patterns, Python syntax,
placeholder path usage, and the SHA-256 manifest.
