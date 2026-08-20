#!/usr/bin/env python3
"""Run the protected FAPROTAX 1.2.7 collapse_table.py under NumPy >=2.

The upstream script imports the removed ``numpy.NaN`` alias.  This wrapper
restores that alias in memory and executes the original script without
copying or modifying any file in the EasyAmplicon installation.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

import numpy


UPSTREAM = Path("/user/link/path/collapse_table.py")  # insert your path


def main() -> None:
    if not UPSTREAM.is_file():
        raise SystemExit(f"Protected upstream FAPROTAX script is missing: {UPSTREAM}")
    # Compatibility only; the upstream source and database remain byte-unchanged.
    if not hasattr(numpy, "NaN"):
        numpy.NaN = numpy.nan  # type: ignore[attr-defined]
    sys.argv[0] = str(UPSTREAM)
    runpy.run_path(str(UPSTREAM), run_name="__main__")


if __name__ == "__main__":
    main()
