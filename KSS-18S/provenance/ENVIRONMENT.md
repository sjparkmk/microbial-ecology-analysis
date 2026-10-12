# Recorded environment, not a newly executed fit

The frozen ecological-context RF run recorded Python 3.12.4, NumPy 1.26.4,
pandas 2.2.2, scikit-learn 1.4.2 and Matplotlib 3.8.4. The depth/breadth code
also uses SciPy; the grouped RF module uses NetworkX. The observed project
versions are SciPy 1.13.1 and NetworkX 3.2.1. Numerical execution should use the
recorded versions where possible, then verify against saved source outputs.

Official R analyses were recorded with R 4.6.1, ANCOMBC 2.14.0, ALDEx2 1.44.0
and vegan 2.7-6. `permute` is also required by restricted permutation stages.
These are historical provenance statements, not a claim that a new environment
was installed or all dependencies have been locked for every platform.

The release does not bundle installed packages or a fabricated dependency lock.
Check-only validation needs standard-library Python; scientific execution needs
the packages above and their upstream dependencies. No package is installed by
the launcher. R substeps are called through `Rscript`, selected explicitly with
`--rscript` if it is not on PATH.

Fixed seeds and settings remain inside the original numerical scripts:
20260923 (comparison), 20260924 (ecological-context RF), 20261001 (breadth,
within-KSS RF and depth sensitivity), and 20261007 (A1 restricted sensitivity).
