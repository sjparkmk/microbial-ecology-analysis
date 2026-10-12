# Release validation — 12 October 2026

## Completed without scientific reanalysis

- 39 release source records checked against original SHA-256 values. Source
  data and scientific code remain unchanged outside this release copy.
- 29 input files pass hash checks. Study sizes: KSS 15, Coastal 75, selected
  EukBank 50; harmonized inputs 140 libraries and 34 source projects.
- KSS raw/fixed-20k tables retain 5,136 original ASV rows; Coastal and EukBank
  ALL_EUK tables have 24,302 and 39,904 ASV rows, respectively. These are table
  row counts, not claims that all rows occur in every subset/library.
- KSS 15-column count and 15-row metadata projections were compared cell-for-
  cell with their source subsets. No retained value changed. Frozen 20k totals
  are unchanged for all 15 libraries.
- Sample uniqueness, count-column membership/order, nonnegative integer counts,
  positive column totals and taxonomy coverage pass structural checks.
- All 9 Python files parse/compile without executing their analysis entry points.
  All 4 R files parse under R 4.6.1 without evaluating statistical expressions.
- Six synthetic launcher tests pass: default check-only behavior, path traversal
  rejection, nested input/output rejection, existing-output refusal, missing
  depth-stage dependency refusal and correct retry dispatch.
- RF numerical functions are AST-identical to their source. Composition
  annotation/aggregation functions are AST-identical; the retained numerical
  main prefix is identical after normalizing two descriptive column names.
  Other unchanged numerical helper functions were compared; changed
  path/guard/documentation statements were reviewed separately.
- No prohibited binary/raw/manuscript files or inspected private path/credential
  patterns were found in the package. Selected EukBank IDs reconcile 50/50,
  with exactly 32 projects. Third-party attribution is documented separately.
- All nine protected manuscript/SI/legend deliverable hashes match their
  existing file-specific manifest. None is included in this public release.

## Limits and checks not run

No scientific stage was executed: no statistical test, model refit, bootstrap,
permutation rerun, rarefaction draw, taxonomic reassignment or figure rendering.
Consequently these checks do not establish end-to-end numerical reproduction.
The missing historical KSS primary inferential source and excluded layout/audit
steps are explicitly described in the README and stage map.

The Python scientific packages are importable at the versions in ENVIRONMENT.md.
The default R session used for syntax checking did not expose vegan, permute,
ANCOMBC or ALDEx2; it also emitted locale startup warnings. No installation or
environment change was attempted. Historical recorded versions are not a claim
of current default-R readiness. Configure an existing compatible R library path
before any optional future R execution.

The first launcher-test attempt could not clean up temporary directories in the
Windows sandbox's default temp area. Running the same tests with a writable,
task-local temporary root passed; the test code was unchanged. An initial local
preservation check mistakenly looked for DOCX hashes among Markdown files. It
was corrected to use the existing file-specific nine-deliverable manifest;
there was no manuscript modification. Input-check feedback also identified the
excluded context row in original KSS metadata; release metadata now uses the
same exact 15-study selection as the count tables.

Run the portable read-only check with:

```sh
python run.py --input-root ./input --output-root ../KSS-18S-output --check
python -B tests/test_launcher.py
```
