# Coastal 16S/18S analysis scripts v1.1.0

This script-only package accompanies the final revised manuscript. The public project location is:

- https://github.com/sjparkmk/microbial-ecology-analysis/tree/main/papers/2024-coastal-16S-18S

## Scope

The scripts implement the common 75-sample workflow, 16S and PR2-primary 18S processing interfaces, fixed-depth rarefaction (20,000 and 50,000 reads per sample), hierarchical exchangeability analyses, Option-C Site-blocked random forest, interaction-aware residual associations, paired differential abundance, distributional breadth, environmental analyses, cross-marker concordance, and the selected R figure presentations.

RF is reported as Supplementary Table S5 only; no RF plotting script is included. NCM and tNST were removed from the submitted manuscript and supplementary figures, so their diagnostic scripts and figure builders are not included here.

## Paths and external dependencies

Use `scripts/config/paths.example.R` or `scripts/config/paths_example.py` and replace every `/user/link/path` value marked `# insert your path`. Raw FASTQs, reference databases, USEARCH, FAPROTAX, and EasyAmplicon2 helper code are external and are not redistributed. The scripts assume the processed inputs and metadata described in the manuscript and Supplementary Tables S1–S4.

## Reproducibility boundary

This package provides the custom scripts and concise parameter/software ledgers. It does not contain data, figures, tables, binaries, reference databases, or third-party source code. Execution-time versions not recoverable from the original tree-building session are disclosed in `provenance/software_versions_v2.tsv`; no current installed version is substituted for an unknown executed version.
