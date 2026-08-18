# Coastal 16S/18S microbial ecology analysis code

Versioned public package: `1.0.0-public`

This directory contains the in-house scripts and provenance needed to inspect
or extend the analysis used for the associated coastal amplicon manuscript.
It is deliberately a code-and-contract release: no raw reads, processed ASV
tables, sample metadata, reference databases, generated outputs, or
third-party executables are included.

## What is included

- `scripts/runnable/` — core alpha/beta, taxonomic-composition,
  cross-marker-concordance, environmental-exploratory, distributional-breadth,
  adapter, rarefaction, and advanced-ecology scripts.
- `scripts/advanced_current/` — the current Option-C random-forest and
  interaction-aware residual-network implementations plus their independent
  validator.
- `scripts/figures/` — selected in-house figure-generation scripts. Generated
  figures are not included.
- `scripts/optional_upstream/` — optional upstream orchestration wrappers that
  require author-supplied FASTQs, databases, and tools.
- `scripts/optional_external/` — optional external-dependency lanes. These
  wrappers do not redistribute FAPROTAX or other third-party resources.
- `config/` — analysis parameters, model contracts, path templates, and status
  gates. Copy `paths.example.yml` and use `/your/data/` for inputs and
  `/outfile/` for outputs.
- `provenance/` — software versions, session records, source exclusions, and
  the public file manifest generated for this release.
- `validation/` — read-only public-release audit and smoke checks.

## Final manuscript analysis contract

The public scripts are aligned to the final analysis decisions:

- common design: 75 samples;
- 16S fixed-depth analysis: 20,000 reads per sample;
- 18S fixed-depth analysis: 50,000 reads per sample with PR2-primary
  taxonomy;
- 16S taxonomy: RDP 16S v18/SINTAX; SILVA 138.2 is retained only for
  reference-chimera screening;
- DADA2 was not executed and is not claimed by this package;
- random forest uses the site-blocked Option-C contract;
- networks are interaction-aware residual association analyses, not evidence
  of direct interaction or causality;
- the full network scripts retain a diagnostic baseline contract, while the
  final displayed subset is limited to 54 robust 16S edges and 12 robust
  cross-marker edges; no 18S network display is asserted;
- the environmental exploratory surface is limited to Spring Sea-sand (`n=7`)
  and does not impute Water or Autumn observations.

## NCM and tNST status

The NCM/tNST code is included because it was run as part of the in-house
advanced-ecology analysis. Its non-informative status is a result of the
analysis, not a missing workflow. The accompanying status gate records that
these outputs are diagnostic/descriptive only. They must not be converted into
positive neutral-assembly, community-assembly, or group-level assembly claims.

## Suggested execution order

1. Obtain the external inputs and databases under the paths in
   `config/input_layout_v1.tsv`.
2. Run the appropriate `--check` or `--dry-run` mode before any write mode.
3. Build the local input adapters, if needed, with
   `scripts/runnable/advanced/00_build_input_adapters.R`.
4. Run fixed-depth rarefaction and the core analysis scripts.
5. Run Option-C random forest and interaction-aware residual-network modules
   only after the protected input contract is satisfied.
6. Run `scripts/runnable/advanced/advanced_ecology_v1.R` for the gated
   advanced-ecology diagnostics, including NCM/tNST. This legacy diagnostic
   module is not the final RF/network evidence surface.
7. Render selected figures into `/outfile/` and retain a local session log and
   output hashes.

Many wrappers intentionally begin from frozen processed inputs. The complete
raw-to-ASV lineage is not redistributed here; upstream reconstruction requires
author-supplied FASTQs, reference databases, and licensed tools.

## External resources and license boundary

The root MIT License applies to author-owned code and documentation in this
repository. USEARCH 12, FAPROTAX, PR2, SILVA, RDP, R/Bioconductor packages,
and other external resources are not relicensed or redistributed here. No
EasyAmplicon2 implementation is copied into this repository; the included
adapter is an in-house table-format adapter only.

Release DOI/Zenodo/Figshare link: **to be added after archival deposit**.
