# KSS-18S — analysis code and selected data

Minimal numerical code and processed inputs supporting the **Korean
Solid-Substrate 18S (KSS-18S)** study.

## Data included

| Collection | Libraries | Included material |
|---|---:|---|
| KSS-18S | 15 | Unrarefied counts, fixed 20,000-read counts, PR2 taxonomy and metadata |
| Coastal (`d2024`) | 75 | Selected ALL_EUK processed counts, taxonomy, metadata and selection record |
| EukBank | 50 | Selected ALL_EUK processed counts, taxonomy, metadata and accessions from 32 projects |
| Harmonized comparison | 140 | Frozen lineage-qualified Genus/Family counts, metadata and existing Bray–Curtis distances |

The Coastal libraries are independently collected biological samples. The
historical label of 25 sampling combinations does **not** reduce these data to
25 pooled libraries. KSS libraries represent independent locations, not repeated
libraries within each location. The harmonized comparison has 34 source projects.

Only the selected external subsets are included, not the complete EukBank
database. There are no raw FASTQ files, database FASTA files, manuscripts,
EndNote libraries, private conversations, machine configuration or credentials.
Consolidated scientific-result workbooks and a three-page Additional Data
figure PDF are available in [supplements/](supplements/README.md).

```
KSS-18S/
  run.py                    # check-only by default; opt-in stage dispatch
  scripts/                  # 6 scientific Python and 4 scientific R scripts
  input/
    kss/                    # 15 study libraries and associated metadata
    d2024/                  # 75 selected Coastal libraries
    eukbank/                # 50 selected libraries + accession/project records
    comparative/            # frozen 140-library aggregate inputs + distances
    reference_results/      # four compact existing result tables, not rerun
  provenance/               # source hashes, portability scope and attribution
  supplements/              # two indexed result workbooks and three AD figures
```

## Start with read-only validation

Python 3.12 is the recorded analysis interpreter. The default input/output roots
are `/input/ur/data/` and `/output/ur/data/`; both can be overridden. The checked-in
`input/` directory supplies the input root. Inputs and outputs must be separate,
non-nested directories. A stage refuses an existing output directory.

From this directory, using a scientific Python environment:

```sh
python run.py --input-root ./input --output-root ../KSS-18S-output --check
```

The check uses Python's standard library: SHA-256, schema, sample identity/order,
integer counts, taxonomy coverage, cohort sizes and the fixed 20k totals. It does
not fit a model, compute a hypothesis test, draw subsamples or generate figures.
On Windows, supply native paths to the same options. No WSL, downloads or runtime
installation are required or performed by this command.

## Optional execution by a future reuser

**The commands below perform analyses. They were NOT run during this release.**
Use a suitable preconfigured Python/R environment; see
[ENVIRONMENT.md](provenance/ENVIRONMENT.md). There is no automatic installation.

```sh
python run.py --input-root ./input --output-root ../KSS-18S-output --run baseline
python run.py --input-root ./input --output-root ../KSS-18S-output --run composition
python run.py --input-root ./input --output-root ../KSS-18S-output --run community
python run.py --input-root ./input --output-root ../KSS-18S-output --run daa
python run.py --input-root ./input --output-root ../KSS-18S-output --run daa-retry
python run.py --input-root ./input --output-root ../KSS-18S-output --run aldex
python run.py --input-root ./input --output-root ../KSS-18S-output --run rf-grouped
python run.py --input-root ./input --output-root ../KSS-18S-output --run rf-lopo
python run.py --input-root ./input --output-root ../KSS-18S-output --run breadth-depth
python run.py --input-root ./input --output-root ../KSS-18S-output --run depth-statistics
python run.py --input-root ./input --output-root ../KSS-18S-output --run a1
```

For R stages, `--rscript "path/to/Rscript"` selects an existing installation.
`depth-statistics` requires `breadth-depth` output. Other stages use the shipped
frozen inputs independently. Optional `--run prepare` rebuilds harmonized input
tables into a **new output folder**; it does not overwrite or silently substitute
the shipped inputs used by other stages.

The `daa-retry` stage preserves a historical estimator-failure exception: the
exact Agaricomycetes lineage named by the zero-coefficient-variance failure is
excluded only from the KSS Family model. It is not selected by significance.
The original `daa` result and the separate retry remain separate. Historical
ANCOM-BC2 outputs may include failed fits; successful process exit is not proof
that every model converged. Inspect `failure.txt`, `completed.txt`, warnings and
model-contract files. The launcher never promotes a failed model to a finding.

## Scientific scope and interpretation

- KSS diversity and recurrence use the fixed 20k table. DAA uses unrarefied
  counts. ASV identifiers are source-local: they are not merged across studies.
- Cross-study comparisons use full lineage-qualified Family/Genus identities
  and existing bootstrap-supported assignments. Historical placeholder labels
  remain unchanged; no taxonomy or lifestyle reassignment was performed.
- `ours`, `d2024` and `eukbank` remain machine identifiers to preserve joins and
  exact source definitions. Their display names are KSS-18S, Coastal and EukBank.
- Lifestyle annotations express taxonomy-based **potential affinities**, not
  measured activity, biomass or obligate nutritional mode.
- Breadth describes distribution over the sampled strata, not physiological
  tolerance. KSS pH midpoints are existing analysis values. Coastal pH records
  are site-level, not 75 independent pH measurements. EukBank pH is not imputed.
- RF project-LOPO performance concerns transfer among the sampled source
  projects. Full-training importance is exploratory, not marker validation.
- The three retained Family marine-water versus land-sediment terms are
  **between-study candidates**, not source-independent habitat indicators:
  no source study spans both categories. Oxytrichidae and Sandonidae have
  nonsingular fits; Litostomatea_XX retains a singular-fit caveat.
- Historical denoising sensitivity is not recurrence validation by repeated
  subsampling. No new recurrence or network test is asserted by this release.

## What this minimal package does not reproduce

This is not a FASTQ-to-publication pipeline or a complete figure-layout archive. It does
not include DADA2/primer processing, complete reference databases, every legacy
exploratory script, tree inference, or all final figure-rendering dependencies.
The exact original KSS alpha-test/PERMANOVA implementation was not located in
the bounded canonical source inventory; existing study results are preserved
in `input/reference_results/`, and no replacement statistical implementation
was invented. `baseline` reproduces descriptive alpha/recurrence definitions,
not those historical inferential tests. `community` is the separate aggregate
Family/Genus comparison and must not be presented as that missing ASV test.

`composition` is an explicitly documented numeric-only extraction of the
existing composition script; private plotting-style dependencies are omitted.
An inherited diagnostic RF plot may be emitted by `rf-lopo`; it is not the final
Figure 5 layout. Refer to [STAGE_MAP.md](provenance/STAGE_MAP.md) for scope.

Full end-to-end execution and numerical equality of regenerated results have
**not** been tested in this release. Source hashes, limited code-diff checks,
syntax checks and structural input validation are distinct from scientific
reproduction. See [VALIDATION.md](provenance/VALIDATION.md).

## Provenance and reuse

[source_manifest.json](provenance/source_manifest.json) records source and
release hashes without private absolute paths. KSS count-table changes are
strict 15-column projections: no retained count is changed. Both the original
frozen analysis metadata and the later recorded KSS metadata are preserved with
their separate roles. [DATA_SOURCES.md](provenance/DATA_SOURCES.md) gives EukBank
attribution, license and snapshot limits. Author-owned code/documentation use
the repository MIT License; third-party data retain their source terms.
