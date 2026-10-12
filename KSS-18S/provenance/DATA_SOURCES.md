# Input provenance and attribution

## KSS-18S and Coastal

The study authors authorized this release of the processed 15-library KSS-18S
and 75-library Coastal inputs. Existing source sample IDs and metadata values
are retained. Contact the repository owner regarding reuse of author-provided
data; this release does not assign a new data license on the authors' behalf.
The MIT License applies to author-owned code and documentation.

KSS raw and fixed-20k tables originally included one extra water-column context
library. This release projects the tables to the 15 study libraries selected by
the existing metadata. All retained values and all original ASV rows are kept;
zero-only rows are not silently discarded. This mechanical projection is not
reanalysis. Original and projected checksums are both recorded.

The original `kss/inputs/metadata_analysis.tsv` supports the frozen comparison;
`kss/metadata_analysis_current.tsv` is the later recorded metadata used by the
existing breadth/geography script. They are not silently interchanged. The
author-confirmed KSS pH procedure uses distilled water at 1:5 soil:water (w/v).
Existing pH point/range records and their analysis midpoints are unchanged.

## Selected EukBank subset

Resource: Berney, C.; Mahé, F.; Henry, N.; Lara, E.; de Vargas, C.; and the
EukBank consortium. *EukBank 18S V4 dataset*. Zenodo,
[doi:10.5281/zenodo.7804946](https://doi.org/10.5281/zenodo.7804946).

The official record was checked on 12 October 2026 via
`https://zenodo.org/api/records/7804946`: open access, license `cc-by-4.0`.
The resource's [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) terms
remain applicable to redistributed EukBank-derived material. Retain attribution,
license information and an indication of selection/processing in downstream use.
The repository's code license does not supersede these data terms.

The release contains a **selected, processed subset**, not the upstream whole
database: 50 libraries from 32 source projects, with their available accessions,
sample metadata and project memberships. Count and taxonomy files are the
existing ALL_EUK subsets used by the analysis. Existing selection records are
included; no new selection or processing was done during release assembly.

The recorded selection used V4-primer eligibility, a read-depth screen, available
coordinates and five operational sample categories. It is purposive, not a
random representative sample of all EukBank projects. The exact original local
download timestamp/version label has not been independently established. The
DOI identifies the resource; the released-file SHA-256 hashes identify the exact
local analysis subset. Do not call this a newly downloaded or verified latest
database snapshot.

The official API reports original-resource MD5 values for reference only:

| Upstream file | MD5 from official record |
|---|---|
| eukbank_18S_V4_counts.tsv.gz | c0f4053cee2146db0ab3566c82b1b560 |
| eukbank_18S_V4_asvs.tsv.gz | f81b4918edb9666f116b450be582d201 |
| eukbank_18S_V4_samples.tsv.gz | cf94771224d7c5d61cd1089f536ee021 |

These are **not** asserted to be checksums of the subset files here. The full
upstream files were not downloaded to assemble this release.
