# KSS-18S consolidated supplementary results

These files consolidate existing study results for browsing and download. No
analysis, model, resampling procedure or statistical test was rerun to create
them. The main manuscript and Supporting Information manuscript are not included.

| File | Contents |
|---|---|
| [Supplementary_Tables_S1-S18.xlsx](Supplementary_Tables_S1-S18.xlsx) | Tables S1–S18; 65 component worksheets plus Contents |
| [Additional_Data_Tables_AD1-AD19.xlsx](Additional_Data_Tables_AD1-AD19.xlsx) | Tables AD1–AD19; 29 component worksheets and 11 summary worksheets plus Contents |
| [Additional_Data_Figures_AD1-AD3.pdf](Additional_Data_Figures_AD1-AD3.pdf) | Three figures, one page each, in numerical order with bookmarks |

Start at **Contents** in either workbook. Table IDs and original scientific row
and column order are retained. TSV-derived sheets have headers at row 5 and data from row 6;
summary sheets retain their source layout. Excel has a 15-significant-digit
numeric guarantee: longer numeric tokens are intentionally stored as text, as
are sample/feature identifiers and missing-value tokens. They are not zeros.
Existing summary-table display precision is retained. These are static result
tables, not live formulas or newly fitted models.

Public-release processing is limited to administrative metadata: local absolute
paths become source-file identifiers; internal editorial labels are removed;
the primary KSS count-table label is neutralized. All numerical table cells,
sample IDs, coordinates, accessions, statistical results and scientific
qualifications are retained. The first figure's footer uses “KSS baseline”.
No plotted value, geometry or statistical conclusion is changed. Historical
source hashes in provenance tables identify original records, not the current
workbooks. Current downloadable-file checksums are in
[SHA256SUMS](SHA256SUMS).

Original operational scripts, private execution logs and manuscript files are
not part of this public supplement folder. Historical record references in AD1
and AD8 do not imply that those private files are downloadable here. Portable
public code remains in [../scripts/](../scripts/) with scope and limitations in
[../provenance/](../provenance/).

## Additional Data figure descriptions

- **AD1 — Denoising sensitivity:** Historical independent versus pseudo-pooled
  DADA2 comparison at 20,000 reads. This is sensitivity to denoising strategy,
  not PCR/field replication or repeated-subsampling validation of recurrence.
- **AD2 — Breadth/null diagnostics:** The fixed-margin null has an incompatible
  zero structure; these diagnostics do not establish ecological generalist or
  specialist classes.
- **AD3 — Coordinate coverage:** Source-reported locations for 15 KSS and 50
  selected EukBank libraries and missing-coordinate coverage for 75 Coastal
  libraries. No coordinate imputation/geocoding and no geographic RF or
  environmental-association test are implied by this figure.

The 15 KSS libraries remain the primary study. External comparisons are
contextual; Family candidates remain between-study associations, not validated
source-independent habitat indicators. See the full scientific scope in
[../README.md](../README.md) and third-party attribution/terms in
[../provenance/DATA_SOURCES.md](../provenance/DATA_SOURCES.md).
