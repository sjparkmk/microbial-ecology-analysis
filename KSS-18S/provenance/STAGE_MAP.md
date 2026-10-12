# Stage-to-evidence map

| Stage | Existing source lineage | Numerical scope / V5.13 association |
|---|---|---|
| prepare | Run A finalization | Frozen Genus/Family input preparation; Figure 3–5 input basis |
| baseline | Run A baseline support | Individual alpha, ASV recurrence and abundance thresholds; Tables S4/S9 |
| composition | Run A follow-up composition | Exact existing trophic annotation and aggregate counts; Table S11; no plotting |
| community | Run A comparison R | Aggregate Bray–Curtis/PCoA/PERMANOVA/PERMDISP; Table S8; historical diagnostics remain labelled |
| daa / daa-retry | Run A comparison R + documented failure retry | ANCOM-BC2 numerical models; finite effects and structural zeros separate; Table S14 basis |
| aldex | Run A comparison R | ALDEx2 two-group sensitivity, including KSS ASV exposure; Table S14 basis |
| rf-grouped | Run A comparison Python | Cluster-grouped RF and diagnostics; Table S17 |
| rf-lopo | 24 September Run A follow-up | Source-project LOPO and fixed-OOF project bootstrap; Table S16 |
| breadth-depth | 1 October additions | Descriptive breadth, KSS RF, fixed depth-sensitivity draws; Tables S2/S13/S18 basis |
| depth-statistics | 1 October R substep | Marginal PERMANOVA/dispersion on six precomputed depth distances; Table S13 |
| a1 | 7 October restricted sensitivity | Fine-block-conditioned sensitivity; current external-comparison interpretation |

These are numerical/source associations, not a claim that each stage alone
regenerates the final edited table or figure. Final layout/export and several
reporting audits are intentionally excluded from this minimal release.

The compact frozen `reference_results/merged_planned_clustered_candidate_support.tsv`
preserves existing candidate interpretation evidence; the mixed-model audit
refitting implementation is not included in this minimal launcher. Its results
must not be relabelled as a new analysis. KSS exact historical alpha/PERMANOVA
tables are also preserved as reference results because their primary source
implementation was not located within the bounded inventory.

Portability changes are limited to roots/import discovery, explicit execution
guards, obsolete generated-document labels, and removal of plotting-only code
from the composition extraction. No statistical formula, model hyperparameter,
threshold, seed, count value or annotation rule is deliberately changed.
