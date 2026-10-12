# Existing numerical contract (release summary)

No new analysis is authorized or reported by assembling this package. This
public summary replaces private absolute provenance paths, not scientific rules.

The breadth/RF/depth implementation retains seed 20261001 and 100 subsampling
draws at both 14,000 and 10,000 reads. It first subsamples existing unrarefied
Eukaryota ASV counts and then aggregates to Genus/Family. Mean distance matrices
represent each library once. The primary 15-library table stays at 20,000 reads.

Class breadth uses bootstrap-supported full lineages and all-eukaryote read
denominators. It averages within observed strata before normalizing a taxon's
stratum means and applying standardized Levins breadth. Taxa must occur in at
least two libraries. Existing pH midpoints and 7 Coastal site records support
descriptive abundance-weighted pH profiles, not inferred tolerances or causal
pH effects. EukBank pH is not imputed. KSS exposure/geography RF is leave-one-
library-out within one study, with training-fold filtering and a dummy baseline;
it is distinct from the 34-project ecological-context LOPO analysis.

Depth statistics use existing vegan implementations and 9,999 within-dataset
permutations. The later A1 stage retains finer conditioning and restricted
permutations, seed 20261007. Original formulas, factor levels, filters and model
settings remain in the provided code. The source manifest identifies that code.

External associations remain
source/design bounded; between-study Family terms are conditional candidates,
not validated ecological indicators. No project-independent causal abiotic
effect is established by these scripts.
