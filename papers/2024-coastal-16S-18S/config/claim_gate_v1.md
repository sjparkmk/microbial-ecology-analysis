# Prespecified network sensitivity and claim gate

This file records the full diagnostic contract used by the in-house network
module. It is not a claim that every diagnostic edge is displayed in the final
manuscript. The final display is narrower: 54 robust 16S edges and 12 robust
cross-marker edges; the 18S within-marker surface is not displayed.

The frozen displayed-edge definition is unchanged: within each network-wide BH family, `q < 0.01`, `|rho| >= 0.70`, LocationID-refit-bootstrap same-sign fraction `>= 0.95`, and a percentile 95% LocationID-refit-bootstrap interval excluding zero.

The following additional rules are applied before recommending continued Supplementary Figure S11 use:

1. The accepted additive baseline must be reproduced from the protected inputs, including its 57/2/15 strict-edge counts and edge identities under the accepted 199-resample residual-only bootstrap.
2. Full additive and interaction model matrices must be estimable. At least 95% of each 999-resample refit bootstrap must be estimable; otherwise the affected model/network fails closed.
3. An interaction-model network requires 3-60 strict edges and strict-edge density no greater than 0.05 to be considered an interpretable supplementary surface.
4. For retained interaction strict edges, Site-cluster sensitivity and seven leave-one-Site-out refits are diagnostics rather than new discovery gates. An edge is classified as Site/LOSO robust when the Site-refit bootstrap meets the frozen same-sign/interval criteria and at least six of seven estimable leave-one-Site-out correlations retain the full-data sign. At least 80% of interaction strict edges must meet this combined diagnostic for the network-level recommendation `retain_SI_association_only`.
5. Any retained claim remains limited to hypothesis-generating residual association. Direct interaction, hub, bridge, causal, biomarker, and cross-marker network-density comparison claims are prohibited.
6. The 18S within-marker surface remains non-interpretable when fewer than three strict edges are retained, regardless of individual edge magnitude.

These rules are fixed before inspecting the new interaction-aware results.
