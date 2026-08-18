from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path.cwd()
PKG = ROOT / "revision/re_analysis/54_revision_niche_breadth_common75_v1"


def levins_ba(values: pd.Series) -> float:
    total = float(values.sum())
    if total <= 0 or len(values) <= 1:
        return math.nan
    b = total * total / float((values * values).sum())
    return (b - 1.0) / (len(values) - 1.0)


summary = pd.read_csv(PKG / "results/distributional_breadth_taxon_summary_v1.tsv", sep="\t")
cells = pd.read_csv(PKG / "source_data/distributional_breadth_cell_means_v1.tsv", sep="\t")
loso = pd.read_csv(PKG / "results/distributional_breadth_leave_one_site_out_v1.tsv", sep="\t")
stability = pd.read_csv(PKG / "results/distributional_breadth_marker_stability_v1.tsv", sep="\t")
examples = pd.read_csv(PKG / "results/distributional_breadth_descriptive_examples_v1.tsv", sep="\t")

checks: dict[str, dict[str, object]] = {}

counts = summary.groupby("Marker").size().to_dict()
checks["marker_taxon_counts"] = {
    "pass": counts == {"16S": 1034, "18S": 146},
    "observed": counts,
}

cell_counts = cells.groupby(["Marker", "Taxon"]).size()
checks["exact_25_cells_per_taxon"] = {
    "pass": bool((cell_counts == 25).all()),
    "minimum": int(cell_counts.min()),
    "maximum": int(cell_counts.max()),
}

calc = (
    cells.groupby(["Marker", "Taxon"])["CellMeanSampleWiseRelativeAbundance"]
    .apply(levins_ba)
    .rename("RecomputedBA")
    .reset_index()
)
joined = summary.merge(calc, on=["Marker", "Taxon"], how="left", validate="one_to_one")
max_ba_error = float((joined["StandardizedLevinsBA"] - joined["RecomputedBA"]).abs().max())
checks["full_ba_recomputed"] = {"pass": max_ba_error < 1e-12, "max_abs_error": max_ba_error}

recalc_loso = []
for (marker, taxon, omitted), part in cells.merge(
    loso[["Marker", "Taxon", "OmittedSite"]].drop_duplicates(),
    on=["Marker", "Taxon"], how="inner"
).groupby(["Marker", "Taxon", "OmittedSite"]):
    keep = part.loc[part["Site"] != omitted, "CellMeanSampleWiseRelativeAbundance"]
    recalc_loso.append((marker, taxon, omitted, levins_ba(keep)))
recalc_loso = pd.DataFrame(
    recalc_loso, columns=["Marker", "Taxon", "OmittedSite", "RecomputedLOSOBA"]
)
loso_join = loso.merge(recalc_loso, on=["Marker", "Taxon", "OmittedSite"], validate="one_to_one")
max_loso_error = float((loso_join["StandardizedLevinsBA"] - loso_join["RecomputedLOSOBA"]).abs().max())
checks["loso_ba_recomputed"] = {"pass": max_loso_error < 1e-12, "max_abs_error": max_loso_error}

checks["filter_contract"] = {
    "pass": bool(((summary["SamplePrevalence"] >= 0.10) & (summary["SitesDetected"] >= 3)).all()),
    "minimum_prevalence": float(summary["SamplePrevalence"].min()),
    "minimum_sites": int(summary["SitesDetected"].min()),
}

forbidden = ["generalist", "specialist", "ocean", "tidal", "island"]
header_text = " ".join(summary.columns).lower()
checks["no_forced_or_geological_axis"] = {
    "pass": not any(term in header_text for term in forbidden),
    "searched_terms": forbidden,
}

checks["support_gate"] = {
    "pass": bool(stability["SupportGate"].all()),
    "a75_a60_spearman": dict(zip(stability["Marker"], stability["A75A60SpearmanAll"])),
    "minimum_loso_spearman": dict(zip(stability["Marker"], stability["MinimumLOSOSpearmanAll"])),
}
checks["a60_domain_sensitivity"] = {
    "pass": bool(
        (stability["MedianAbsoluteA60A75Delta"] <= 0.10).all()
        and (stability["A60EvaluableTaxa"] > 0).all()
    ),
    "evaluable_taxa": dict(zip(stability["Marker"], stability["A60EvaluableTaxa"])),
    "median_absolute_delta": dict(zip(stability["Marker"], stability["MedianAbsoluteA60A75Delta"])),
    "contract": "A60 detected in >=6/60 samples and >=3/20 cells; excludes both Ocean Sites",
}

valid_boot = summary["SiteStratifiedPositionBootstrapValid"]
checks["position_bootstrap_completeness"] = {
    "pass": bool((valid_boot >= 1999).all()),
    "all_2000": int((valid_boot == 2000).sum()),
    "one_or_more_invalid": int((valid_boot < 2000).sum()),
    "minimum_valid": int(valid_boot.min()),
    "note": "One 16S genus has 1999/2000 finite estimates because one resample was all-zero; its interval uses the 1999 finite estimates.",
}

checks["named_example_gate"] = {
    "pass": bool(
        examples["EligibleAsDescriptiveTailExample"].all()
        and examples["A60SameBreadthDecile"].all()
        and (examples["LOSOSameBreadthDecileCount"] >= 6).all()
        and (~examples["OppositeHalfReversal"]).all()
        and examples.groupby(["Marker", "FullDataBreadthDecile"]).size().between(4, 5).all()
    ),
    "rows": int(len(examples)),
    "note": "Up to five examples per marker and tail were requested; only prevalence>=20%-passing taxa are displayed.",
    "groups": {
        "|".join(key): int(value)
        for key, value in examples.groupby(["Marker", "FullDataBreadthDecile"]).size().items()
    },
}

checks["ba_bounds"] = {
    "pass": bool(summary["StandardizedLevinsBA"].between(0, 1).all()),
    "minimum": float(summary["StandardizedLevinsBA"].min()),
    "maximum": float(summary["StandardizedLevinsBA"].max()),
}

all_pass = all(bool(item["pass"]) for item in checks.values())
report = {
    "status": "PASS" if all_pass else "FAIL",
    "analysis": "common-75 continuous distributional breadth",
    "checks": checks,
}

out = PKG / "qc/distributional_breadth_independent_validation_v1.json"
out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
if not all_pass:
    raise SystemExit("Independent validation failed; inspect " + str(out))
print(json.dumps(report, indent=2, ensure_ascii=False))
