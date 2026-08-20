#!/usr/bin/env Rscript

# Metadata-corrected v2 wrapper for the immutable v1 hierarchical analysis.
# The biological computations and exact transformations are inherited from the
# hash-tracked v1 implementation. v2 corrects only model-specific complement-
# symmetry metadata: A60 Freedman-Lane transformations have a 1/32 nominal
# exact grid and no validated 2/32 complement-symmetry floor.

current_script_path <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve v2 runner path", call. = FALSE)
  normalizePath(sub("^--file=", "", hit), winslash = "/", mustWork = TRUE)
}

V1_RUNNER <- normalizePath(file.path(dirname(current_script_path()), "..", "..",
  "01_hierarchical_exchangeability_v1", "scripts", "run_hierarchical_exchangeability_v1.R"),
  winslash = "/", mustWork = TRUE)
v1_lines <- readLines(V1_RUNNER, warn = FALSE, encoding = "UTF-8")
entry <- grep("^tryCatch\\(main\\(\\), error = function", v1_lines)
if (length(entry) != 1L || entry[[1L]] <= 1L) stop("Cannot isolate immutable v1 implementation", call. = FALSE)
eval(parse(text = v1_lines[seq_len(entry[[1L]] - 1L)], keep.source = TRUE), envir = .GlobalEnv)

SCRIPT_VERSION <- "2.0.0"

input_paths_v1 <- input_paths
input_paths <- function(src) {
  p <- input_paths_v1(src)
  c(p, v1_implementation_runner = V1_RUNNER)
}

output_manifest <- function(stage) {
  files <- list.files(stage, recursive = TRUE, full.names = TRUE)
  files <- files[!file.info(files)$isdir]
  files <- files[!grepl("^output_sha256_manifest_v[12]\\.tsv$", basename(files))]
  root <- normalizePath(stage, winslash = "/", mustWork = TRUE)
  rel <- substring(normalizePath(files, winslash = "/", mustWork = TRUE), nchar(root) + 2L)
  ord <- order(rel)
  data.frame(RelativePath = rel[ord], SizeBytes = as.numeric(file.info(files[ord])$size),
             SHA256 = vapply(files[ord], sha256, character(1L)), stringsAsFactors = FALSE)
}

beta_analysis_v1 <- beta_analysis
beta_analysis <- function(inputs) {
  b <- beta_analysis_v1(inputs)
  is_a60 <- b$results$ModelID == "BETA_TYPE_CAMPAIGN"
  b$results$SymmetryBound[is_a60] <- NA_real_
  b$results$ComplementSymmetryValidated <- !is_a60
  b$results$ExactResolutionNote <- ifelse(is_a60,
    "nominal exact grid 1/32; no complement-symmetry floor applies to the Freedman-Lane A60 statistic",
    sprintf("nominal exact grid 1/%d; complement/global-reversal symmetry validates a 2/%d statistic floor",
            b$results$ExactPatterns, b$results$ExactPatterns))

  b$validation$ComplementSymmetryValidated <- b$validation$ModelID != "BETA_TYPE_CAMPAIGN"
  b$validation$SymmetryMetadata <- ifelse(b$validation$ComplementSymmetryValidated,
    "complement/global-reversal PERMANOVA F symmetry validated across marker and normalization surfaces",
    "A60 Freedman-Lane complement/global-reversal PERMANOVA F symmetry is not valid; no 2/N floor")

  # Verify the symmetry classification directly from every saved exact F surface.
  b$validation$ComplementSymmetryMaxFDifference <- NA_real_
  for (model_id in MODELS) {
    z <- b$pattern_stats[b$pattern_stats$ModelID == model_id, , drop = FALSE]
    diffs <- numeric()
    for (marker in MARKERS) for (surface in SURFACES) {
      zz <- z[z$Marker == marker & z$Surface == surface, , drop = FALSE]
      zz <- zz[order(zz$PatternID), , drop = FALSE]
      diffs <- c(diffs, abs(zz$PERMANOVA_F - rev(zz$PERMANOVA_F)))
    }
    maxdiff <- max(diffs)
    b$validation$ComplementSymmetryMaxFDifference[b$validation$ModelID == model_id] <- maxdiff
    if (model_id == "BETA_TYPE_CAMPAIGN") {
      assert(maxdiff > 1e-8, "A60 complement asymmetry was not reproduced")
    } else {
      assert(maxdiff <= 1e-8, sprintf("Complement symmetry failed for %s", model_id))
    }
  }

  b$gate$MandatoryQualifier <- ifelse(b$gate$ModelID == "BETA_TYPE_CAMPAIGN",
    "Stricter Site-synchronized exact sensitivity; A60 has a 1/32 nominal grid and no complement-symmetry floor; matched four-cell omnibus PERMDISP uses identical transformation IDs",
    "Stricter Site-synchronized exact sensitivity; complement symmetry bounds the F-tail floor at 2/N; matched PERMDISP uses identical transformation IDs")
  b
}

readme_lines <- function() c(
  "# Hierarchical exchangeability sensitivity v2", "",
  "This package supersedes `01_hierarchical_exchangeability_v1` for all use. Biological statistics are unchanged; v2 corrects the model-specific exact-resolution metadata after an independent complement-symmetry audit.", "",
  "## Alpha", "",
  "- A72 uses all 128 one-sign-per-Site patterns while retaining local Sea-sand/Water paired contrasts.",
  "- Both the original 12 Site x Season block weighting and an equal-Site estimand are reported.",
  "- Each strict A72 P value replaces, rather than augments, the A72 row in its original five-contrast BH family.", "",
  "## Beta", "",
  "- A72, A60, A36, and A30 use 128, 32, 64, and 32 unique Site-synchronized transformations, respectively.",
  "- PERMANOVA uses a Freedman-Lane reduced-model residual permutation in the Gower-centered Bray-Curtis kernel.",
  "- A72, A36, and A30 have validated complement/global-reversal F symmetry; their attainable F-tail floor is 2/N.",
  "- A60 does not have complement/global-reversal symmetry under its Freedman-Lane statistic. Its nominal exact resolution is 1/32 and no 2/32 symmetry floor is claimed.",
  "- The A60 transformation reverses the Autumn substrate contrast coherently for all fixed locations within a selected Site. Raw permutation is invalid because it does not preserve the reduced main-effect space, so residual-kernel permutation is mandatory.",
  "- Matched PERMDISP uses the same saved transformation IDs. A60 PERMDISP is a four-cell SampleType x Season omnibus audit, not an interaction-specific dispersion test.",
  "- A42 Site is fixed-site/descriptive only because Site is both the target and highest exchangeability cluster.", "",
  "## Interpretation", "",
  "When accepted pair-restricted and stricter Site-synchronized evidence differ, the stricter result governs manuscript claims. PERMANOVA is multivariate and has no scalar signed direction; LOSO tables therefore report F/R2 stability and explicit R2 changes.", "",
  "Run `Rscript scripts/run_hierarchical_exchangeability_v2.R --check` before use."
)

write_stage_v1 <- write_stage
write_stage <- function(stage, inputs, analyses) {
  write_stage_v1(stage, inputs, analyses)
  old_manifest <- file.path(stage, "provenance", "output_sha256_manifest_v1.tsv")
  assert(file.exists(old_manifest), "Expected transient v1 manifest is missing")
  assert(file.remove(old_manifest), "Could not remove transient staging manifest")

  files <- list.files(stage, recursive = TRUE, full.names = TRUE)
  rename_from <- files[grepl("_v1(?=\\.[^.]+$)", basename(files), perl = TRUE)]
  for (from in rename_from) {
    to <- file.path(dirname(from), sub("_v1(?=\\.[^.]+$)", "_v2", basename(from), perl = TRUE))
    assert(!file.exists(to), sprintf("v2 rename collision: %s", to))
    assert(file.rename(from, to), sprintf("Could not version output: %s", from))
  }

  prior <- normalizePath(file.path(package_root(), "..", "01_hierarchical_exchangeability_v1"),
                         winslash = "/", mustWork = TRUE)
  ledger <- data.frame(
    PriorPackage = prior,
    CurrentPackage = normalizePath(package_root(), winslash = "/", mustWork = TRUE),
    PriorStatus = "SUPERSEDED_METADATA_ONLY",
    CurrentStatus = "CURRENT_ACCEPTED_CANDIDATE_PENDING_PARENT_AUDIT",
    Reason = "Independent audit showed that A60 Freedman-Lane complement/global reversal is not symmetric",
    CorrectedField = "BETA_TYPE_CAMPAIGN SymmetryBound and complement-symmetry interpretation",
    PriorValue = "2/32 symmetry floor claimed",
    CurrentValue = "nominal grid 1/32; no complement-symmetry floor",
    BiologicalNumericChange = "none; F, R2, P, q, PERMDISP, alpha, LOSO, and transformation mappings unchanged",
    Replaces = "01_hierarchical_exchangeability_v1 in full for downstream use",
    stringsAsFactors = FALSE)
  write_tsv(ledger, file.path(stage, "provenance", "v1_to_v2_supersession_ledger_v2.tsv"))
  write_tsv(output_manifest(stage), file.path(stage, "provenance", "output_sha256_manifest_v2.tsv"))
}

verify_manifest <- function(root) {
  manifest <- read_tsv(file.path(root, "provenance", "output_sha256_manifest_v2.tsv"))
  assert_schema(manifest, c("RelativePath", "SizeBytes", "SHA256"), "v2 output manifest")
  for (i in seq_len(nrow(manifest))) {
    p <- file.path(root, manifest$RelativePath[[i]])
    assert(file.exists(p), sprintf("Manifest target missing: %s", p))
    assert(as.numeric(file.info(p)$size) == as.numeric(manifest$SizeBytes[[i]]), sprintf("Size mismatch: %s", p))
    assert(sha256(p) == tolower(manifest$SHA256[[i]]), sprintf("SHA256 mismatch: %s", p))
  }
  actual <- list.files(root, recursive = TRUE, full.names = TRUE)
  actual <- actual[!file.info(actual)$isdir]
  rootn <- normalizePath(root, winslash = "/", mustWork = TRUE)
  rel <- substring(normalizePath(actual, winslash = "/", mustWork = TRUE), nchar(rootn) + 2L)
  expected <- c(manifest$RelativePath, "provenance/output_sha256_manifest_v2.tsv",
                "scripts/run_hierarchical_exchangeability_v2.R")
  assert(setequal(rel, expected), sprintf("Manifest coverage mismatch; extra=%s missing=%s",
    paste(setdiff(rel, expected), collapse = "|"), paste(setdiff(expected, rel), collapse = "|")))
}

check_existing <- function(inputs) {
  root <- package_root(); verify_manifest(root)
  inman <- read_tsv(file.path(root, "provenance", "input_sha256_manifest_v2.tsv"))
  for (i in seq_len(nrow(inman))) {
    p <- inman$AbsolutePath[[i]]
    assert(file.exists(p) && sha256(p) == tolower(inman$SHA256[[i]]), sprintf("Input drift: %s", p))
  }
  qc <- read_tsv(file.path(root, "qc", "hierarchical_exchangeability_validation_v2.tsv"))
  assert(all(qc$Status == "PASS"), "Stored QC contains failure")
  fresh <- analyze_all(inputs)
  saved_a <- read_tsv(file.path(root, "results", "alpha_A72_site_synchronized_results_v2.tsv"))
  compare_tables(fresh$alpha$results, saved_a, c("Marker", "Metric", "Estimand"),
                 c("EffectMeanContrast", "PExactSiteSynchronized"), "v2 alpha sensitivity")
  saved_b <- read_tsv(file.path(root, "results", "beta_hierarchical_permanova_v2.tsv"))
  compare_tables(fresh$beta$results, saved_b, c("Marker", "Surface", "ModelID"),
                 c("F", "R2", "PExactSiteSynchronized", "PAdjustedBHFourHierarchicalTargets"), "v2 beta PERMANOVA")
  saved_d <- read_tsv(file.path(root, "results", "beta_hierarchical_permdisp_v2.tsv"))
  compare_tables(fresh$beta$dispersion, saved_d, c("Marker", "Surface", "ModelID"),
                 c("F", "PExactSameTransformations", "PAdjustedBHFourHierarchicalTargets"), "v2 beta PERMDISP")
  a60 <- saved_b$ModelID == "BETA_TYPE_CAMPAIGN"
  assert(all(is.na(saved_b$SymmetryBound[a60])) && all(!saved_b$ComplementSymmetryValidated[a60]),
         "A60 symmetry metadata correction is absent")
  assert(all(abs(saved_b$SymmetryBound[!a60] - 2/saved_b$ExactPatterns[!a60]) <= TOL) &&
           all(saved_b$ComplementSymmetryValidated[!a60]), "Validated model symmetry metadata drift")

  # Prove that all biological numerical outputs are unchanged from v1.
  prior <- file.path(root, "..", "01_hierarchical_exchangeability_v1", "results")
  v1a <- read_tsv(file.path(prior, "alpha_A72_site_synchronized_results_v1.tsv"))
  compare_tables(v1a, saved_a, c("Marker", "Metric", "Estimand"),
                 c("EffectMeanContrast", "PExactSiteSynchronized"), "v1-to-v2 alpha")
  v1b <- read_tsv(file.path(prior, "beta_hierarchical_permanova_v1.tsv"))
  compare_tables(v1b, saved_b, c("Marker", "Surface", "ModelID"),
                 c("F", "R2", "PExactSiteSynchronized", "PAdjustedBHFourHierarchicalTargets"), "v1-to-v2 beta")
  v1d <- read_tsv(file.path(prior, "beta_hierarchical_permdisp_v1.tsv"))
  compare_tables(v1d, saved_d, c("Marker", "Surface", "ModelID"),
                 c("F", "PExactSameTransformations", "PAdjustedBHFourHierarchicalTargets"), "v1-to-v2 PERMDISP")
  cat(sprintf("PASS: %s\n", root))
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L, runLast = FALSE)
})
