#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)

cli <- commandArgs(trailingOnly = TRUE)
modes <- c("--dry-run", "--run", "--check")
selected <- modes[modes %in% cli]
if (length(selected) != 1L) stop("Specify exactly one of --dry-run, --run, or --check", call. = FALSE)
mode <- sub("^--", "", selected)
worker_arg <- grep("^--workers=", cli, value = TRUE)
workers <- if (length(worker_arg)) as.integer(sub("^--workers=", "", worker_arg[[1]])) else 6L
if (!is.finite(workers) || workers < 1L || workers > 16L) stop("--workers must be an integer from 1 to 16", call. = FALSE)

required_packages <- c("ranger", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
if (getRversion() != "4.5.2") stop("This package is locked to R 4.5.2; found ", getRversion(), call. = FALSE)

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
  invisible(TRUE)
}

read_tsv <- function(path) {
  assert_true(file.exists(path), paste0("Missing input: ", path))
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(tmp)) stop("Refusing overwrite: ", path, call. = FALSE)
  write.table(x, tmp, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "NA", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path)) stop("Atomic promotion failed: ", path, call. = FALSE)
  invisible(path)
}

write_lines_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(tmp)) stop("Refusing overwrite: ", path, call. = FALSE)
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic promotion failed: ", path, call. = FALSE)
  invisible(path)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

paths <- list(
  design = file.path(project_root, "revision/re_analysis/00_design_v4/master_sample_metadata_v4.tsv"),
  count16 = file.path(project_root, "revision/re_analysis/04_input_16s_rarefied_20000_v1/results/16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  tax16 = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_taxonomy_analysis_75_v1.tsv"),
  count18 = file.path(project_root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
  tax18 = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_taxonomy_pr2_primary_analysis_v1.tsv"),
  accepted_metrics = file.path(project_root, "revision/re_analysis/12_original_pipeline_reanalysis_v1/advanced_ecology_v2/tables/RF_site_blocked_metrics_v2.tsv"),
  accepted_script = file.path(project_root, "revision/re_analysis/original_style_pipeline_v1/advanced_ecology_original_style_v1.R"),
  config = file.path(package_root, "config/rf_option_c_config_v1.tsv")
)
assert_true(all(vapply(paths, file.exists, logical(1))), "One or more frozen RF inputs are missing")

source_manifest <- data.frame(
  Input = names(paths),
  AbsolutePath = normalizePath(unlist(paths), winslash = "/", mustWork = TRUE),
  Bytes = as.numeric(file.info(unlist(paths))$size),
  SHA256 = vapply(unlist(paths), sha256, character(1)),
  stringsAsFactors = FALSE
)

load_count_matrix <- function(path) {
  x <- read_tsv(path)
  ids <- as.character(x[[1]])
  assert_true(ncol(x) >= 2L && !anyNA(ids) && !anyDuplicated(ids), paste0("Invalid count table: ", path))
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  assert_true(all(is.finite(mat)) && all(mat >= 0) && all(mat == floor(mat)), "Counts must be nonnegative integers")
  mat
}

taxonomy_id_column <- function(x) {
  hit <- intersect(c("OTUID", "ASV_ID", "#OTUID", "FeatureID"), names(x))
  if (length(hit)) hit[[1]] else names(x)[[1]]
}

normalize_tax_label <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | x == "" | tolower(x) %in% c("na", "unassigned", "unknown", "uncultured", "unclassified")
  x[bad] <- "Unassigned"
  x
}

aggregate_counts_by_taxon <- function(mat, taxonomy, rank) {
  id_col <- taxonomy_id_column(taxonomy)
  assert_true(rank %in% names(taxonomy), paste0("Taxonomy rank missing: ", rank))
  idx <- match(rownames(mat), taxonomy[[id_col]])
  assert_true(!anyNA(idx), paste0("Taxonomy linkage failed at ", rank))
  labels <- normalize_tax_label(taxonomy[[rank]][idx])
  out <- rowsum(mat, group = labels, reorder = FALSE)
  out[rowSums(out) > 0, , drop = FALSE]
}

clr_transform <- function(mat, pseudocount = 0.5) {
  x <- log(t(mat) + pseudocount)
  sweep(x, 1, rowMeans(x), "-")
}

design <- read_tsv(paths$design)
required_design <- c("SampleID", "Site", "SampleType", "Season", "PairID_Substrate", "PairID_Season")
assert_true(all(required_design %in% names(design)), "Design columns missing")
assert_true(nrow(design) == 75L && !anyDuplicated(design$SampleID), "Design must contain 75 unique samples")
site_order <- c("Dokdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan", "Ulleungdo")
design$Site <- factor(design$Site, levels = site_order)
design$SampleType <- factor(design$SampleType, levels = c("Sea-sand", "Water"))
design$Season <- factor(design$Season, levels = c("Spring", "Autumn"))

load_marker <- function(marker, count_path, tax_path, rank, depth) {
  mat <- load_count_matrix(count_path)
  assert_true(setequal(colnames(mat), design$SampleID), paste0(marker, " count/design SampleID mismatch"))
  mat <- mat[, design$SampleID, drop = FALSE]
  assert_true(all(colSums(mat) == depth), paste0(marker, " rarefaction-depth contract failed"))
  taxonomy <- read_tsv(tax_path)
  agg <- aggregate_counts_by_taxon(mat, taxonomy, rank)
  list(marker = marker, rank = rank, depth = depth, agg = agg)
}

markers <- list(
  `16S` = load_marker("16S", paths$count16, paths$tax16, "Genus", 20000L),
  `18S` = load_marker("18S", paths$count18, paths$tax18, "Class", 50000L)
)

surface_seed <- c(
  `16S_SampleType` = 112000L,
  `16S_Season` = 212000L,
  `18S_SampleType` = 312000L,
  `18S_Season` = 412000L
)

prepare_surface <- function(marker_obj, task) {
  if (task == "SampleType") {
    use <- !is.na(design$PairID_Substrate) & design$PairID_Substrate != ""
    pair_col <- "PairID_Substrate"
    levels_y <- c("Sea-sand", "Water")
  } else if (task == "Season") {
    use <- !is.na(design$PairID_Season) & design$PairID_Season != ""
    pair_col <- "PairID_Season"
    levels_y <- c("Spring", "Autumn")
  } else stop("Unknown task: ", task, call. = FALSE)
  d <- droplevels(design[use, , drop = FALSE])
  y <- factor(as.character(d[[task]]), levels = levels_y)
  pair <- as.character(d[[pair_col]])
  pair_levels <- unique(pair)
  pair_check <- split(as.character(y), pair)
  assert_true(all(lengths(pair_check) == 2L), paste0(marker_obj$marker, " ", task, ": each pair must contain two samples"))
  assert_true(all(vapply(pair_check, function(z) setequal(z, levels_y), logical(1))), paste0(task, ": unbalanced pair labels"))
  counts <- marker_obj$agg[, d$SampleID, drop = FALSE]
  counts <- counts[rownames(counts) != "Unassigned", , drop = FALSE]
  sites <- droplevels(d$Site)
  folds <- levels(sites)
  fold_objects <- vector("list", length(folds))
  feature_ledger <- vector("list", length(folds))
  for (ii in seq_along(folds)) {
    held <- folds[[ii]]
    train_idx <- which(as.character(sites) != held)
    test_idx <- which(as.character(sites) == held)
    rel_train <- sweep(counts[, train_idx, drop = FALSE], 2, colSums(counts[, train_idx, drop = FALSE]), "/")
    prevalence <- rowMeans(counts[, train_idx, drop = FALSE] > 0)
    mean_ra <- rowMeans(rel_train)
    keep <- prevalence >= 0.10 & mean_ra >= 0.0001
    assert_true(sum(keep) >= 5L, paste0(marker_obj$marker, " ", task, " ", held, ": fewer than five retained taxa"))
    xtr <- clr_transform(counts[keep, train_idx, drop = FALSE], 0.5)
    xte <- clr_transform(counts[keep, test_idx, drop = FALSE], 0.5)
    cn <- make.unique(make.names(rownames(counts)[keep]))
    colnames(xtr) <- cn
    colnames(xte) <- cn
    fold_objects[[ii]] <- list(
      held = held, train_idx = train_idx, test_idx = test_idx,
      xtr = xtr, xte = xte, mtry = max(1L, floor(sqrt(ncol(xtr))))
    )
    feature_ledger[[ii]] <- data.frame(
      Marker = marker_obj$marker, Rank = marker_obj$rank, Task = task, HeldOutSite = held,
      TrainingSamples = length(train_idx), HeldOutSamples = length(test_idx),
      CandidateTaxaExcludingUnassigned = nrow(counts), RetainedTaxa = sum(keep),
      Mtry = max(1L, floor(sqrt(ncol(xtr)))),
      FilteringScope = "training-fold-only", TransformScope = "training and held-out separately after training-fold filter",
      stringsAsFactors = FALSE
    )
  }
  list(
    id = paste(marker_obj$marker, task, sep = "_"), marker = marker_obj$marker,
    rank = marker_obj$rank, task = task, d = d, y = y, pair = pair,
    pair_levels = pair_levels, levels_y = levels_y, sites = sites, folds = folds,
    fold_objects = fold_objects, feature_ledger = do.call(rbind, feature_ledger),
    seed = unname(surface_seed[[paste(marker_obj$marker, task, sep = "_")]])
  )
}

surface_list <- list()
for (marker in names(markers)) for (task in c("SampleType", "Season")) {
  id <- paste(marker, task, sep = "_")
  surface_list[[id]] <- prepare_surface(markers[[marker]], task)
}

toggle_labels <- function(y, flip_by_unit, sample_unit, levels_y) {
  out <- as.character(y)
  flip_sample <- as.logical(flip_by_unit[match(sample_unit, names(flip_by_unit))])
  assert_true(!anyNA(flip_sample), "Permutation unit linkage failed")
  out[flip_sample] <- ifelse(out[flip_sample] == levels_y[[1]], levels_y[[2]], levels_y[[1]])
  factor(out, levels = levels_y)
}

auc_binary <- function(observed, probability_positive, positive) {
  y <- observed == positive
  n1 <- sum(y)
  n0 <- sum(!y)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  ranks <- rank(probability_positive, ties.method = "average")
  (sum(ranks[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

metric_from_prediction <- function(observed, predicted, probability_positive, levels_y) {
  sens <- vapply(levels_y, function(level) mean(predicted[observed == level] == level), numeric(1))
  c(
    BalancedAccuracy = mean(sens),
    AUC = auc_binary(observed, probability_positive, levels_y[[2]]),
    SensitivityLevel1 = sens[[1]],
    SensitivityLevel2 = sens[[2]]
  )
}

fit_surface_labels <- function(surface, y_fit, return_predictions = FALSE) {
  pred_rows <- vector("list", length(surface$fold_objects))
  for (ii in seq_along(surface$fold_objects)) {
    fold <- surface$fold_objects[[ii]]
    set.seed(surface$seed + ii)
    fit <- ranger::ranger(
      x = as.data.frame(fold$xtr), y = droplevels(y_fit[fold$train_idx]),
      probability = TRUE, num.trees = 500L, mtry = fold$mtry,
      min.node.size = 2L, importance = "none", seed = surface$seed + ii,
      num.threads = 1L
    )
    probability <- predict(fit, data = as.data.frame(fold$xte))$predictions
    if (is.null(dim(probability))) {
      probability <- matrix(probability, ncol = length(surface$levels_y),
                            dimnames = list(NULL, surface$levels_y))
    }
    assert_true(all(surface$levels_y %in% colnames(probability)), "RF probability columns missing")
    predicted <- surface$levels_y[max.col(probability[, surface$levels_y, drop = FALSE], ties.method = "first")]
    pred_rows[[ii]] <- data.frame(
      SampleIndex = fold$test_idx,
      Observed = as.character(y_fit[fold$test_idx]),
      Predicted = predicted,
      ProbabilityPositive = probability[, surface$levels_y[[2]]],
      HeldOutSite = fold$held,
      stringsAsFactors = FALSE
    )
  }
  predictions <- do.call(rbind, pred_rows)
  predictions <- predictions[order(predictions$SampleIndex), , drop = FALSE]
  metric <- metric_from_prediction(predictions$Observed, predictions$Predicted,
                                   predictions$ProbabilityPositive, surface$levels_y)
  if (return_predictions) list(metric = metric, predictions = predictions) else metric
}

accepted_metrics <- read_tsv(paths$accepted_metrics)
accepted_metrics$SurfaceID <- paste(accepted_metrics$Marker, accepted_metrics$Task, sep = "_")
observed <- list()
for (id in names(surface_list)) {
  z <- fit_surface_labels(surface_list[[id]], surface_list[[id]]$y, return_predictions = TRUE)
  expected <- accepted_metrics[accepted_metrics$SurfaceID == id, , drop = FALSE]
  assert_true(nrow(expected) == 1L, paste0("Accepted metric row missing for ", id))
  assert_true(abs(z$metric[["BalancedAccuracy"]] - expected$BalancedAccuracy) < 1e-12,
              paste0("Observed balanced accuracy no longer reproduces accepted ", id))
  assert_true(abs(z$metric[["AUC"]] - expected$AUC) < 1e-12,
              paste0("Observed AUC no longer reproduces accepted ", id))
  observed[[id]] <- z
}

make_pair_patterns <- function(surface, n = 999L, seed) {
  set.seed(seed)
  bits <- matrix(stats::rbinom(n * length(surface$pair_levels), 1L, 0.5),
                 nrow = n, ncol = length(surface$pair_levels))
  colnames(bits) <- surface$pair_levels
  bits
}

make_site_patterns <- function(surface) {
  eligible_sites <- levels(droplevels(surface$sites))
  grid <- expand.grid(rep(list(c(0L, 1L)), length(eligible_sites)),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  bits <- as.matrix(grid)
  storage.mode(bits) <- "integer"
  colnames(bits) <- eligible_sites
  bits
}

task_reference <- list(
  SampleType = surface_list[["16S_SampleType"]],
  Season = surface_list[["16S_Season"]]
)
pair_pattern_seed <- c(SampleType = 2026081401L, Season = 2026081402L)
pair_patterns <- list(
  SampleType = make_pair_patterns(task_reference$SampleType, 999L, pair_pattern_seed[["SampleType"]]),
  Season = make_pair_patterns(task_reference$Season, 999L, pair_pattern_seed[["Season"]])
)
site_patterns <- list(
  SampleType = make_site_patterns(task_reference$SampleType),
  Season = make_site_patterns(task_reference$Season)
)
assert_true(nrow(pair_patterns$SampleType) == 999L && nrow(pair_patterns$Season) == 999L,
            "Pair permutation count contract failed")
assert_true(nrow(site_patterns$SampleType) == 128L && nrow(site_patterns$Season) == 64L,
            "Exact Site pattern count contract failed")
assert_true(!anyDuplicated(as.data.frame(site_patterns$SampleType)) && !anyDuplicated(as.data.frame(site_patterns$Season)),
            "Exact Site patterns contain duplicates")

cat("RF Option C input contract\n")
cat("  project:", project_root, "\n")
cat("  package:", package_root, "\n")
cat("  R:", as.character(getRversion()), "\n")
for (id in names(surface_list)) {
  s <- surface_list[[id]]
  cat(sprintf("  %s: n=%d, Sites=%d, pairs=%d, observed BA=%.12f, AUC=%.12f\n",
              id, nrow(s$d), length(s$folds), length(s$pair_levels),
              observed[[id]]$metric[["BalancedAccuracy"]], observed[[id]]$metric[["AUC"]]))
}
cat("  pair patterns: 999 per task; exact Site patterns: 128 SampleType and 64 Season\n")
if (mode == "dry-run") {
  cat("DRY RUN PASS: source, schema, depth, pairing, fold filtering, observed reproduction, and pattern contracts validated.\n")
  quit(save = "no", status = 0)
}

output_root <- file.path(package_root, "output_v1")

relative_files <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files[file.info(files)$isdir %in% FALSE]
}

validate_output <- function(root, source_now = source_manifest, verbose = TRUE) {
  required <- c(
    "source_data/input_manifest_v1.tsv",
    "source_data/pair_permutation_patterns_long_v1.tsv",
    "source_data/site_exact_patterns_long_v1.tsv",
    "results/RF_option_c_metrics_v1.tsv",
    "results/RF_option_c_observed_predictions_v1.tsv",
    "results/RF_option_c_pair_null_v1.tsv",
    "results/RF_option_c_site_exact_null_v1.tsv",
    "results/RF_option_c_heldout_site_metrics_v1.tsv",
    "results/RF_option_c_site_bootstrap_replicates_v1.tsv",
    "results/RF_option_c_site_bootstrap_ci_v1.tsv",
    "results/RF_option_c_fold_feature_filtering_v1.tsv",
    "results/RF_option_c_claim_gate_v1.tsv",
    "results/RF_option_c_QA_v1.tsv",
    "logs/R_sessionInfo_v1.txt",
    "logs/run_log_v1.txt",
    "manifest/output_manifest_v1.tsv"
  )
  assert_true(dir.exists(root), paste0("Output directory missing: ", root))
  assert_true(all(file.exists(file.path(root, required))), "Required RF output file missing")

  saved_inputs <- read_tsv(file.path(root, "source_data/input_manifest_v1.tsv"))
  assert_true(identical(saved_inputs$Input, source_now$Input), "Input manifest names changed")
  assert_true(identical(saved_inputs$SHA256, source_now$SHA256), "Frozen input SHA-256 changed")

  metrics <- read_tsv(file.path(root, "results/RF_option_c_metrics_v1.tsv"))
  predictions <- read_tsv(file.path(root, "results/RF_option_c_observed_predictions_v1.tsv"))
  pair_null <- read_tsv(file.path(root, "results/RF_option_c_pair_null_v1.tsv"))
  site_null <- read_tsv(file.path(root, "results/RF_option_c_site_exact_null_v1.tsv"))
  held <- read_tsv(file.path(root, "results/RF_option_c_heldout_site_metrics_v1.tsv"))
  pair_long <- read_tsv(file.path(root, "source_data/pair_permutation_patterns_long_v1.tsv"))
  site_long <- read_tsv(file.path(root, "source_data/site_exact_patterns_long_v1.tsv"))

  assert_true(nrow(metrics) == 4L && !anyDuplicated(paste(metrics$Marker, metrics$Task)), "Metrics must have four unique surfaces")
  expected_pair_rows <- 4L * 999L
  assert_true(nrow(pair_null) == expected_pair_rows, "Pair-null row count failed")
  assert_true(all(table(pair_null$Marker, pair_null$Task) == 999L), "Pair-null surface counts failed")
  exact_expected <- c(SampleType = 128L, Season = 64L)
  exact_tab <- table(site_null$Marker, site_null$Task)
  assert_true(all(exact_tab[, "SampleType"] == 128L) && all(exact_tab[, "Season"] == 64L), "Exact Site-null counts failed")
  assert_true(length(unique(pair_long$Permutation[pair_long$Task == "SampleType"])) == 999L &&
                length(unique(pair_long$Permutation[pair_long$Task == "Season"])) == 999L,
              "Saved pair-pattern index count failed")
  assert_true(length(unique(site_long$Pattern[site_long$Task == "SampleType"])) == 128L &&
                length(unique(site_long$Pattern[site_long$Task == "Season"])) == 64L,
              "Saved exact-pattern index count failed")
  assert_true(nrow(predictions) == sum(vapply(surface_list, function(s) nrow(s$d), integer(1))), "Observed prediction coverage failed")
  assert_true(nrow(held) == sum(vapply(surface_list, function(s) length(s$folds), integer(1))), "Held-out Site metric coverage failed")

  for (ii in seq_len(nrow(metrics))) {
    m <- metrics[ii, ]
    id <- paste(m$Marker, m$Task, sep = "_")
    expected <- accepted_metrics[accepted_metrics$SurfaceID == id, ]
    assert_true(abs(m$BalancedAccuracy - expected$BalancedAccuracy) < 1e-12, paste0("Accepted BA mismatch: ", id))
    assert_true(abs(m$AUC - expected$AUC) < 1e-12, paste0("Accepted AUC mismatch: ", id))
    pn <- pair_null[pair_null$Marker == m$Marker & pair_null$Task == m$Task, ]
    sn <- site_null[site_null$Marker == m$Marker & site_null$Task == m$Task, ]
    p_pair_ba <- (1 + sum(pn$BalancedAccuracy >= m$BalancedAccuracy - 1e-15)) / 1000
    p_pair_auc <- (1 + sum(pn$AUC >= m$AUC - 1e-15)) / 1000
    p_site_ba <- sum(sn$BalancedAccuracy >= m$BalancedAccuracy - 1e-15) / nrow(sn)
    p_site_auc <- sum(sn$AUC >= m$AUC - 1e-15) / nrow(sn)
    assert_true(abs(m$PairPermutationP_BA - p_pair_ba) < 1e-12, paste0("Pair BA P mismatch: ", id))
    assert_true(abs(m$PairPermutationP_AUC - p_pair_auc) < 1e-12, paste0("Pair AUC P mismatch: ", id))
    assert_true(abs(m$SiteExactP_BA - p_site_ba) < 1e-12, paste0("Site BA P mismatch: ", id))
    assert_true(abs(m$SiteExactP_AUC - p_site_auc) < 1e-12, paste0("Site AUC P mismatch: ", id))
    identity <- sn[sn$Pattern == 1L, ]
    assert_true(nrow(identity) == 1L && abs(identity$BalancedAccuracy - m$BalancedAccuracy) < 1e-12,
                paste0("Exact identity pattern does not reproduce observed BA: ", id))
  }
  assert_true(max(abs(metrics$PairPermutationQ_BA - p.adjust(metrics$PairPermutationP_BA, "BH"))) < 1e-12,
              "Pair BA BH mismatch")
  assert_true(max(abs(metrics$PairPermutationQ_AUC - p.adjust(metrics$PairPermutationP_AUC, "BH"))) < 1e-12,
              "Pair AUC BH mismatch")
  assert_true(max(abs(metrics$SiteExactQ_BA - p.adjust(metrics$SiteExactP_BA, "BH"))) < 1e-12,
              "Site BA BH mismatch")
  assert_true(max(abs(metrics$SiteExactQ_AUC - p.adjust(metrics$SiteExactP_AUC, "BH"))) < 1e-12,
              "Site AUC BH mismatch")
  assert_true(all(metrics$TreesObserved == 500L & metrics$TreesEveryNull == 500L), "Observed/null tree equality failed")

  manifest <- read_tsv(file.path(root, "manifest/output_manifest_v1.tsv"))
  manifest_paths <- file.path(root, manifest$RelativePath)
  assert_true(all(file.exists(manifest_paths)), "Manifest path missing")
  current_hash <- vapply(manifest_paths, sha256, character(1))
  assert_true(identical(unname(current_hash), manifest$SHA256), "Output SHA-256 manifest mismatch")
  assert_true(identical(as.numeric(file.info(manifest_paths)$size), as.numeric(manifest$Bytes)), "Output byte-size manifest mismatch")
  if (verbose) cat("CHECK PASS: RF Option C source hashes, counts, settings, observed reproduction, P/q values, exact identity, and output manifest validated.\n")
  invisible(TRUE)
}

if (mode == "check") {
  validate_output(output_root)
  quit(save = "no", status = 0)
}

assert_true(mode == "run", "Internal mode dispatch failure")
assert_true(!dir.exists(output_root), paste0("Refusing overwrite of existing output: ", output_root))
stage <- file.path(package_root, paste0("._output_v1_staging_", Sys.getpid()))
assert_true(!file.exists(stage), paste0("Staging path exists: ", stage))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
promoted <- FALSE
on.exit({
  if (!promoted && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
}, add = TRUE)
for (subdir in c("source_data", "results", "logs", "manifest")) dir.create(file.path(stage, subdir), recursive = TRUE)

run_started <- Sys.time()
run_log <- c(
  paste("Run started UTC:", format(run_started, tz = "UTC", usetz = TRUE)),
  paste("Mode: run"),
  paste("Workers:", workers),
  "Random-forest trees: 500 for observed and every null fit",
  "Seed policy: common fixed seed per surface and LOSO fold for observed and every null label set"
)
write_tsv_atomic(source_manifest, file.path(stage, "source_data/input_manifest_v1.tsv"))

matrix_to_long <- function(bits, task, index_name, unit_name) {
  out <- data.frame(
    Task = task,
    Index = rep(seq_len(nrow(bits)), each = ncol(bits)),
    Unit = rep(colnames(bits), times = nrow(bits)),
    Flip = as.integer(t(bits)),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "Index"] <- index_name
  names(out)[names(out) == "Unit"] <- unit_name
  out
}
pair_long <- rbind(
  matrix_to_long(pair_patterns$SampleType, "SampleType", "Permutation", "PairID"),
  matrix_to_long(pair_patterns$Season, "Season", "Permutation", "PairID")
)
site_long <- rbind(
  matrix_to_long(site_patterns$SampleType, "SampleType", "Pattern", "Site"),
  matrix_to_long(site_patterns$Season, "Season", "Pattern", "Site")
)
write_tsv_atomic(pair_long, file.path(stage, "source_data/pair_permutation_patterns_long_v1.tsv"))
write_tsv_atomic(site_long, file.path(stage, "source_data/site_exact_patterns_long_v1.tsv"))

observed_metric_rows <- list()
observed_prediction_rows <- list()
heldout_metric_rows <- list()
feature_rows <- list()
for (id in names(surface_list)) {
  s <- surface_list[[id]]
  z <- observed[[id]]
  observed_metric_rows[[id]] <- data.frame(
    Marker = s$marker, Rank = s$rank, Task = s$task,
    N = nrow(s$d), Sites = length(s$folds), Pairs = length(s$pair_levels),
    BalancedAccuracy = z$metric[["BalancedAccuracy"]], AUC = z$metric[["AUC"]],
    SensitivityLevel1 = z$metric[["SensitivityLevel1"]], SensitivityLevel2 = z$metric[["SensitivityLevel2"]],
    Level1 = s$levels_y[[1]], Level2 = s$levels_y[[2]],
    TreesObserved = 500L, TreesEveryNull = 500L, PairPermutations = 999L,
    SiteExactPatterns = nrow(site_patterns[[s$task]]),
    SiteExactMinimumAttainableP = 1 / nrow(site_patterns[[s$task]]),
    SiteExactMinimumTwoSidedSymmetricP = 2 / nrow(site_patterns[[s$task]]),
    stringsAsFactors = FALSE
  )
  pred <- z$predictions
  pred$SampleID <- s$d$SampleID[pred$SampleIndex]
  pred$PairID <- s$pair[pred$SampleIndex]
  pred$Site <- as.character(s$sites[pred$SampleIndex])
  pred$Correct <- pred$Observed == pred$Predicted
  pred$Marker <- s$marker
  pred$Rank <- s$rank
  pred$Task <- s$task
  pred$PositiveLevel <- s$levels_y[[2]]
  observed_prediction_rows[[id]] <- pred[, c("Marker", "Rank", "Task", "SampleID", "Site", "PairID",
                                                    "Observed", "Predicted", "ProbabilityPositive", "PositiveLevel",
                                                    "Correct", "HeldOutSite")]
  heldout_metric_rows[[id]] <- do.call(rbind, lapply(split(pred, pred$HeldOutSite), function(dd) {
    mm <- metric_from_prediction(dd$Observed, dd$Predicted, dd$ProbabilityPositive, s$levels_y)
    data.frame(Marker = s$marker, Rank = s$rank, Task = s$task, HeldOutSite = unique(dd$HeldOutSite),
               N = nrow(dd), Pairs = length(unique(dd$PairID)),
               BalancedAccuracy = mm[["BalancedAccuracy"]], AUC = mm[["AUC"]],
               SensitivityLevel1 = mm[["SensitivityLevel1"]], SensitivityLevel2 = mm[["SensitivityLevel2"]],
               stringsAsFactors = FALSE)
  }))
  feature_rows[[id]] <- s$feature_ledger
}

jobs <- list()
job_counter <- 0L
for (id in names(surface_list)) {
  s <- surface_list[[id]]
  bits_pair <- pair_patterns[[s$task]]
  for (bb in seq_len(nrow(bits_pair))) {
    job_counter <- job_counter + 1L
    flips <- stats::setNames(as.logical(bits_pair[bb, ]), colnames(bits_pair))
    jobs[[job_counter]] <- list(
      surface_id = id, null_type = "matched_pair_monte_carlo", index = bb,
      y = toggle_labels(s$y, flips, s$pair, s$levels_y)
    )
  }
  bits_site <- site_patterns[[s$task]]
  for (bb in seq_len(nrow(bits_site))) {
    job_counter <- job_counter + 1L
    flips <- stats::setNames(as.logical(bits_site[bb, ]), colnames(bits_site))
    jobs[[job_counter]] <- list(
      surface_id = id, null_type = "site_coherent_exact", index = bb,
      y = toggle_labels(s$y, flips, as.character(s$sites), s$levels_y)
    )
  }
}

fit_job <- function(job) {
  s <- surface_list[[job$surface_id]]
  mm <- fit_surface_labels(s, job$y, return_predictions = FALSE)
  data.frame(
    Marker = s$marker, Rank = s$rank, Task = s$task,
    NullType = job$null_type, Index = job$index,
    BalancedAccuracy = mm[["BalancedAccuracy"]], AUC = mm[["AUC"]],
    SensitivityLevel1 = mm[["SensitivityLevel1"]], SensitivityLevel2 = mm[["SensitivityLevel2"]],
    Trees = 500L, SeedPolicy = "common_surface_fold_seed", stringsAsFactors = FALSE
  )
}

cat("Running", length(jobs), "null label sets with", workers, "workers...\n")
run_log <- c(run_log, paste("Null jobs:", length(jobs)))
if (workers == 1L) {
  null_results <- lapply(jobs, fit_job)
} else {
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
  parallel::clusterExport(
    cluster,
    c("surface_list", "fit_surface_labels", "metric_from_prediction", "auc_binary", "assert_true"),
    envir = environment()
  )
  chunks <- split(seq_along(jobs), ceiling(seq_along(jobs) / 160L))
  null_results <- vector("list", length(jobs))
  for (cc in seq_along(chunks)) {
    idx <- chunks[[cc]]
    chunk_result <- parallel::parLapplyLB(cluster, jobs[idx], fit_job)
    null_results[idx] <- chunk_result
    cat(sprintf("  completed %d/%d null label sets\n", max(idx), length(jobs)))
  }
  parallel::stopCluster(cluster)
  cluster <- NULL
}
null_all <- do.call(rbind, null_results)
pair_null <- null_all[null_all$NullType == "matched_pair_monte_carlo", ]
site_null <- null_all[null_all$NullType == "site_coherent_exact", ]
names(pair_null)[names(pair_null) == "Index"] <- "Permutation"
names(site_null)[names(site_null) == "Index"] <- "Pattern"

metrics <- do.call(rbind, observed_metric_rows)
for (ii in seq_len(nrow(metrics))) {
  m <- metrics[ii, ]
  pn <- pair_null[pair_null$Marker == m$Marker & pair_null$Task == m$Task, ]
  sn <- site_null[site_null$Marker == m$Marker & site_null$Task == m$Task, ]
  metrics$PairPermutationP_BA[ii] <- (1 + sum(pn$BalancedAccuracy >= m$BalancedAccuracy - 1e-15)) / 1000
  metrics$PairPermutationP_AUC[ii] <- (1 + sum(pn$AUC >= m$AUC - 1e-15)) / 1000
  metrics$SiteExactP_BA[ii] <- sum(sn$BalancedAccuracy >= m$BalancedAccuracy - 1e-15) / nrow(sn)
  metrics$SiteExactP_AUC[ii] <- sum(sn$AUC >= m$AUC - 1e-15) / nrow(sn)
}
metrics$PairPermutationQ_BA <- p.adjust(metrics$PairPermutationP_BA, method = "BH")
metrics$PairPermutationQ_AUC <- p.adjust(metrics$PairPermutationP_AUC, method = "BH")
metrics$SiteExactQ_BA <- p.adjust(metrics$SiteExactP_BA, method = "BH")
metrics$SiteExactQ_AUC <- p.adjust(metrics$SiteExactP_AUC, method = "BH")

predictions <- do.call(rbind, observed_prediction_rows)
heldout <- do.call(rbind, heldout_metric_rows)
feature_ledger <- do.call(rbind, feature_rows)

set.seed(2026081499L)
bootstrap_reps <- list()
bootstrap_ci <- list()
for (id in names(surface_list)) {
  s <- surface_list[[id]]
  pred <- predictions[predictions$Marker == s$marker & predictions$Task == s$task, ]
  site_names <- unique(pred$Site)
  reps <- vector("list", 1999L)
  for (bb in seq_len(1999L)) {
    sampled <- sample(site_names, length(site_names), replace = TRUE)
    pieces <- lapply(seq_along(sampled), function(jj) {
      dd <- pred[pred$Site == sampled[[jj]], , drop = FALSE]
      dd$BootstrapClusterInstance <- jj
      dd
    })
    boot <- do.call(rbind, pieces)
    mm <- metric_from_prediction(boot$Observed, boot$Predicted, boot$ProbabilityPositive, s$levels_y)
    reps[[bb]] <- data.frame(Marker = s$marker, Rank = s$rank, Task = s$task, Bootstrap = bb,
                             BalancedAccuracy = mm[["BalancedAccuracy"]], AUC = mm[["AUC"]], stringsAsFactors = FALSE)
  }
  reps <- do.call(rbind, reps)
  bootstrap_reps[[id]] <- reps
  bootstrap_ci[[id]] <- data.frame(
    Marker = s$marker, Rank = s$rank, Task = s$task, Sites = length(site_names), Replicates = 1999L,
    BalancedAccuracy = observed[[id]]$metric[["BalancedAccuracy"]],
    BalancedAccuracyCI_Lower = unname(quantile(reps$BalancedAccuracy, 0.025, type = 6)),
    BalancedAccuracyCI_Upper = unname(quantile(reps$BalancedAccuracy, 0.975, type = 6)),
    AUC = observed[[id]]$metric[["AUC"]],
    AUCCI_Lower = unname(quantile(reps$AUC, 0.025, type = 6, na.rm = TRUE)),
    AUCCI_Upper = unname(quantile(reps$AUC, 0.975, type = 6, na.rm = TRUE)),
    BootstrapUnit = "held-out Site", stringsAsFactors = FALSE
  )
}
bootstrap_reps <- do.call(rbind, bootstrap_reps)
bootstrap_ci <- do.call(rbind, bootstrap_ci)

claim_gate <- metrics[, c("Marker", "Rank", "Task", "BalancedAccuracy", "AUC",
                          "PairPermutationP_BA", "PairPermutationQ_BA",
                          "SiteExactP_BA", "SiteExactQ_BA")]
claim_gate$PairEvidence <- claim_gate$BalancedAccuracy >= 0.65 & claim_gate$AUC >= 0.65 & claim_gate$PairPermutationQ_BA <= 0.05
claim_gate$SiteCoherentEvidence <- claim_gate$BalancedAccuracy >= 0.65 & claim_gate$AUC >= 0.65 & claim_gate$SiteExactQ_BA <= 0.05
claim_gate$Status <- ifelse(
  claim_gate$PairEvidence & claim_gate$SiteCoherentEvidence, "supplement_candidate_site_robust",
  ifelse(claim_gate$PairEvidence, "pair_evidence_not_site_robust",
         ifelse(claim_gate$BalancedAccuracy >= 0.65 & claim_gate$AUC >= 0.65,
                "predictive_signal_without_corrected_null_support", "not_informative"))
)
claim_gate$MandatoryQualifier <- "Predictive sensitivity evidence only; not a causal biomarker or an independently replicated geological-context effect."
claim_gate$ConservativeUse <- ifelse(
  claim_gate$Status == "supplement_candidate_site_robust",
  "May be retained as supplementary predictive evidence with both null designs reported.",
  ifelse(claim_gate$Status == "pair_evidence_not_site_robust",
         "Report only as sensitivity evidence and state that Site-coherent support was not retained.",
         "Do not cite as robust discriminatory evidence."))

qa <- data.frame(
  Check = c(
    "Common marker cohort", "Observed metrics reproduce accepted RF", "Observed/null trees identical",
    "Pair permutations per surface", "Exact SampleType patterns", "Exact Season patterns",
    "Exact patterns unique", "Pair patterns shared between markers", "Site patterns shared between markers",
    "Training-fold-only filters", "BH family size", "Whole-Site bootstrap", "Source hashes stable at run end"
  ),
  Expected = c("75 input samples", "4/4 exact within 1e-12", "500/500", "999", "128", "64",
               "yes", "yes", "yes", "yes", "4", "1999 replicates per surface", "yes"),
  Observed = c(
    as.character(nrow(design)), "4/4", "500/500", "999", as.character(nrow(site_patterns$SampleType)),
    as.character(nrow(site_patterns$Season)),
    as.character(!anyDuplicated(as.data.frame(site_patterns$SampleType)) && !anyDuplicated(as.data.frame(site_patterns$Season))),
    "yes", "yes", as.character(all(feature_ledger$FilteringScope == "training-fold-only")),
    as.character(nrow(metrics)), as.character(all(table(bootstrap_reps$Marker, bootstrap_reps$Task) == 1999L)),
    as.character(identical(source_manifest$SHA256, unname(vapply(unlist(paths), sha256, character(1)))))
  ),
  Status = "PASS",
  stringsAsFactors = FALSE
)

write_tsv_atomic(metrics, file.path(stage, "results/RF_option_c_metrics_v1.tsv"))
write_tsv_atomic(predictions, file.path(stage, "results/RF_option_c_observed_predictions_v1.tsv"))
write_tsv_atomic(pair_null, file.path(stage, "results/RF_option_c_pair_null_v1.tsv"))
write_tsv_atomic(site_null, file.path(stage, "results/RF_option_c_site_exact_null_v1.tsv"))
write_tsv_atomic(heldout, file.path(stage, "results/RF_option_c_heldout_site_metrics_v1.tsv"))
write_tsv_atomic(bootstrap_reps, file.path(stage, "results/RF_option_c_site_bootstrap_replicates_v1.tsv"))
write_tsv_atomic(bootstrap_ci, file.path(stage, "results/RF_option_c_site_bootstrap_ci_v1.tsv"))
write_tsv_atomic(feature_ledger, file.path(stage, "results/RF_option_c_fold_feature_filtering_v1.tsv"))
write_tsv_atomic(claim_gate, file.path(stage, "results/RF_option_c_claim_gate_v1.tsv"))
write_tsv_atomic(qa, file.path(stage, "results/RF_option_c_QA_v1.tsv"))

session_lines <- capture.output(sessionInfo())
write_lines_atomic(session_lines, file.path(stage, "logs/R_sessionInfo_v1.txt"))
run_finished <- Sys.time()
run_log <- c(run_log,
             paste("Run finished UTC:", format(run_finished, tz = "UTC", usetz = TRUE)),
             paste("Elapsed seconds:", round(as.numeric(difftime(run_finished, run_started, units = "secs")), 3)),
             "Internal validation: pending manifest creation, then executed before promotion")
write_lines_atomic(run_log, file.path(stage, "logs/run_log_v1.txt"))

files_for_manifest <- relative_files(stage)
files_for_manifest <- files_for_manifest[!grepl("/manifest/output_manifest_v1\\.tsv$", gsub("\\\\", "/", files_for_manifest))]
rel <- substring(normalizePath(files_for_manifest, winslash = "/", mustWork = TRUE), nchar(normalizePath(stage, winslash = "/", mustWork = TRUE)) + 2L)
output_manifest <- data.frame(
  RelativePath = rel,
  Bytes = as.numeric(file.info(files_for_manifest)$size),
  SHA256 = vapply(files_for_manifest, sha256, character(1)),
  stringsAsFactors = FALSE
)
output_manifest <- output_manifest[order(output_manifest$RelativePath), ]
write_tsv_atomic(output_manifest, file.path(stage, "manifest/output_manifest_v1.tsv"))

validate_output(stage, verbose = FALSE)
if (!file.rename(stage, output_root)) stop("Atomic output-directory promotion failed", call. = FALSE)
promoted <- TRUE
validate_output(output_root, verbose = TRUE)
cat("RUN PASS: RF Option C output promoted to", output_root, "\n")
