#!/usr/bin/env Rscript

# Figure-free cross-marker concordance analysis for the protected shared
# 75-sample 16S/18S revision data.

options(stringsAsFactors = FALSE, warn = 1, digits = 15)

SCRIPT_VERSION <- "1.0.0"
TARGET_NAME <- "09_cross_marker_concordance_v1"
N_PERMUTATIONS <- 9999L
N_BOOTSTRAPS <- 9999L
CI_LEVEL <- 0.95
TOLERANCE <- 1e-12
SUBSET_ORDER <- c("A75", "A72", "Sea-sand", "Water", "First", "Second")
EXPECTED_SUBSET_N <- c(A75 = 75L, A72 = 72L, `Sea-sand` = 39L,
                       Water = 36L, First = 42L, Second = 33L)
ALPHA_METRICS <- c("Observed", "Shannon", "Simpson")

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
assert_true <- function(value, message) if (!isTRUE(value)) stop(message, call. = FALSE)

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(normalizePath(".", winslash = "/", mustWork = TRUE))
  normalizePath(sub("^--file=", "", hit[[1L]]), winslash = "/", mustWork = TRUE)
}

default_base <- function() {
  normalizePath(file.path(dirname(script_path()), ".."), winslash = "/", mustWork = TRUE)
}

usage <- function() {
  cat(paste0(
    "run_cross_marker_concordance_v1.R ", SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript run_cross_marker_concordance_v1.R --dry-run [--base PATH]\n",
    "  Rscript run_cross_marker_concordance_v1.R --run     [--base PATH]\n",
    "  Rscript run_cross_marker_concordance_v1.R --check   [--base PATH]\n\n",
    "Exactly one mode is required. --dry-run validates inputs without writing.\n",
    "--run refuses to overwrite 09_cross_marker_concordance_v1, builds in a\n",
    "sibling staging directory, and commits by atomic directory rename.\n",
    "--check validates schemas, hashes, row identities, and numerical ledgers.\n"
  ))
}

parse_cli <- function(args) {
  if (any(args %in% c("--help", "-h"))) return(list(help = TRUE))
  modes <- intersect(args, c("--dry-run", "--run", "--check"))
  if (length(modes) != 1L) stopf("Specify exactly one of --dry-run, --run, or --check")
  base <- default_base()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--dry-run", "--run", "--check")) {
      i <- i + 1L
    } else if (arg == "--base") {
      if (i == length(args)) stopf("--base requires a path")
      base <- args[[i + 1L]]
      i <- i + 2L
    } else if (startsWith(arg, "--base=")) {
      base <- sub("^--base=", "", arg)
      i <- i + 1L
    } else stopf("Unknown argument: %s", arg)
  }
  assert_true(dir.exists(base), sprintf("Base directory does not exist: %s", base))
  list(help = FALSE, mode = sub("^--", "", modes[[1L]]),
       base = normalizePath(base, winslash = "/", mustWork = TRUE))
}

check_packages <- function() {
  needed <- c("vegan", "permute", "digest")
  missing <- needed[!vapply(needed, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L))]
  if (length(missing)) stopf("Missing required R package(s): %s", paste(missing, collapse = ", "))
}

read_tsv <- function(path) {
  assert_true(file.exists(path), sprintf("Required file is missing: %s", path))
  tryCatch(utils::read.delim(path, header = TRUE, sep = "\t", quote = "",
    comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = "NA", fileEncoding = "UTF-8-BOM"),
    error = function(e) stopf("Cannot read %s: %s", path, conditionMessage(e)))
}

assert_schema <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stopf("%s is missing column(s): %s", label, paste(missing, collapse = ", "))
}

assert_unique <- function(x, label) {
  if (anyNA(x) || any(!nzchar(as.character(x)))) stopf("%s contains blank/NA identifiers", label)
  dup <- unique(x[duplicated(x)])
  if (length(dup)) stopf("%s contains duplicate identifier(s): %s", label,
                         paste(head(dup, 10L), collapse = ", "))
}

sha256_file <- function(path) {
  assert_true(file.exists(path), sprintf("Cannot hash missing file: %s", path))
  tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
file_size <- function(path) as.numeric(file.info(path)$size)

verify_package_manifest <- function(root, manifest_rel) {
  manifest_path <- file.path(root, manifest_rel)
  manifest <- read_tsv(manifest_path)
  assert_schema(manifest, c("RelativePath", "SizeBytes", "SHA256"),
                paste(basename(root), "manifest"))
  assert_unique(manifest$RelativePath, paste(basename(root), "manifest RelativePath"))
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(root, manifest$RelativePath[[i]])
    assert_true(file.exists(path), sprintf("Upstream manifest target missing: %s", path))
    assert_true(file_size(path) == as.numeric(manifest$SizeBytes[[i]]),
                sprintf("Upstream size mismatch: %s", path))
    assert_true(identical(sha256_file(path), tolower(manifest$SHA256[[i]])),
                sprintf("Upstream SHA256 mismatch: %s", path))
  }
  actual <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  actual <- actual[file.info(actual)$isdir %in% FALSE]
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  actual_rel <- substring(normalizePath(actual, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L)
  expected <- c(manifest$RelativePath, gsub("\\\\", "/", manifest_rel))
  assert_true(setequal(actual_rel, expected), sprintf("Upstream manifest coverage drift: %s", root))
  invisible(TRUE)
}

input_paths <- function(base) c(
  analysis_script = script_path(),
  master_metadata = file.path(base, "00_design_v2", "master_sample_metadata_v2.tsv"),
  set_membership = file.path(base, "00_design_v2", "analysis_set_membership_v2.tsv"),
  alpha_metrics = file.path(base, "05_alpha_diversity_v1", "results", "alpha_sample_metrics_v1.tsv"),
  alpha_manifest = file.path(base, "05_alpha_diversity_v1", "provenance", "output_sha256_manifest_v1.tsv"),
  beta_pcoa = file.path(base, "06_beta_diversity_v2", "results", "beta_pcoa_coordinates_lingoes_v2.tsv"),
  beta_registry = file.path(base, "06_beta_diversity_v2", "results", "beta_surface_registry_v2.tsv"),
  beta_manifest = file.path(base, "06_beta_diversity_v2", "provenance", "output_sha256_manifest_v2.tsv"),
  counts_16s = file.path(base, "04_input_16s_rarefied_20000_v1", "results",
                         "16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  counts_18s = file.path(base, "03_18s_rarefied_50000_v1", "results",
                         "18s_asv_counts_rarefied_50000_analysis_v1.tsv")
)

read_count_table <- function(path, expected_samples, marker, expected_depth) {
  x <- read_tsv(path)
  assert_true(ncol(x) >= 2L, sprintf("%s count table has fewer than two columns", marker))
  assert_true(names(x)[[1L]] %in% c("#OTUID", "OTUID", "#OTU ID", "ASVID", "FeatureID"),
              sprintf("Unexpected feature-ID header in %s", path))
  assert_unique(as.character(x[[1L]]), paste(marker, "feature IDs"))
  sample_cols <- names(x)[-1L]
  assert_unique(sample_cols, paste(marker, "sample columns"))
  assert_true(setequal(sample_cols, expected_samples),
              sprintf("%s sample columns do not match master metadata", marker))
  m <- as.matrix(x[, match(expected_samples, sample_cols) + 1L, drop = FALSE])
  suppressWarnings(storage.mode(m) <- "numeric")
  assert_true(!anyNA(m) && all(is.finite(m)) && all(m >= 0) &&
                all(abs(m - round(m)) <= TOLERANCE),
              sprintf("%s table must contain finite nonnegative integer counts", marker))
  m <- t(round(m)); rownames(m) <- expected_samples
  depths <- rowSums(m)
  assert_true(all(depths == expected_depth),
              sprintf("%s table is not exactly %d reads/sample", marker, expected_depth))
  m
}

subset_definitions <- function(master, membership) {
  a72 <- membership[as.character(membership$SetID) == "A72_SUBSTRATE_PAIRED_CELLS" &
                      as.character(membership$UnitType) == "sample" &
                      as.character(membership$Included) == "TRUE", "SampleID"]
  a72 <- master$SampleID[master$SampleID %in% a72]
  ids <- list(
    A75 = master$SampleID,
    A72 = a72,
    `Sea-sand` = master$SampleID[master$SampleType == "Sea-sand"],
    Water = master$SampleID[master$SampleType == "Water"],
    First = master$SampleID[master$Campaign == "First"],
    Second = master$SampleID[master$Campaign == "Second"]
  )
  defs <- c(
    A75 = "all protected shared 75 samples",
    A72 = "A72 complete substrate-comparison cells; no suffix-derived PairIDs",
    `Sea-sand` = "A75 samples with SampleType=Sea-sand",
    Water = "A75 samples with SampleType=Water",
    First = "A75 samples with Campaign=First",
    Second = "A75 samples with Campaign=Second"
  )
  for (nm in SUBSET_ORDER) {
    assert_unique(ids[[nm]], paste(nm, "subset SampleID"))
    assert_true(length(ids[[nm]]) == EXPECTED_SUBSET_N[[nm]],
                sprintf("Subset %s count is %d, expected %d", nm, length(ids[[nm]]), EXPECTED_SUBSET_N[[nm]]))
  }
  list(ids = ids, definitions = defs)
}

load_inputs <- function(base) {
  paths <- input_paths(base)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stopf("Required input(s) missing:\n%s",
                             paste(sprintf("  %s: %s", names(missing), missing), collapse = "\n"))
  verify_package_manifest(file.path(base, "05_alpha_diversity_v1"),
                          file.path("provenance", "output_sha256_manifest_v1.tsv"))
  verify_package_manifest(file.path(base, "06_beta_diversity_v2"),
                          file.path("provenance", "output_sha256_manifest_v2.tsv"))
  master <- read_tsv(paths[["master_metadata"]])
  membership <- read_tsv(paths[["set_membership"]])
  alpha <- read_tsv(paths[["alpha_metrics"]])
  beta_pcoa <- read_tsv(paths[["beta_pcoa"]])
  beta_registry <- read_tsv(paths[["beta_registry"]])
  assert_schema(master, c("SampleID", "Site", "GeologicalHabitat", "SampleType", "Campaign",
                          "Included16S", "Included18S", "PairID_Substrate", "PairID_Season"), "master metadata")
  assert_unique(master$SampleID, "master SampleID")
  assert_true(nrow(master) == 75L, "Master metadata must contain exactly 75 samples")
  assert_true(all(as.character(master$Included16S) == "TRUE") &&
                all(as.character(master$Included18S) == "TRUE"), "Master inclusion flags are not all TRUE")
  assert_true(setequal(master$Campaign, c("First", "Second")), "Campaign must be First/Second")
  assert_true(setequal(master$SampleType, c("Sea-sand", "Water")), "SampleType must be Sea-sand/Water")
  assert_true(!any(grepl("Autumn", unlist(master), fixed = TRUE), na.rm = TRUE), "Autumn label is prohibited")
  assert_true(all(is.na(master$PairID_Substrate) | !nzchar(master$PairID_Substrate)) &&
                all(is.na(master$PairID_Season) | !nzchar(master$PairID_Season)),
              "Suffix-derived PairIDs are prohibited")
  assert_schema(membership, c("SetID", "UnitType", "SampleID", "Included"), "set membership")
  subsets <- subset_definitions(master, membership)
  assert_schema(alpha, c("Marker", "SampleID", "Site", "GeologicalHabitat", "SampleType",
                         "Campaign", "Metric", "Value", "CountSurface", "RarefactionDepth"), "alpha metrics")
  alpha_use <- alpha[alpha$Metric %in% ALPHA_METRICS, , drop = FALSE]
  assert_true(nrow(alpha_use) == 2L * 75L * length(ALPHA_METRICS), "Alpha source row count drift")
  assert_true(setequal(alpha_use$Marker, c("16S", "18S")), "Alpha source marker drift")
  assert_true(all(is.finite(alpha_use$Value)), "Alpha source contains nonfinite values")
  assert_true(!any(duplicated(alpha_use[c("Marker", "SampleID", "Metric")])), "Duplicate alpha source key")
  assert_true(setequal(unique(alpha_use$SampleID), master$SampleID), "Alpha sample identity drift")
  assert_schema(beta_pcoa, c("SampleID", "Marker", "Surface", "SetID", "Axis", "Coordinate",
                             "Correction", "Site", "SampleType", "Campaign"), "beta v2 PCoA")
  assert_schema(beta_registry, c("Marker", "Surface", "SetID", "SampleN", "Distance",
                                 "Normalization", "Correction"), "beta v2 registry")
  beta_a75 <- beta_registry[beta_registry$Surface == "rarefied_primary" &
                              beta_registry$SetID == "A75_FULL", , drop = FALSE]
  assert_true(nrow(beta_a75) == 2L && setequal(beta_a75$Marker, c("16S", "18S")) &&
                all(beta_a75$SampleN == 75L) && all(beta_a75$Distance == "Bray-Curtis") &&
                all(beta_a75$Correction == "Lingoes"), "Canonical beta v2 A75 registry drift")
  counts <- list(
    `16S` = read_count_table(paths[["counts_16s"]], master$SampleID, "16S", 20000L),
    `18S` = read_count_table(paths[["counts_18s"]], master$SampleID, "18S", 50000L)
  )
  list(paths = paths, master = master, membership = membership, alpha = alpha_use,
       subsets = subsets, counts = counts)
}

cluster_bootstrap <- function(x, y, cluster, nboot, seed) {
  lev <- unique(as.character(cluster))
  set.seed(seed)
  vals <- rep(NA_real_, nboot)
  for (b in seq_len(nboot)) {
    draw <- sample(lev, length(lev), replace = TRUE)
    idx <- unlist(lapply(draw, function(z) which(as.character(cluster) == z)), use.names = FALSE)
    vals[[b]] <- suppressWarnings(stats::cor(x[idx], y[idx], method = "spearman"))
  }
  valid <- vals[is.finite(vals)]
  assert_true(length(valid) >= floor(0.95 * nboot), "Too few valid site-cluster bootstrap replicates")
  alpha <- (1 - CI_LEVEL) / 2
  q <- as.numeric(stats::quantile(valid, probs = c(alpha, 1 - alpha), names = FALSE,
                                  type = 6, na.rm = TRUE))
  list(lower = q[[1L]], upper = q[[2L]], valid = length(valid))
}

compute_alpha <- function(inputs) {
  pair_parts <- list(); stat_parts <- list(); k <- 0L; s <- 0L
  for (subset_index in seq_along(SUBSET_ORDER)) {
    subset <- SUBSET_ORDER[[subset_index]]; ids <- inputs$subsets$ids[[subset]]
    meta <- inputs$master[match(ids, inputs$master$SampleID), , drop = FALSE]
    for (metric_index in seq_along(ALPHA_METRICS)) {
      metric <- ALPHA_METRICS[[metric_index]]
      a16 <- inputs$alpha[inputs$alpha$Marker == "16S" & inputs$alpha$Metric == metric, , drop = FALSE]
      a18 <- inputs$alpha[inputs$alpha$Marker == "18S" & inputs$alpha$Metric == metric, , drop = FALSE]
      x <- a16$Value[match(ids, a16$SampleID)]; y <- a18$Value[match(ids, a18$SampleID)]
      assert_true(all(is.finite(x)) && all(is.finite(y)), "Alpha pairing produced nonfinite values")
      k <- k + 1L
      pair_parts[[k]] <- data.frame(
        AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
        SampleID = ids, Site = meta$Site, GeologicalHabitat = meta$GeologicalHabitat,
        SampleType = meta$SampleType, Campaign = meta$Campaign, Metric = metric,
        Value16S = x, Value18S = y, CountSurface16S = "rarefied_primary_20000",
        CountSurface18S = "rarefied_primary_50000",
        InferenceScope = "descriptive cross-marker association; no causal or season claim",
        stringsAsFactors = FALSE)
      seed <- 73100L + subset_index * 100L + metric_index
      boot <- cluster_bootstrap(x, y, meta$Site, N_BOOTSTRAPS, seed)
      test <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
      s <- s + 1L
      stat_parts[[s]] <- data.frame(
        AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
        Metric = metric, SampleN = length(ids), SiteN = length(unique(meta$Site)),
        RhoObserved = unname(test$estimate), PAsymptoticDescriptive = test$p.value,
        BootstrapUnit = "Site", BootstrapReplicatesRequested = N_BOOTSTRAPS,
        BootstrapReplicatesValid = boot$valid, BootstrapSeed = seed,
        BootstrapCILevel = CI_LEVEL, BootstrapCILower = boot$lower, BootstrapCIUpper = boot$upper,
        InferenceScope = "descriptive; site-cluster bootstrap CI; asymptotic p is not cluster-adjusted",
        stringsAsFactors = FALSE)
    }
  }
  list(pairs = do.call(rbind, pair_parts), stats = do.call(rbind, stat_parts))
}

total_admissible <- function(site) {
  sizes <- as.numeric(table(site)); log_total <- sum(lgamma(sizes + 1))
  total <- exp(log_total)
  nonidentity_numeric <- if (is.finite(total)) total - 1 else Inf
  label <- if (is.finite(nonidentity_numeric) && nonidentity_numeric <= 9e15) {
    format(round(nonidentity_numeric), scientific = FALSE, trim = TRUE)
  } else format(signif(nonidentity_numeric, 15), scientific = TRUE, trim = TRUE)
  list(log_total = log_total, nonidentity_numeric = nonidentity_numeric, label = label)
}

restricted_permutations <- function(site, requested, seed) {
  site <- factor(site, levels = unique(site))
  admissible <- total_admissible(site)
  control <- permute::how(within = permute::Within(type = "free"), blocks = site,
                          observed = FALSE, nperm = requested)
  identity <- seq_along(site)
  if (is.finite(admissible$nonidentity_numeric) && admissible$nonidentity_numeric <= requested) {
    p <- as.matrix(permute::allPerms(length(site), control = control))
    mode <- "exhaustive_all_nonidentity"
  } else {
    set.seed(seed)
    p <- as.matrix(permute::shuffleSet(length(site), nset = requested, control = control,
                                       quietly = TRUE))
    mode <- "monte_carlo_without_identity"
  }
  if (!nrow(p)) stopf("No restricted permutations generated")
  p <- p[!apply(p, 1L, function(z) identical(as.integer(z), identity)), , drop = FALSE]
  assert_true(nrow(p) > 0L, "Only identity permutation was generated")
  assert_true(all(apply(p, 1L, function(z) all(as.character(site[z]) == as.character(site)))),
              "Permutation matrix violates Site restriction")
  if (mode == "monte_carlo_without_identity")
    assert_true(nrow(p) == requested, "Restricted Monte Carlo permutation count drift")
  list(matrix = p, mode = mode, total = admissible$label)
}

compute_beta_concordance <- function(inputs) {
  full_dist <- lapply(inputs$counts, function(m) vegan::vegdist(m, method = "bray"))
  pair_parts <- list(); bray_stats <- list(); proc_scores <- list(); proc_stats <- list()
  for (i in seq_along(SUBSET_ORDER)) {
    subset <- SUBSET_ORDER[[i]]; ids <- inputs$subsets$ids[[subset]]
    meta <- inputs$master[match(ids, inputs$master$SampleID), , drop = FALSE]
    d16 <- as.dist(as.matrix(full_dist$`16S`)[ids, ids, drop = FALSE])
    d18 <- as.dist(as.matrix(full_dist$`18S`)[ids, ids, drop = FALSE])
    m16 <- as.matrix(d16); m18 <- as.matrix(d18)
    cmb <- utils::combn(seq_along(ids), 2L)
    v16 <- m16[cbind(cmb[1L, ], cmb[2L, ])]
    v18 <- m18[cbind(cmb[1L, ], cmb[2L, ])]
    pair_parts[[i]] <- data.frame(
      AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
      SampleID1 = ids[cmb[1L, ]], SampleID2 = ids[cmb[2L, ]],
      Site1 = meta$Site[cmb[1L, ]], Site2 = meta$Site[cmb[2L, ]],
      GeologicalHabitat1 = meta$GeologicalHabitat[cmb[1L, ]],
      GeologicalHabitat2 = meta$GeologicalHabitat[cmb[2L, ]],
      SampleType1 = meta$SampleType[cmb[1L, ]], SampleType2 = meta$SampleType[cmb[2L, ]],
      Campaign1 = meta$Campaign[cmb[1L, ]], Campaign2 = meta$Campaign[cmb[2L, ]],
      PairScopeSite = ifelse(meta$Site[cmb[1L, ]] == meta$Site[cmb[2L, ]], "WithinSite", "BetweenSite"),
      Bray16S = v16, Bray18S = v18,
      DistanceSurface16S = "Bray-Curtis; rarefied_primary_20000",
      DistanceSurface18S = "Bray-Curtis; rarefied_primary_50000",
      stringsAsFactors = FALSE)

    bray_seed <- 84100L + i
    bp <- restricted_permutations(meta$Site, N_PERMUTATIONS, bray_seed)
    # A label permutation only reorders the pairwise distances, so pre-ranking
    # once is exactly equivalent to recomputing Spearman ranks 9,999 times and
    # avoids repeated sorting of thousands of distance pairs.
    rank16 <- rank(v16, ties.method = "average")
    rank18 <- rank(v18, ties.method = "average")
    rank18_matrix <- matrix(0, nrow = length(ids), ncol = length(ids))
    rank18_matrix[cbind(cmb[1L, ], cmb[2L, ])] <- rank18
    rank18_matrix[cbind(cmb[2L, ], cmb[1L, ])] <- rank18
    observed <- suppressWarnings(stats::cor(rank16, rank18, method = "pearson"))
    null <- apply(bp$matrix, 1L, function(p) {
      pm <- rank18_matrix[p, p, drop = FALSE]
      suppressWarnings(stats::cor(rank16, pm[lower.tri(pm)], method = "pearson"))
    })
    assert_true(all(is.finite(null)) && is.finite(observed), "Nonfinite Bray concordance statistic")
    p_two <- (1 + sum(abs(null) >= abs(observed) - TOLERANCE)) / (1 + length(null))
    bray_stats[[i]] <- data.frame(
      AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
      SampleN = length(ids), PairN = length(v16), SiteN = length(unique(meta$Site)),
      RhoObserved = observed, Alternative = "two.sided",
      PermutationRestriction = "sample labels permuted only within Site",
      PermutationMode = bp$mode, PermutationsRequested = N_PERMUTATIONS,
      PermutationsEvaluated = nrow(bp$matrix), TotalAdmissibleNonIdentity = bp$total,
      Seed = bray_seed, PPermutation = p_two, PAdjustedBHAcrossSubsets = NA_real_,
      InferenceScope = "cross-marker Bray association; overlapping subset surfaces; no causal or season claim",
      stringsAsFactors = FALSE)

    p16 <- vegan::wcmdscale(d16, k = 2L, eig = TRUE, add = "lingoes")$points
    p18 <- vegan::wcmdscale(d18, k = 2L, eig = TRUE, add = "lingoes")$points
    p16 <- as.matrix(p16); p18 <- as.matrix(p18)
    rownames(p16) <- ids; rownames(p18) <- ids
    proc_seed <- 94100L + i
    pp <- restricted_permutations(meta$Site, N_PERMUTATIONS, proc_seed)
    fit <- vegan::protest(p16, p18, permutations = pp$matrix, symmetric = TRUE)
    residual <- sqrt(rowSums((fit$X - fit$Yrot)^2))
    proc_scores[[i]] <- data.frame(
      AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
      SampleID = ids, Site = meta$Site, GeologicalHabitat = meta$GeologicalHabitat,
      SampleType = meta$SampleType, Campaign = meta$Campaign,
      Axis1_16S = fit$X[, 1L], Axis2_16S = fit$X[, 2L],
      Axis1_18S_Procrustes = fit$Yrot[, 1L], Axis2_18S_Procrustes = fit$Yrot[, 2L],
      Residual = residual, AxesUsed = 2L, OrdinationCorrection = "Lingoes",
      DistanceSurface = "Bray-Curtis on marker-specific fixed-depth rarefied counts",
      stringsAsFactors = FALSE)
    proc_stats[[i]] <- data.frame(
      AnalysisSubset = subset, SubsetDefinition = inputs$subsets$definitions[[subset]],
      SampleN = length(ids), SiteN = length(unique(meta$Site)), AxesUsed = 2L,
      ProcrustesM2 = fit$ss, ProcrustesCorrelation = fit$t0, Alternative = "greater",
      PermutationRestriction = "sample labels permuted only within Site",
      PermutationMode = pp$mode, PermutationsRequested = N_PERMUTATIONS,
      PermutationsEvaluated = nrow(pp$matrix), TotalAdmissibleNonIdentity = pp$total,
      Seed = proc_seed, PPermutation = fit$signif, PAdjustedBHAcrossSubsets = NA_real_,
      OrdinationCorrection = "Lingoes",
      DistanceSurface = "Bray-Curtis on marker-specific fixed-depth rarefied counts",
      InferenceScope = "two-axis symmetric Procrustes/protest; overlapping subsets; no causal or season claim",
      stringsAsFactors = FALSE)
  }
  bs <- do.call(rbind, bray_stats); ps <- do.call(rbind, proc_stats)
  bs$PAdjustedBHAcrossSubsets <- stats::p.adjust(bs$PPermutation, method = "BH")
  ps$PAdjustedBHAcrossSubsets <- stats::p.adjust(ps$PPermutation, method = "BH")
  list(pairwise = do.call(rbind, pair_parts), bray_stats = bs,
       procrustes_scores = do.call(rbind, proc_scores), procrustes_stats = ps)
}

write_tsv_new <- function(x, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA", fileEncoding = "UTF-8", eol = "\n")
}

write_text_new <- function(lines, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(lines), con = path, sep = "\n", useBytes = TRUE)
}

input_manifest <- function(inputs) data.frame(
  Role = names(inputs$paths),
  AbsolutePath = normalizePath(unname(inputs$paths), winslash = "/", mustWork = TRUE),
  SizeBytes = vapply(unname(inputs$paths), file_size, numeric(1L)),
  SHA256 = vapply(unname(inputs$paths), sha256_file, character(1L)),
  stringsAsFactors = FALSE)

schema_rows <- function(named_tables) {
  out <- list(); k <- 0L
  for (table_name in names(named_tables)) for (column_name in names(named_tables[[table_name]])) {
    k <- k + 1L
    out[[k]] <- data.frame(Table = table_name,
      ColumnOrder = match(column_name, names(named_tables[[table_name]])), Column = column_name,
      RStorageClass = paste(class(named_tables[[table_name]][[column_name]]), collapse = "|"),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

make_output_manifest <- function(root) {
  manifest_path <- file.path(root, "provenance", "output_sha256_manifest_v1.tsv")
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[normalizePath(files, winslash = "/", mustWork = TRUE) !=
                   normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  rel <- substring(normalizePath(files, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L)
  ord <- order(rel)
  data.frame(RelativePath = rel[ord], SizeBytes = vapply(files[ord], file_size, numeric(1L)),
             SHA256 = vapply(files[ord], sha256_file, character(1L)), stringsAsFactors = FALSE)
}

analysis_config <- function() data.frame(
  Key = c("script", "script_version", "markers", "subsets", "alpha_metrics",
          "alpha_inference", "distance", "ordination", "bray_permutation_alternative",
          "protest_alternative", "permutation_restriction", "permutations",
          "bootstrap_replicates", "campaign_labels", "pair_ids", "figure_backend"),
  Value = c(basename(script_path()), SCRIPT_VERSION, "16S|18S",
            paste(SUBSET_ORDER, collapse = "|"), paste(ALPHA_METRICS, collapse = "|"),
            "descriptive Spearman with Site-cluster percentile bootstrap CI",
            "Bray-Curtis on marker-specific fixed-depth rarefied counts",
            "two-axis Lingoes PCoA; symmetric Procrustes/protest",
            "two.sided", "greater", "within Site", as.character(N_PERMUTATIONS),
            as.character(N_BOOTSTRAPS), "First|Second", "not used", "none (numerical module only)"),
  stringsAsFactors = FALSE)

build_qc <- function(inputs, alpha, beta) {
  expected_pairs <- sum(vapply(EXPECTED_SUBSET_N, function(n) choose(n, 2), numeric(1L)))
  rows <- data.frame(
    Check = c("master_sample_count", "subset_count_vector", "alpha_pair_rows", "alpha_summary_rows",
              "bray_pair_rows", "bray_summary_rows", "procrustes_score_rows",
              "procrustes_summary_rows", "pair_ids_used", "campaign_labels",
              "permutation_restriction", "all_primary_values_finite"),
    Observed = c(nrow(inputs$master), paste(vapply(inputs$subsets$ids, length, integer(1L)), collapse = "|"),
      nrow(alpha$pairs), nrow(alpha$stats), nrow(beta$pairwise), nrow(beta$bray_stats),
      nrow(beta$procrustes_scores), nrow(beta$procrustes_stats), 0,
      paste(sort(unique(inputs$master$Campaign)), collapse = "|"),
      paste(unique(beta$bray_stats$PermutationRestriction), collapse = "|"),
      all(is.finite(c(alpha$pairs$Value16S, alpha$pairs$Value18S, beta$pairwise$Bray16S,
                      beta$pairwise$Bray18S, beta$procrustes_scores$Residual)))),
    Expected = c(75, paste(EXPECTED_SUBSET_N, collapse = "|"),
      sum(EXPECTED_SUBSET_N) * length(ALPHA_METRICS), length(SUBSET_ORDER) * length(ALPHA_METRICS),
      expected_pairs, length(SUBSET_ORDER), sum(EXPECTED_SUBSET_N), length(SUBSET_ORDER),
      0, "First|Second", "sample labels permuted only within Site", TRUE),
    stringsAsFactors = FALSE)
  rows$Status <- ifelse(rows$Observed == rows$Expected, "PASS", "FAIL")
  rows
}

readme_lines <- function() c(
  "# Cross-marker concordance package v1", "",
  "This immutable numerical package compares the protected 16S and 18S rarefied datasets across six prespecified, overlapping sample subsets.", "",
  "- Alpha: sample-wise Observed, Shannon, and Simpson associations are descriptive Spearman correlations with percentile intervals from resampling whole Sites.",
  "- Beta: marker-specific Bray-Curtis pairwise vectors are compared by Spearman correlation. Sample labels are permuted only within Site (9,999 deterministic permutations unless exhaustive enumeration is smaller).",
  "- Ordination: two-axis Lingoes-corrected PCoA configurations are compared by symmetric Procrustes/protest with the same within-Site restriction.",
  "- A72 membership is read from the canonical design table; replicate suffixes are never treated as PairIDs.",
  "- Campaign is retained as First/Second. Outputs do not support causal or season claims.", "",
  "Run `Rscript scripts/run_cross_marker_concordance_v1.R --check` before downstream plotting or reporting."
)

write_package <- function(stage, inputs, alpha, beta) {
  tables <- list(
    concordance_alpha_pairs_v1 = alpha$pairs,
    concordance_alpha_spearman_v1 = alpha$stats,
    concordance_bray_pairwise_v1 = beta$pairwise,
    concordance_bray_spearman_v1 = beta$bray_stats,
    concordance_procrustes_scores_v1 = beta$procrustes_scores,
    concordance_procrustes_summary_v1 = beta$procrustes_stats)
  for (nm in names(tables)) write_tsv_new(tables[[nm]], file.path(stage, "results", paste0(nm, ".tsv")))
  qc <- build_qc(inputs, alpha, beta)
  assert_true(all(qc$Status == "PASS"), "Concordance QC failed before finalization")
  write_text_new(readme_lines(), file.path(stage, "README_v1.md"))
  write_tsv_new(input_manifest(inputs), file.path(stage, "provenance", "input_sha256_manifest_v1.tsv"))
  write_tsv_new(analysis_config(), file.path(stage, "provenance", "analysis_config_v1.tsv"))
  write_text_new(capture.output(utils::sessionInfo()), file.path(stage, "provenance", "r_session_info_v1.txt"))
  write_tsv_new(qc, file.path(stage, "qc", "concordance_qc_summary_v1.tsv"))
  tables$concordance_qc_summary_v1 <- qc
  write_tsv_new(schema_rows(tables), file.path(stage, "provenance", "output_schema_v1.tsv"))
  write_tsv_new(make_output_manifest(stage), file.path(stage, "provenance", "output_sha256_manifest_v1.tsv"))
}

safe_remove_stage <- function(path, base) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  base_norm <- normalizePath(base, winslash = "/", mustWork = TRUE)
  allowed <- paste0(base_norm, "/.", TARGET_NAME, ".staging_")
  if (dir.exists(path) && startsWith(path_norm, allowed)) unlink(path, recursive = TRUE, force = TRUE)
}

run_analysis <- function(inputs, base) {
  final <- file.path(base, TARGET_NAME)
  if (dir.exists(final) || file.exists(final)) stopf("Refusing to overwrite existing output: %s", final)
  stage <- file.path(base, paste0(".", TARGET_NAME, ".staging_", Sys.getpid(), "_",
                                  format(Sys.time(), "%Y%m%dT%H%M%S")))
  assert_true(!dir.exists(stage) && !file.exists(stage), "Staging path collision")
  dir.create(stage, recursive = FALSE, showWarnings = FALSE)
  success <- FALSE
  on.exit(if (!success) safe_remove_stage(stage, base), add = TRUE)
  alpha <- compute_alpha(inputs)
  beta <- compute_beta_concordance(inputs)
  write_package(stage, inputs, alpha, beta)
  if (!file.rename(stage, final)) stopf("Could not commit staging directory: %s", final)
  success <- TRUE
  cat(sprintf("Created %s\n", final))
}

verify_output_manifest <- function(root) {
  path <- file.path(root, "provenance", "output_sha256_manifest_v1.tsv")
  manifest <- read_tsv(path)
  assert_schema(manifest, c("RelativePath", "SizeBytes", "SHA256"), "output manifest")
  assert_unique(manifest$RelativePath, "output manifest RelativePath")
  for (i in seq_len(nrow(manifest))) {
    target <- file.path(root, manifest$RelativePath[[i]])
    assert_true(file.exists(target), sprintf("Manifest target missing: %s", target))
    assert_true(file_size(target) == as.numeric(manifest$SizeBytes[[i]]), sprintf("Size mismatch: %s", target))
    assert_true(identical(sha256_file(target), tolower(manifest$SHA256[[i]])), sprintf("SHA256 mismatch: %s", target))
  }
  actual <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  actual <- actual[file.info(actual)$isdir %in% FALSE]
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  actual_rel <- substring(normalizePath(actual, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L)
  expected <- c(manifest$RelativePath, "provenance/output_sha256_manifest_v1.tsv")
  assert_true(setequal(actual_rel, expected), "Output manifest coverage drift")
}

verify_input_manifest <- function(root) {
  x <- read_tsv(file.path(root, "provenance", "input_sha256_manifest_v1.tsv"))
  assert_schema(x, c("Role", "AbsolutePath", "SizeBytes", "SHA256"), "input manifest")
  for (i in seq_len(nrow(x))) {
    path <- x$AbsolutePath[[i]]
    assert_true(file.exists(path), sprintf("Current input missing: %s", path))
    assert_true(file_size(path) == as.numeric(x$SizeBytes[[i]]), sprintf("Current input size drift: %s", path))
    assert_true(identical(sha256_file(path), tolower(x$SHA256[[i]])), sprintf("Current input SHA256 drift: %s", path))
  }
}

check_existing <- function(base) {
  root <- file.path(base, TARGET_NAME)
  assert_true(dir.exists(root), sprintf("Output package is missing: %s", root))
  verify_output_manifest(root); verify_input_manifest(root)
  files <- c("concordance_alpha_pairs_v1.tsv", "concordance_alpha_spearman_v1.tsv",
             "concordance_bray_pairwise_v1.tsv", "concordance_bray_spearman_v1.tsv",
             "concordance_procrustes_scores_v1.tsv", "concordance_procrustes_summary_v1.tsv")
  x <- lapply(setNames(files, files), function(f) read_tsv(file.path(root, "results", f)))
  expected_rows <- c(sum(EXPECTED_SUBSET_N) * length(ALPHA_METRICS),
    length(SUBSET_ORDER) * length(ALPHA_METRICS),
    sum(vapply(EXPECTED_SUBSET_N, function(n) choose(n, 2), numeric(1L))),
    length(SUBSET_ORDER), sum(EXPECTED_SUBSET_N), length(SUBSET_ORDER))
  assert_true(identical(as.integer(vapply(x, nrow, integer(1L))), as.integer(expected_rows)),
              "Concordance output row-count drift")
  required <- list(
    c("AnalysisSubset", "SubsetDefinition", "SampleID", "Site", "GeologicalHabitat", "SampleType", "Campaign", "Metric", "Value16S", "Value18S", "CountSurface16S", "CountSurface18S", "InferenceScope"),
    c("AnalysisSubset", "SubsetDefinition", "Metric", "SampleN", "SiteN", "RhoObserved", "PAsymptoticDescriptive", "BootstrapUnit", "BootstrapReplicatesRequested", "BootstrapReplicatesValid", "BootstrapSeed", "BootstrapCILevel", "BootstrapCILower", "BootstrapCIUpper", "InferenceScope"),
    c("AnalysisSubset", "SubsetDefinition", "SampleID1", "SampleID2", "Site1", "Site2", "GeologicalHabitat1", "GeologicalHabitat2", "SampleType1", "SampleType2", "Campaign1", "Campaign2", "PairScopeSite", "Bray16S", "Bray18S", "DistanceSurface16S", "DistanceSurface18S"),
    c("AnalysisSubset", "SubsetDefinition", "SampleN", "PairN", "SiteN", "RhoObserved", "Alternative", "PermutationRestriction", "PermutationMode", "PermutationsRequested", "PermutationsEvaluated", "TotalAdmissibleNonIdentity", "Seed", "PPermutation", "PAdjustedBHAcrossSubsets", "InferenceScope"),
    c("AnalysisSubset", "SubsetDefinition", "SampleID", "Site", "GeologicalHabitat", "SampleType", "Campaign", "Axis1_16S", "Axis2_16S", "Axis1_18S_Procrustes", "Axis2_18S_Procrustes", "Residual", "AxesUsed", "OrdinationCorrection", "DistanceSurface"),
    c("AnalysisSubset", "SubsetDefinition", "SampleN", "SiteN", "AxesUsed", "ProcrustesM2", "ProcrustesCorrelation", "Alternative", "PermutationRestriction", "PermutationMode", "PermutationsRequested", "PermutationsEvaluated", "TotalAdmissibleNonIdentity", "Seed", "PPermutation", "PAdjustedBHAcrossSubsets", "OrdinationCorrection", "DistanceSurface", "InferenceScope"))
  for (i in seq_along(x)) {
    assert_true(identical(names(x[[i]]), required[[i]]), sprintf("Exact schema drift in %s", files[[i]]))
    assert_true(setequal(unique(x[[i]]$AnalysisSubset), SUBSET_ORDER), sprintf("Subset drift in %s", files[[i]]))
  }
  assert_true(!any(grepl("Autumn|Season", unlist(x), ignore.case = FALSE), na.rm = TRUE),
              "Prohibited Autumn/Season label found in concordance outputs")
  qc <- read_tsv(file.path(root, "qc", "concordance_qc_summary_v1.tsv"))
  assert_true(all(qc$Status == "PASS"), "Stored concordance QC contains non-PASS row")
  cat(sprintf("PASS: %s\n", root))
}

dry_run_report <- function(inputs, base) {
  target <- file.path(base, TARGET_NAME)
  cat("DRY-RUN PASS: canonical design, alpha, beta v2, rarefied counts, hashes, schemas, and six sample subsets validated.\n")
  cat(sprintf("Subset counts: %s\n", paste(sprintf("%s=%d", names(EXPECTED_SUBSET_N), EXPECTED_SUBSET_N), collapse = ", ")))
  cat(sprintf("16S: %d samples x %d ASVs at 20,000 reads/sample\n", nrow(inputs$counts$`16S`), ncol(inputs$counts$`16S`)))
  cat(sprintf("18S: %d samples x %d ASVs at 50,000 reads/sample\n", nrow(inputs$counts$`18S`), ncol(inputs$counts$`18S`)))
  cat(sprintf("Planned immutable output: %s\n", target))
  if (dir.exists(target) || file.exists(target)) cat("NOTE: --run would refuse because the final output exists.\n")
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) { usage(); return(invisible(0L)) }
  check_packages()
  if (cli$mode == "check") check_existing(cli$base) else {
    inputs <- load_inputs(cli$base)
    if (cli$mode == "dry-run") dry_run_report(inputs, cli$base) else run_analysis(inputs, cli$base)
  }
  invisible(0L)
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L, runLast = FALSE)
})
