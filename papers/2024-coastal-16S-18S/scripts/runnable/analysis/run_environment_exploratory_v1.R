#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(digest)
})

SCRIPT_VERSION <- "1.0.0"
TARGET_NAME <- "08_environment_exploratory_v1"
SCOPE <- list(
  SetID = "A7_SITE_CENTROIDS",
  Campaign = "First",
  SeasonDisplay = "Spring",
  SampleType = "Sea-sand",
  InferenceScope = "exploratory_site_level_n7"
)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert_true <- function(x, msg) if (!isTRUE(x)) stopf("%s", msg)

script_path <- function() {
  x <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(x)) stopf("Cannot resolve script path")
  normalizePath(sub("^--file=", "", x[[1L]]), winslash = "/", mustWork = TRUE)
}

usage <- function() {
  cat("Usage: Rscript run_environment_exploratory_v1.R --dry-run|--run|--check [--base PATH]\n")
}

parse_cli <- function(args) {
  if (any(args %in% c("-h", "--help"))) return(list(help = TRUE))
  modes <- intersect(args, c("--dry-run", "--run", "--check"))
  assert_true(length(modes) == 1L, "Choose exactly one mode: --dry-run, --run, or --check")
  base <- normalizePath(file.path(dirname(script_path()), ".."), winslash = "/", mustWork = TRUE)
  hit <- which(args == "--base")
  if (length(hit)) {
    assert_true(length(hit) == 1L && hit < length(args), "--base requires one path")
    base <- normalizePath(args[[hit + 1L]], winslash = "/", mustWork = TRUE)
  }
  list(help = FALSE, mode = sub("^--", "", modes), base = base)
}

read_tsv <- function(path) {
  assert_true(file.exists(path), sprintf("Required input missing: %s", path))
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
}

write_tsv_new <- function(x, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
              na = "NA", eol = "\n", fileEncoding = "UTF-8")
}

write_lines_new <- function(x, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, path, useBytes = TRUE)
}

sha256_file <- function(path) tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
file_size <- function(path) unname(file.info(path)$size)

input_paths <- function(base) list(
  master_metadata = file.path(base, "00_design_v2", "master_sample_metadata_v2.tsv"),
  environment_link = file.path(base, "00_design_v2", "environment_site_metadata_link_v2.tsv"),
  counts_16s = file.path(base, "04_input_16s_rarefied_20000_v1", "results", "16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  taxonomy_16s = file.path(base, "04_input_16s_rarefied_20000_v1", "results", "16s_taxonomy_rarefied_20000_analysis_v1.tsv"),
  counts_18s = file.path(base, "03_18s_rarefied_50000_v1", "results", "18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
  taxonomy_18s = file.path(base, "03_18s_rarefied_50000_v1", "results", "18s_taxonomy_pr2_primary_rarefied_50000_analysis_v1.tsv"),
  environment_csv = file.path(base, "ori_metadata_env", "env.csv"),
  runner = script_path()
)

read_count <- function(path, expected_depth) {
  x <- read_tsv(path)
  assert_true(ncol(x) == 76L, sprintf("Expected ID + 75 samples in %s", path))
  ids <- as.character(x[[1L]])
  assert_true(!anyDuplicated(ids) && all(nzchar(ids)), "Feature IDs must be unique and nonblank")
  m <- t(as.matrix(x[-1L]))
  storage.mode(m) <- "numeric"
  rownames(m) <- names(x)[-1L]
  colnames(m) <- ids
  assert_true(all(is.finite(m)) && all(m >= 0) && all(m == round(m)), "Count table contains invalid cells")
  assert_true(all(rowSums(m) == expected_depth), sprintf("Rarefied sample depth differs from %d", expected_depth))
  m
}

validate_inputs <- function(base) {
  paths <- input_paths(base)
  missing <- paths[!file.exists(unlist(paths))]
  assert_true(!length(missing), sprintf("Missing input(s): %s", paste(names(missing), collapse = ", ")))
  master <- read_tsv(paths$master_metadata)
  env <- read_tsv(paths$environment_link)
  assert_true(nrow(master) == 75L && !anyDuplicated(master$SampleID), "Master metadata must contain 75 unique samples")
  eligible <- as.character(master$EnvironmentEligible) == "TRUE"
  assert_true(sum(eligible) == 21L, "Environment eligibility must contain exactly 21 samples")
  assert_true(all(master$Campaign[eligible] == "First") && all(master$SampleType[eligible] == "Sea-sand"),
              "Environment-eligible samples must be First Sea-sand")
  assert_true(nrow(env) == 7L && !anyDuplicated(env$Site), "Environment link must contain seven unique sites")
  assert_true(all(env$EnvironmentCampaignScope == "First") && all(env$EnvironmentSeasonLabel == "Spring") &&
                all(env$EnvironmentSampleTypeScope == "Sea-sand"), "Environment link scope drift")
  assert_true(as.numeric(env$TOC[env$SiteCodeOriginal == "SH"]) == 0,
              "Expected legacy Sehwa TOC=0 representation for documented ND caveat")
  c16 <- read_count(paths$counts_16s, 20000)
  c18 <- read_count(paths$counts_18s, 50000)
  assert_true(setequal(rownames(c16), master$SampleID), "16S sample set differs from master metadata")
  assert_true(setequal(rownames(c18), master$SampleID), "18S sample set differs from master metadata")
  c16 <- c16[master$SampleID, , drop = FALSE]
  c18 <- c18[master$SampleID, , drop = FALSE]
  list(paths = paths, master = master, env = env, counts16 = c16, counts18 = c18)
}

scope_columns <- function(n) data.frame(
  SetID = rep(SCOPE$SetID, n), Campaign = rep(SCOPE$Campaign, n),
  SeasonDisplay = rep(SCOPE$SeasonDisplay, n), SampleType = rep(SCOPE$SampleType, n),
  InferenceScope = rep(SCOPE$InferenceScope, n), stringsAsFactors = FALSE
)

add_scope <- function(x) cbind(x, scope_columns(nrow(x)), stringsAsFactors = FALSE)

site_order <- function(env) as.character(env$Site)

environment_matrix <- function(env, include_toc = FALSE, sh_toc = NULL) {
  vars <- c("pH", "NaCl", "T-N", "T-S", "Ca", "Mg", "K", "Na", "T-P", "water")
  raw <- as.data.frame(lapply(vars, function(v) as.numeric(env[[v]])), check.names = FALSE)
  names(raw) <- vars
  if (include_toc) {
    toc <- as.numeric(env$TOC)
    if (!is.null(sh_toc)) toc[env$SiteCodeOriginal == "SH"] <- sh_toc
    raw$TOC <- toc
  }
  log_vars <- c("NaCl", "T-N", "T-S", "Ca", "Mg", "K", "Na", "T-P")
  transformed <- raw
  transformed[log_vars] <- lapply(transformed[log_vars], log10)
  scaled <- scale(transformed)
  rownames(scaled) <- env$Site
  list(raw = raw, transformed = transformed, scaled = scaled, variables = names(raw))
}

aggregate_site_ra <- function(counts, master, env) {
  eligible_ids <- master$SampleID[as.character(master$EnvironmentEligible) == "TRUE"]
  x <- counts[eligible_ids, , drop = FALSE]
  x <- x / rowSums(x)
  out <- do.call(rbind, lapply(env$Site, function(s) {
    ids <- master$SampleID[as.character(master$EnvironmentEligible) == "TRUE" & master$Site == s]
    assert_true(length(ids) == 3L, sprintf("Expected three eligible observations for %s", s))
    colMeans(x[ids, , drop = FALSE])
  }))
  rownames(out) <- env$Site
  colnames(out) <- colnames(x)
  assert_true(max(abs(rowSums(out) - 1)) < 1e-10, "Site-mean relative abundance does not sum to one")
  out
}

aggregate_taxon <- function(site_asv, taxonomy_path, id_col, rank_col) {
  tax <- read_tsv(taxonomy_path)
  assert_true(all(c(id_col, rank_col) %in% names(tax)), "Taxonomy schema mismatch")
  idx <- match(colnames(site_asv), tax[[id_col]])
  assert_true(!anyNA(idx), "Count features missing from taxonomy")
  label <- trimws(as.character(tax[[rank_col]][idx]))
  label[is.na(label) | !nzchar(label)] <- "Unassigned"
  t(rowsum(t(site_asv), group = label, reorder = FALSE))
}

all_permutations <- function(n) {
  p <- permute::allPerms(n, control = permute::how(within = permute::Within(type = "free")))
  assert_true(nrow(p) == factorial(n) - 1L, sprintf("Permutation count drift for n=%d", n))
  p
}

extract_anova_row <- function(x, row, marker, level, type, term, r2, adjr2) {
  d <- as.data.frame(x)
  nms <- names(d)
  variance_col <- intersect(c("Variance", "SumOfSqs"), nms)
  f_col <- intersect(c("F", "F.Model"), nms)
  p_col <- grep("Pr\\(>F\\)", nms, value = TRUE)
  data.frame(
    Marker = marker, FeatureLevel = level, AnalysisSurface = "rarefied_primary",
    TestType = type, Term = term,
    Df = as.numeric(d[row, "Df"]),
    Variance = if (length(variance_col)) as.numeric(d[row, variance_col[[1L]]]) else NA_real_,
    F = if (length(f_col)) as.numeric(d[row, f_col[[1L]]]) else NA_real_,
    P = if (length(p_col)) as.numeric(d[row, p_col[[1L]]]) else NA_real_,
    Permutations = 5039L, R2 = r2, AdjustedR2 = adjr2,
    stringsAsFactors = FALSE
  )
}

run_rda <- function(marker, feature_level, site_taxon, pc_scores, env, permutations) {
  y <- vegan::decostand(site_taxon, method = "hellinger")
  dat <- data.frame(EnvPC1 = pc_scores[, 1], EnvPC2 = pc_scores[, 2], row.names = rownames(pc_scores))
  model <- vegan::rda(y ~ EnvPC1 + EnvPC2, data = dat)
  rs <- vegan::RsquareAdj(model)
  overall <- anova(model, permutations = permutations)
  terms <- anova(model, by = "term", permutations = permutations)
  tests <- rbind(
    extract_anova_row(overall, 1L, marker, feature_level, "overall", "Model", rs$r.squared, rs$adj.r.squared),
    extract_anova_row(terms, 1L, marker, feature_level, "term", "EnvPC1", rs$r.squared, rs$adj.r.squared),
    extract_anova_row(terms, 2L, marker, feature_level, "term", "EnvPC2", rs$r.squared, rs$adj.r.squared)
  )
  ss <- scores(model, display = "sites", choices = 1:2, scaling = 2)
  assert_true(ncol(ss) >= 2L, sprintf("%s RDA has fewer than two constrained axes", marker))
  eig <- model$CCA$eig
  pct <- 100 * eig / sum(eig)
  site_rows <- data.frame(
    Marker = marker, FeatureLevel = feature_level, AnalysisSurface = "rarefied_primary",
    SiteCodeOriginal = env$SiteCodeOriginal, Site = env$Site, GeologicalHabitat = env$GeologicalHabitat,
    RDA1 = ss[env$Site, 1], RDA2 = ss[env$Site, 2],
    RDA1PercentConstrained = pct[[1L]], RDA2PercentConstrained = pct[[2L]], stringsAsFactors = FALSE
  )
  bp <- scores(model, display = "bp", choices = 1:2, scaling = 2)
  biplot_rows <- data.frame(
    Marker = marker, FeatureLevel = feature_level, ScoreType = "environment_PC_predictor",
    Label = rownames(bp), RDA1 = bp[, 1], RDA2 = bp[, 2],
    MeanRelativeAbundance = NA_real_, AbundanceRank = NA_integer_, EnvfitR2 = NA_real_, EnvfitP = NA_real_,
    EnvfitPAdjustedBH = NA_real_, ArrowSelected = TRUE, stringsAsFactors = FALSE
  )
  list(model = model, tests = add_scope(tests), sites = add_scope(site_rows), biplot = add_scope(biplot_rows))
}

mantel_one <- function(site_asv, env_scaled, method = "spearman", permutations = NULL) {
  cd <- vegan::vegdist(site_asv, method = "bray")
  ed <- dist(env_scaled, method = "euclidean")
  vegan::mantel(cd, ed, method = method, permutations = permutations, na.rm = TRUE)
}

mantel_outputs <- function(marker, site_asv, env_scaled, env, perm7) {
  full <- mantel_one(site_asv, env_scaled, "spearman", perm7)
  summary <- data.frame(
    Marker = marker, Method = "Spearman", AnalysisSurface = "ASV_site_mean_relative_abundance_Bray_vs_environment_scaled_Euclidean",
    Nsites = 7L, Statistic = unname(full$statistic), P = unname(full$signif), Permutations = 5039L,
    PAdjustedBHPrimaryAcrossMarkers = NA_real_, stringsAsFactors = FALSE
  )
  loso <- do.call(rbind, lapply(seq_len(nrow(env)), function(i) {
    keep <- setdiff(seq_len(nrow(env)), i)
    p6 <- all_permutations(6L)
    fit <- mantel_one(site_asv[keep, , drop = FALSE], env_scaled[keep, , drop = FALSE], "spearman", p6)
    stat <- unname(fit$statistic)
    data.frame(
      Marker = marker, Method = "Spearman", LeftOutSiteCode = env$SiteCodeOriginal[[i]], LeftOutSite = env$Site[[i]],
      Nsites = 6L, Statistic = stat, P = unname(fit$signif), Permutations = 719L,
      Sign = ifelse(stat > 0, "positive", ifelse(stat < 0, "negative", "zero")),
      FullStatistic = unname(full$statistic), SignMatchesFull = sign(stat) == sign(unname(full$statistic)),
      GeologicalHabitat = env$GeologicalHabitat[[i]], stringsAsFactors = FALSE
    )
  }))
  pearson <- mantel_one(site_asv, env_scaled, "pearson", perm7)
  pearson_row <- data.frame(
    Marker = marker, Method = "Pearson_sensitivity", Nsites = 7L,
    Statistic = unname(pearson$statistic), P = unname(pearson$signif), Permutations = 5039L,
    stringsAsFactors = FALSE
  )
  list(summary = add_scope(summary), loso = add_scope(loso), pearson = add_scope(pearson_row))
}

build_results <- function(inputs) {
  env <- inputs$env
  em <- environment_matrix(env, include_toc = FALSE)
  pca <- prcomp(em$scaled, center = FALSE, scale. = FALSE)
  prop <- pca$sdev^2 / sum(pca$sdev^2)
  station_label <- ifelse(env$SiteCodeOriginal == "DD", "Dokdo-1", env$SiteCodeOriginal)
  station_status <- ifelse(
    env$SiteCodeOriginal == "DD",
    "environment_station_Dokdo-1_to_microbiome_site_exact_link_not_documented",
    "site_code_link_from_original_environment_table"
  )
  pca_sites <- data.frame(
    SiteCodeOriginal = env$SiteCodeOriginal, Site = env$Site, GeologicalHabitat = env$GeologicalHabitat,
    EnvironmentStationLabel = station_label, StationMappingStatus = station_status,
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    PC1ProportionVariance = prop[[1L]], PC2ProportionVariance = prop[[2L]], stringsAsFactors = FALSE
  )
  units <- c(pH = "unitless", NaCl = "%", `T-N` = "%", `T-S` = "mg/kg", Ca = "mg/kg", Mg = "mg/kg",
             K = "mg/kg", Na = "mg/kg", `T-P` = "mg/kg", water = "%")
  transforms <- ifelse(em$variables %in% c("pH", "water"), "z-score", "log10_then_z-score")
  pca_loadings <- data.frame(
    Variable = em$variables, Unit = unname(units[em$variables]), Transform = transforms,
    PC1 = pca$rotation[em$variables, 1], PC2 = pca$rotation[em$variables, 2],
    PC1ProportionVariance = prop[[1L]], PC2ProportionVariance = prop[[2L]], stringsAsFactors = FALSE
  )

  site16 <- aggregate_site_ra(inputs$counts16, inputs$master, env)
  site18 <- aggregate_site_ra(inputs$counts18, inputs$master, env)
  tax16 <- aggregate_taxon(site16, inputs$paths$taxonomy_16s, "OTUID", "Genus")
  tax18 <- aggregate_taxon(site18, inputs$paths$taxonomy_18s, "ASV_ID", "Class")
  perm7 <- all_permutations(7L)
  rda16 <- run_rda("16S", "Genus", tax16, pca$x[, 1:2, drop = FALSE], env, perm7)
  rda18 <- run_rda("18S", "Class", tax18, pca$x[, 1:2, drop = FALSE], env, perm7)
  rda_sites <- rbind(rda16$sites, rda18$sites)
  rda_biplot <- rbind(rda16$biplot, rda18$biplot)
  rda_tests <- rbind(rda16$tests, rda18$tests)

  m16 <- mantel_outputs("16S", site16, em$scaled, env, perm7)
  m18 <- mantel_outputs("18S", site18, em$scaled, env, perm7)
  mantel <- rbind(m16$summary, m18$summary)
  mantel$PAdjustedBHPrimaryAcrossMarkers <- p.adjust(mantel$P, method = "BH")
  mantel_loso <- rbind(m16$loso, m18$loso)
  pearson <- rbind(m16$pearson, m18$pearson)

  rda_asv <- do.call(rbind, lapply(list(`16S` = site16, `18S` = site18), function(x) {
    marker <- if (identical(x, site16)) "16S" else "18S"
    model <- vegan::rda(vegan::decostand(x, "hellinger") ~ EnvPC1 + EnvPC2,
                        data = data.frame(EnvPC1 = pca$x[, 1], EnvPC2 = pca$x[, 2], row.names = env$Site))
    a <- anova(model, permutations = perm7)
    rs <- vegan::RsquareAdj(model)
    data.frame(Marker = marker, FeatureLevel = "ASV", TestType = "overall", Term = "Model",
               F = as.numeric(as.data.frame(a)[1, "F"]), P = as.numeric(as.data.frame(a)[1, "Pr(>F)"]),
               Permutations = 5039L, R2 = rs$r.squared, AdjustedR2 = rs$adj.r.squared,
               AnalysisSurface = "secondary_ASV_sensitivity", stringsAsFactors = FALSE)
  }))
  rda_asv <- add_scope(rda_asv)

  toc_sensitivity <- do.call(rbind, lapply(c(0, 0.06, 0.12), function(toc_value) {
    et <- environment_matrix(env, include_toc = TRUE, sh_toc = toc_value)
    do.call(rbind, lapply(list(`16S` = site16, `18S` = site18), function(x) {
      marker <- if (identical(x, site16)) "16S" else "18S"
      fit <- mantel_one(x, et$scaled, "spearman", perm7)
      data.frame(Marker = marker, SehwaTOCSubstitutionPercent = toc_value, Method = "Spearman",
                 Nsites = 7L, Statistic = unname(fit$statistic), P = unname(fit$signif), Permutations = 5039L,
                 TOCTransform = "z-score_raw_TOC_sensitivity_only", stringsAsFactors = FALSE)
    }))
  }))
  toc_sensitivity <- add_scope(toc_sensitivity)

  dates <- c(PS = "2024-03-07", SS = "2024-03-07", SH = "2024-03-07", GS = "2024-03-12",
             DD = "2024-05-01", ULD = "2024-05-02", PH = "2024-05-03")
  source_ledger <- cbind(
    env[, c("EnvironmentRecordID", "SiteCodeOriginal", "Site", "GeologicalHabitat")],
    SamplingDate = unname(dates[env$SiteCodeOriginal]), EnvironmentStationLabel = station_label,
    StationMappingStatus = station_status,
    PrimaryVariables = "pH|NaCl|T-N|T-S|Ca|Mg|K|Na|T-P|water",
    PrimaryTOCPolicy = "excluded_Sehwa_ND_legacy_zero_without_LOD_LOQ",
    stringsAsFactors = FALSE
  )
  source_ledger <- add_scope(source_ledger)

  list(
    pca_sites = add_scope(pca_sites), pca_loadings = add_scope(pca_loadings),
    rda_sites = rda_sites, rda_biplot = rda_biplot, rda_tests = rda_tests,
    mantel = mantel, mantel_loso = mantel_loso, pearson = pearson,
    rda_asv = rda_asv, toc_sensitivity = toc_sensitivity, source_ledger = source_ledger
  )
}

validate_results <- function(x, env) {
  assert_true(nrow(x$pca_sites) == 7L && nrow(x$pca_loadings) == 10L, "PCA row-count contract failed")
  assert_true(!"TOC" %in% x$pca_loadings$Variable, "TOC leaked into primary PCA")
  assert_true(nrow(x$rda_sites) == 14L && nrow(x$rda_tests) == 6L, "RDA row-count contract failed")
  assert_true(nrow(x$mantel) == 2L && nrow(x$mantel_loso) == 14L, "Mantel row-count contract failed")
  assert_true(all(x$mantel$Nsites == 7L) && all(x$mantel_loso$Nsites == 6L), "Mantel independent-unit drift")
  required_scope <- c("SetID", "Campaign", "SeasonDisplay", "SampleType", "InferenceScope")
  for (item in x) {
    assert_true(all(required_scope %in% names(item)), "Scope columns missing")
    assert_true(all(item$SetID == SCOPE$SetID & item$Campaign == "First" & item$SeasonDisplay == "Spring" &
                      item$SampleType == "Sea-sand" & item$InferenceScope == SCOPE$InferenceScope), "Scope value drift")
    assert_true(!any(grepl("Autumn", unlist(item), fixed = TRUE), na.rm = TRUE), "Prohibited Autumn label found")
  }
  assert_true(setequal(x$pca_sites$Site, env$Site), "PCA site set drift")
  TRUE
}

table_schema <- function(named_tables) {
  do.call(rbind, lapply(names(named_tables), function(nm) {
    x <- named_tables[[nm]]
    data.frame(Table = nm, ColumnOrder = seq_along(x), Column = names(x),
               RStorageClass = vapply(x, function(z) class(z)[[1L]], character(1L)), stringsAsFactors = FALSE)
  }))
}

write_package <- function(results, inputs, stage) {
  dir.create(file.path(stage, "results"), recursive = TRUE)
  dir.create(file.path(stage, "qc"), recursive = TRUE)
  dir.create(file.path(stage, "provenance"), recursive = TRUE)
  result_tables <- list(
    environment_pca_site_scores_v1 = results$pca_sites,
    environment_pca_variable_loadings_v1 = results$pca_loadings,
    environment_rda_site_scores_v1 = results$rda_sites,
    environment_rda_biplot_scores_v1 = results$rda_biplot,
    environment_rda_model_tests_v1 = results$rda_tests,
    environment_mantel_summary_v1 = results$mantel,
    environment_mantel_loso_v1 = results$mantel_loso
  )
  for (nm in names(result_tables)) write_tsv_new(result_tables[[nm]], file.path(stage, "results", paste0(nm, ".tsv")))
  qc_tables <- list(
    environment_mantel_pearson_sensitivity_v1 = results$pearson,
    environment_rda_asv_sensitivity_v1 = results$rda_asv,
    environment_toc_sensitivity_v1 = results$toc_sensitivity
  )
  for (nm in names(qc_tables)) write_tsv_new(qc_tables[[nm]], file.path(stage, "qc", paste0(nm, ".tsv")))
  write_tsv_new(results$source_ledger, file.path(stage, "provenance", "environment_source_ledger_v1.tsv"))

  input_manifest <- do.call(rbind, lapply(names(inputs$paths), function(role) {
    p <- normalizePath(inputs$paths[[role]], winslash = "/", mustWork = TRUE)
    data.frame(Role = role, AbsolutePath = p, SizeBytes = file_size(p), SHA256 = sha256_file(p), stringsAsFactors = FALSE)
  }))
  write_tsv_new(input_manifest, file.path(stage, "provenance", "input_sha256_manifest_v1.tsv"))
  config <- data.frame(
    Key = c("script_version", "independent_unit", "sample_observations", "site_centroids", "primary_environment_variables",
            "primary_toc_policy", "community_aggregation", "rda_primary", "mantel_primary", "permutations", "claim_scope"),
    Value = c(SCRIPT_VERSION, "site", "21 First-campaign Spring Sea-sand samples", "7",
              "pH|NaCl|T-N|T-S|Ca|Mg|K|Na|T-P|water", "TOC excluded; Sehwa ND was legacy-encoded zero",
              "sample-wise relative abundance then arithmetic mean of three samples per site",
              "16S Genus and 18S PR2 Class Hellinger ~ EnvPC1 + EnvPC2",
              "ASV site-mean Bray vs 10-variable scaled Euclidean; Spearman", "5039 nonidentity site permutations",
              "exploratory site-level association n=7; no causal or habitat-replicated claim"), stringsAsFactors = FALSE
  )
  write_tsv_new(config, file.path(stage, "provenance", "analysis_config_v1.tsv"))
  schema <- table_schema(c(result_tables, qc_tables, list(environment_source_ledger_v1 = results$source_ledger)))
  write_tsv_new(schema, file.path(stage, "provenance", "output_schema_v1.tsv"))
  qc_summary <- data.frame(
    Check = c("eligible_sample_count", "site_count", "pca_primary_variable_count", "primary_toc_rows", "rda_site_rows", "mantel_full_rows", "mantel_loso_rows", "nonidentity_permutations"),
    Observed = c(21, 7, nrow(results$pca_loadings), sum(results$pca_loadings$Variable == "TOC"), nrow(results$rda_sites), nrow(results$mantel), nrow(results$mantel_loso), 5039),
    Expected = c(21, 7, 10, 0, 14, 2, 14, 5039), Status = "PASS", stringsAsFactors = FALSE
  )
  write_tsv_new(qc_summary, file.path(stage, "qc", "environment_qc_summary_v1.tsv"))
  write_lines_new(capture.output(sessionInfo()), file.path(stage, "provenance", "r_session_info_v1.txt"))
  readme <- c(
    "# Environment exploratory package v1", "",
    "This package is restricted to 21 First-campaign Spring Sea-sand sample observations aggregated to seven independent site centroids.",
    "Environmental values are never copied to Water or Second-campaign samples.", "",
    "Primary PCA uses ten variables and excludes TOC because Sehwa was ND in the laboratory workbook but legacy env.csv encoded zero without LOD/LOQ.",
    "RDA uses 16S Genus and 18S PR2 Class Hellinger composition constrained by environmental PC1 and PC2.",
    "Mantel uses ASV site-mean Bray-Curtis against scaled environmental Euclidean distance with all 5,039 nonidentity site permutations.",
    "Pearson, ASV-level RDA, and TOC-substitution results are sensitivity/QC outputs only.",
    "All inference is exploratory site-level n=7; causal driver and habitat-replicated claims are prohibited."
  )
  write_lines_new(readme, file.path(stage, "README_v1.md"))
  manifest_path <- file.path(stage, "provenance", "output_sha256_manifest_v1.tsv")
  files <- list.files(stage, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE & normalizePath(files, winslash = "/", mustWork = TRUE) != normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
  root_norm <- normalizePath(stage, winslash = "/", mustWork = TRUE)
  manifest <- data.frame(
    RelativePath = substring(normalizePath(files, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L),
    SizeBytes = vapply(files, file_size, numeric(1L)), SHA256 = vapply(files, sha256_file, character(1L)), stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$RelativePath), , drop = FALSE]
  write_tsv_new(manifest, manifest_path)
}

verify_manifest <- function(root) {
  path <- file.path(root, "provenance", "output_sha256_manifest_v1.tsv")
  m <- read_tsv(path)
  assert_true(identical(names(m), c("RelativePath", "SizeBytes", "SHA256")), "Output manifest schema drift")
  for (i in seq_len(nrow(m))) {
    p <- file.path(root, m$RelativePath[[i]])
    assert_true(file.exists(p) && file_size(p) == as.numeric(m$SizeBytes[[i]]) && sha256_file(p) == tolower(m$SHA256[[i]]),
                sprintf("Output manifest mismatch: %s", p))
  }
  actual <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  actual <- actual[file.info(actual)$isdir %in% FALSE]
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  rel <- substring(normalizePath(actual, winslash = "/", mustWork = TRUE), nchar(root_norm) + 2L)
  assert_true(setequal(rel, c(m$RelativePath, "provenance/output_sha256_manifest_v1.tsv")), "Output manifest coverage drift")
}

verify_input_manifest <- function(root) {
  m <- read_tsv(file.path(root, "provenance", "input_sha256_manifest_v1.tsv"))
  for (i in seq_len(nrow(m))) {
    p <- m$AbsolutePath[[i]]
    assert_true(file.exists(p) && file_size(p) == as.numeric(m$SizeBytes[[i]]) && sha256_file(p) == tolower(m$SHA256[[i]]),
                sprintf("Input drift: %s", p))
  }
}

check_existing <- function(base) {
  root <- file.path(base, TARGET_NAME)
  assert_true(dir.exists(root), sprintf("Output package missing: %s", root))
  verify_manifest(root); verify_input_manifest(root)
  x <- list(
    pca_sites = read_tsv(file.path(root, "results", "environment_pca_site_scores_v1.tsv")),
    pca_loadings = read_tsv(file.path(root, "results", "environment_pca_variable_loadings_v1.tsv")),
    rda_sites = read_tsv(file.path(root, "results", "environment_rda_site_scores_v1.tsv")),
    rda_biplot = read_tsv(file.path(root, "results", "environment_rda_biplot_scores_v1.tsv")),
    rda_tests = read_tsv(file.path(root, "results", "environment_rda_model_tests_v1.tsv")),
    mantel = read_tsv(file.path(root, "results", "environment_mantel_summary_v1.tsv")),
    mantel_loso = read_tsv(file.path(root, "results", "environment_mantel_loso_v1.tsv")),
    pearson = read_tsv(file.path(root, "qc", "environment_mantel_pearson_sensitivity_v1.tsv")),
    rda_asv = read_tsv(file.path(root, "qc", "environment_rda_asv_sensitivity_v1.tsv")),
    toc_sensitivity = read_tsv(file.path(root, "qc", "environment_toc_sensitivity_v1.tsv")),
    source_ledger = read_tsv(file.path(root, "provenance", "environment_source_ledger_v1.tsv"))
  )
  env <- read_tsv(file.path(base, "00_design_v2", "environment_site_metadata_link_v2.tsv"))
  validate_results(x, env)
  cat(sprintf("PASS: %s\n", root))
}

run_analysis <- function(inputs, base) {
  final <- file.path(base, TARGET_NAME)
  assert_true(!dir.exists(final) && !file.exists(final), sprintf("Refusing to overwrite existing output: %s", final))
  stage <- file.path(base, sprintf(".%s.stage_%d", TARGET_NAME, Sys.getpid()))
  assert_true(!dir.exists(stage) && !file.exists(stage), "Staging path collision")
  dir.create(stage)
  ok <- FALSE
  on.exit(if (!ok && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  results <- build_results(inputs)
  validate_results(results, inputs$env)
  write_package(results, inputs, stage)
  verify_manifest(stage)
  assert_true(file.rename(stage, final), "Atomic stage promotion failed")
  ok <- TRUE
  check_existing(base)
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) { usage(); return(invisible(0L)) }
  if (cli$mode == "check") return(check_existing(cli$base))
  inputs <- validate_inputs(cli$base)
  if (cli$mode == "dry-run") {
    cat(sprintf("DRY-RUN PASS: 21 eligible samples -> 7 sites; 16S=%d ASVs, 18S=%d ASVs; primary PCA=10 variables excluding TOC; target=%s\n",
                ncol(inputs$counts16), ncol(inputs$counts18), file.path(cli$base, TARGET_NAME)))
  } else run_analysis(inputs, cli$base)
  invisible(0L)
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L, runLast = FALSE)
})
