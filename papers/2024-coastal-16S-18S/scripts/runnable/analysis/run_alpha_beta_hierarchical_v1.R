#!/usr/bin/env Rscript

# Hierarchical alpha- and beta-diversity analysis for the protected 75-sample
# 16S/18S revision data.  This script is intentionally figure-free: plotting
# is handled by a separate R-only figure layer after the numerical outputs have
# passed --check.

options(stringsAsFactors = FALSE, warn = 1, digits = 15)

SCRIPT_VERSION <- "3.0.0"
EXPECTED_SET_COUNTS <- c(
  A75_FULL = 75L,
  A72_SUBSTRATE_PAIRED_CELLS = 72L,
  A60_COMPLETE_FACTORIAL = 60L,
  A36_SAND_SEASON = 36L,
  A30_WATER_SEASON = 30L,
  A42_FIRST_SEASON_SITE = 42L
)
SET_ORDER <- names(EXPECTED_SET_COUNTS)
METRIC_ORDER <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
MARKER_ORDER <- c("16S", "18S")
SURFACE_ORDER <- c("rarefied_primary", "unrarefied_tss_sensitivity")
N_PERMUTATIONS <- 9999L
BASE_SEED <- 54101L
TOLERANCE <- 1e-12

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}

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
    "run_alpha_beta_hierarchical_v1.R ", SCRIPT_VERSION, "\n\n",
    "Usage:\n",
    "  Rscript run_alpha_beta_hierarchical_v1.R --dry-run [--base PATH]\n",
    "  Rscript run_alpha_beta_hierarchical_v1.R --run     [--base PATH]\n",
    "  Rscript run_alpha_beta_hierarchical_v1.R --check   [--base PATH]\n\n",
    "Modes are mutually exclusive. --run refuses to overwrite either final\n",
    "output directory and builds both packages in sibling staging folders.\n",
    "--dry-run validates all inputs and schemas without writing any file.\n",
    "--check validates the existing packages, SHA256 manifests, core schemas,\n",
    "set membership, and current input hashes without modifying outputs.\n\n",
    "Outputs:\n",
    "  05_alpha_diversity_v3\n",
    "  06_beta_diversity_v4\n"
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
    } else {
      stopf("Unknown argument: %s", arg)
    }
  }
  assert_true(dir.exists(base), sprintf("Base directory does not exist: %s", base))
  list(help = FALSE, mode = sub("^--", "", modes[[1L]]),
       base = normalizePath(base, winslash = "/", mustWork = TRUE))
}

required_packages <- function() c("vegan", "permute", "digest")

check_packages <- function() {
  missing <- required_packages()[!vapply(required_packages(), requireNamespace,
                                          quietly = TRUE, FUN.VALUE = logical(1L))]
  if (length(missing)) stopf("Missing required R package(s): %s", paste(missing, collapse = ", "))
}

input_paths <- function(base) {
  c(
    master_metadata = file.path(base, "00_design_v4", "master_sample_metadata_v4.tsv"),
    set_membership = file.path(base, "00_design_v4", "analysis_set_membership_v4.tsv"),
    model_registry = file.path(base, "00_design_v4", "model_registry_v4.tsv"),
    counts_16s_rarefied = file.path(base, "04_input_16s_rarefied_20000_v1", "results",
                                    "16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
    counts_18s_rarefied = file.path(base, "03_18s_rarefied_50000_v1", "results",
                                    "18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
    counts_16s_unrarefied = file.path(base, "01_input_16s_75_v1",
                                      "16s_asv_counts_analysis_75_v1.tsv"),
    counts_18s_unrarefied = file.path(base, "02_18s_unoise_pr2_v1", "results",
                                      "18s_asv_counts_pr2_analysis_v1.tsv")
  )
}

output_paths <- function(base) {
  c(alpha = file.path(base, "05_alpha_diversity_v3"),
    beta = file.path(base, "06_beta_diversity_v4"))
}

read_tsv <- function(path) {
  assert_true(file.exists(path), sprintf("Required file is missing: %s", path))
  tryCatch(
    utils::read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                      check.names = FALSE, stringsAsFactors = FALSE,
                      na.strings = c("NA"), fileEncoding = "UTF-8-BOM"),
    error = function(e) stopf("Cannot read %s: %s", path, conditionMessage(e))
  )
}

assert_schema <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stopf("%s is missing column(s): %s", label, paste(missing, collapse = ", "))
}

assert_unique <- function(x, label) {
  if (anyNA(x) || any(!nzchar(x))) stopf("%s contains blank/NA identifiers", label)
  duplicates <- unique(x[duplicated(x)])
  if (length(duplicates)) stopf("%s contains duplicate identifier(s): %s",
                                label, paste(head(duplicates, 10L), collapse = ", "))
}

read_count_table <- function(path, expected_samples, marker, surface) {
  table <- read_tsv(path)
  assert_true(ncol(table) >= 2L, sprintf("%s %s count table has fewer than two columns", marker, surface))
  id_header <- names(table)[[1L]]
  assert_true(id_header %in% c("#OTUID", "OTUID", "#OTU ID", "ASVID", "FeatureID"),
              sprintf("Unexpected feature-ID header in %s: %s", path, id_header))
  feature_ids <- as.character(table[[1L]])
  assert_unique(feature_ids, sprintf("%s %s feature IDs", marker, surface))
  observed_samples <- names(table)[-1L]
  assert_unique(observed_samples, sprintf("%s %s sample columns", marker, surface))
  if (!setequal(observed_samples, expected_samples)) {
    stopf("%s %s sample columns differ from master metadata; missing=%s; extra=%s",
          marker, surface,
          paste(setdiff(expected_samples, observed_samples), collapse = "|"),
          paste(setdiff(observed_samples, expected_samples), collapse = "|"))
  }
  values <- as.matrix(table[, match(expected_samples, observed_samples) + 1L, drop = FALSE])
  suppressWarnings(storage.mode(values) <- "numeric")
  if (anyNA(values) || any(!is.finite(values))) stopf("%s %s count table has nonnumeric/nonfinite values", marker, surface)
  if (any(values < 0) || any(abs(values - round(values)) > TOLERANCE)) {
    stopf("%s %s count table must contain nonnegative integer counts", marker, surface)
  }
  values <- round(values)
  rownames(values) <- feature_ids
  colnames(values) <- expected_samples
  sample_by_feature <- t(values)
  depths <- rowSums(sample_by_feature)
  if (any(depths <= 0)) stopf("%s %s contains an empty sample", marker, surface)
  list(marker = marker, surface = surface, path = path, id_header = id_header,
       matrix = sample_by_feature, feature_count = ncol(sample_by_feature), depths = depths)
}

validate_design_inputs <- function(master, membership, registry) {
  assert_schema(master,
    c("SampleID", "Site", "GeologicalHabitat", "SampleType", "Campaign", "Season",
      "FieldReplicate", "ReplicateStatus", "LocationID", "PairID_Substrate", "PairID_Season"),
    "master metadata")
  assert_unique(master$SampleID, "master SampleID")
  assert_true(nrow(master) == 75L, sprintf("Master metadata must have 75 rows, found %d", nrow(master)))
  assert_true(setequal(master$SampleType, c("Sea-sand", "Water")), "SampleType must be Sea-sand/Water")
  assert_true(setequal(master$Campaign, c("First", "Second")), "Campaign must be First/Second")
  assert_true(setequal(master$Season, c("Spring", "Autumn")), "Season must be Spring/Autumn")
  assert_true(all(master$Season[master$Campaign == "First"] == "Spring") &&
                all(master$Season[master$Campaign == "Second"] == "Autumn"),
              "Campaign-to-season mapping drift")
  location_counts <- table(master$LocationID)
  assert_true(length(location_counts) == 21L &&
                identical(sort(as.integer(location_counts)), sort(c(rep(4L, 15L), rep(3L, 3L), rep(2L, 3L)))),
              "LocationID must encode 21 fixed locations with the expected 4/3/2-observation pattern")
  pair_values <- master$PairID_Substrate[!is.na(master$PairID_Substrate) & nzchar(master$PairID_Substrate)]
  pair_counts <- table(pair_values)
  assert_true(length(pair_counts) == 36L && all(pair_counts == 2L),
              "PairID_Substrate must encode exactly 36 two-sample Sea-sand/Water pairs")
  unpaired <- master$SampleID[is.na(master$PairID_Substrate) | !nzchar(master$PairID_Substrate)]
  assert_true(setequal(unpaired, c("DD_s_S1", "DD_s_S2", "DD_s_S3")),
              "Only Dokdo Autumn Sea-sand samples may lack a substrate PairID")
  paired_rows <- !is.na(master$PairID_Substrate) & nzchar(master$PairID_Substrate)
  pair_types <- split(master$SampleType[paired_rows], master$PairID_Substrate[paired_rows])
  assert_true(all(vapply(pair_types, function(z) setequal(z, c("Sea-sand", "Water")), logical(1L))),
              "Each PairID_Substrate must contain one Sea-sand and one Water sample")
  season_values <- master$PairID_Season[!is.na(master$PairID_Season) & nzchar(master$PairID_Season)]
  season_counts <- table(season_values)
  assert_true(length(season_counts) == 33L && all(season_counts == 2L),
              "PairID_Season must encode exactly 33 two-season same-substrate pairs")
  season_rows <- !is.na(master$PairID_Season) & nzchar(master$PairID_Season)
  season_groups <- split(master[season_rows, , drop = FALSE], master$PairID_Season[season_rows])
  season_ok <- vapply(season_groups, function(x) {
    nrow(x) == 2L && length(unique(x$LocationID)) == 1L &&
      length(unique(x$SampleType)) == 1L && setequal(x$Season, c("Spring", "Autumn"))
  }, logical(1L))
  assert_true(all(season_ok), "Each PairID_Season must link one fixed-location Spring/Autumn pair of one substrate")
  season_unpaired <- master$SampleID[!season_rows]
  expected_season_unpaired <- c(paste0("ULD_s_F", 1:3), paste0("ULD_w_F", 1:3), paste0("DD_w_F", 1:3))
  assert_true(setequal(season_unpaired, expected_season_unpaired),
              "Across-season unmatched sample set drift")
  assert_true(!any(grepl("SS_BS", master$SampleID, fixed = TRUE)), "Excluded SS_BS sample found in master")
  master$SiteSeasonBlock <- paste(master$Site, master$Season, sep = "|")
  master$SiteCampaignBlock <- master$SiteSeasonBlock
  master$SiteSampleTypeBlock <- paste(master$Site, master$SampleType, sep = "|")

  assert_schema(membership,
    c("SetID", "UnitType", "SampleID", "Included", "Site", "SampleType", "Campaign"),
    "analysis-set membership")
  member_samples <- membership[membership$UnitType == "sample", , drop = FALSE]
  observed <- vapply(names(EXPECTED_SET_COUNTS), function(set_id) {
    rows <- member_samples[member_samples$SetID == set_id & member_samples$Included == "TRUE", , drop = FALSE]
    assert_unique(rows$SampleID, sprintf("membership %s SampleID", set_id))
    if (!all(rows$SampleID %in% master$SampleID)) stopf("Membership %s contains unknown SampleID", set_id)
    nrow(rows)
  }, integer(1L))
  if (!identical(unname(observed), unname(EXPECTED_SET_COUNTS))) {
    stopf("Analysis-set counts differ; observed=%s expected=%s",
          paste(sprintf("%s:%d", names(observed), observed), collapse = ","),
          paste(sprintf("%s:%d", names(EXPECTED_SET_COUNTS), EXPECTED_SET_COUNTS), collapse = ","))
  }

  assert_schema(registry,
    c("ModelID", "AnalysisFamily", "SetID", "Formula", "ClusterOrStrata", "PairIDsUsed"),
    "model registry")
  needed <- c("ALPHA_A72_BLOCK_CONTRAST", "ALPHA_A60_TYPE_BY_CAMPAIGN",
              "ALPHA_A36_SAND_CAMPAIGN", "ALPHA_A30_WATER_CAMPAIGN",
              "ALPHA_A42_FIRST_SAMPLETYPE", "BETA_SAMPLETYPE", "BETA_TYPE_CAMPAIGN",
              "BETA_SAND_CAMPAIGN", "BETA_WATER_CAMPAIGN", "BETA_SITE")
  if (!all(needed %in% registry$ModelID)) {
    stopf("Model registry is missing required ModelID(s): %s", paste(setdiff(needed, registry$ModelID), collapse = ", "))
  }
  paired_models <- c("ALPHA_A72_BLOCK_CONTRAST", "ALPHA_A60_TYPE_BY_CAMPAIGN",
                     "ALPHA_A36_SAND_CAMPAIGN", "ALPHA_A30_WATER_CAMPAIGN",
                     "ALPHA_A42_FIRST_SAMPLETYPE", "BETA_SAMPLETYPE",
                     "BETA_TYPE_CAMPAIGN", "BETA_SAND_CAMPAIGN", "BETA_WATER_CAMPAIGN")
  assert_true(all(registry$PairIDsUsed[match(paired_models, registry$ModelID)] == "TRUE"),
              "Confirmed substrate-paired models are not marked as using PairID_Substrate")
  beta_a72 <- registry[registry$ModelID == "BETA_SAMPLETYPE", , drop = FALSE]
  assert_true(nrow(beta_a72) == 1L &&
                identical(beta_a72$Formula[[1L]], "Distance ~ SiteSeasonBlock + SampleType") &&
                identical(beta_a72$ClusterOrStrata[[1L]], "PairID_Substrate"),
              "BETA_SAMPLETYPE registry is not aligned to the confirmed paired A72 PERMANOVA")
  beta_location <- registry[registry$ModelID %in%
    c("BETA_TYPE_CAMPAIGN", "BETA_SAND_CAMPAIGN", "BETA_WATER_CAMPAIGN"), , drop = FALSE]
  assert_true(nrow(beta_location) == 3L && all(beta_location$ClusterOrStrata == "LocationID"),
              "Season-sensitive beta models must restrict permutations within fixed LocationID")
  list(master = master, set_counts = observed)
}

load_inputs <- function(base) {
  paths <- input_paths(base)
  missing <- paths[!file.exists(paths)]
  if (length(missing)) stopf("Required input(s) missing:\n%s", paste(sprintf("  %s: %s", names(missing), missing), collapse = "\n"))
  master <- read_tsv(paths[["master_metadata"]])
  membership <- read_tsv(paths[["set_membership"]])
  registry <- read_tsv(paths[["model_registry"]])
  design <- validate_design_inputs(master, membership, registry)
  sample_ids <- design$master$SampleID
  counts <- list(
    `16S` = list(
      rarefied_primary = read_count_table(paths[["counts_16s_rarefied"]], sample_ids, "16S", "rarefied_primary"),
      unrarefied_tss_sensitivity = read_count_table(paths[["counts_16s_unrarefied"]], sample_ids, "16S", "unrarefied_tss_sensitivity")
    ),
    `18S` = list(
      rarefied_primary = read_count_table(paths[["counts_18s_rarefied"]], sample_ids, "18S", "rarefied_primary"),
      unrarefied_tss_sensitivity = read_count_table(paths[["counts_18s_unrarefied"]], sample_ids, "18S", "unrarefied_tss_sensitivity")
    )
  )
  assert_true(all(counts$`16S`$rarefied_primary$depths == 20000), "16S primary table is not exactly 20,000 reads/sample")
  assert_true(all(counts$`18S`$rarefied_primary$depths == 50000), "18S primary table is not exactly 50,000 reads/sample")
  for (marker in MARKER_ORDER) {
    a <- counts[[marker]]$rarefied_primary$matrix
    b <- counts[[marker]]$unrarefied_tss_sensitivity$matrix
    assert_true(identical(rownames(a), rownames(b)), sprintf("%s primary/sensitivity sample order differs", marker))
  }
  list(paths = paths, master = design$master, membership = membership,
       registry = registry, counts = counts, set_counts = design$set_counts)
}

sha256_file <- function(path) {
  assert_true(file.exists(path), sprintf("Cannot hash missing file: %s", path))
  tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

file_size <- function(path) as.numeric(file.info(path)$size)

write_tsv_new <- function(x, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA", fileEncoding = "UTF-8", eol = "\n")
}

write_text_new <- function(lines, path) {
  if (file.exists(path)) stopf("Refusing to overwrite file: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(lines), con = path, sep = "\n", useBytes = TRUE)
}

input_manifest <- function(inputs) {
  data.frame(
    Role = names(inputs$paths),
    AbsolutePath = normalizePath(unname(inputs$paths), winslash = "/", mustWork = TRUE),
    SizeBytes = vapply(unname(inputs$paths), file_size, numeric(1L)),
    SHA256 = vapply(unname(inputs$paths), sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
}

make_output_manifest <- function(root) {
  manifest_name <- if (grepl("alpha_diversity_v3$", root)) "output_sha256_manifest_v3.tsv" else "output_sha256_manifest_v4.tsv"
  manifest_path <- file.path(root, "provenance", manifest_name)
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[normalizePath(files, winslash = "/", mustWork = TRUE) !=
                   normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
  rel <- substring(normalizePath(files, winslash = "/", mustWork = TRUE),
                   nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L)
  ord <- order(rel)
  data.frame(RelativePath = rel[ord], SizeBytes = vapply(files[ord], file_size, numeric(1L)),
             SHA256 = vapply(files[ord], sha256_file, character(1L)), stringsAsFactors = FALSE)
}

analysis_config <- function(inputs, package_name) {
  data.frame(
    Key = c("script", "script_version", "package", "analysis_date", "r_version",
            "vegan_version", "permute_version", "digest_version", "markers",
            "alpha_count_surface", "beta_primary_surface", "beta_sensitivity_surface",
            "permutations", "base_seed", "sample_pairing", "season_labels", "campaign_provenance", "figure_backend"),
    Value = c(basename(script_path()), SCRIPT_VERSION, package_name,
              format(Sys.Date(), "%Y-%m-%d"), R.version.string,
              as.character(utils::packageVersion("vegan")),
              as.character(utils::packageVersion("permute")),
              as.character(utils::packageVersion("digest")), "16S|18S",
              "marker-specific fixed-depth rarefied counts",
              "Bray-Curtis on fixed-depth rarefied counts",
              "Bray-Curtis on sample-wise total-sum-scaled unrarefied counts",
              as.character(N_PERMUTATIONS), as.character(BASE_SEED),
              "confirmed fixed LocationID revisited in Spring/Autumn plus same-time Sea-sand/Water pairs; within-location differences precede site/block aggregation",
              "Spring|Autumn", "First|Second", "none (numerical module only)"),
    stringsAsFactors = FALSE
  )
}

schema_rows <- function(named_tables) {
  rows <- list(); k <- 0L
  for (table_name in names(named_tables)) {
    table <- named_tables[[table_name]]
    for (column_name in names(table)) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        Table = table_name,
        ColumnOrder = match(column_name, names(table)),
        Column = column_name,
        RStorageClass = paste(class(table[[column_name]]), collapse = "|"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

set_sample_ids <- function(inputs, set_id) {
  rows <- inputs$membership[
    inputs$membership$UnitType == "sample" & inputs$membership$SetID == set_id &
      inputs$membership$Included == "TRUE", , drop = FALSE]
  ids <- inputs$master$SampleID[inputs$master$SampleID %in% rows$SampleID]
  assert_true(length(ids) == EXPECTED_SET_COUNTS[[set_id]], sprintf("Runtime set count drift for %s", set_id))
  ids
}

alpha_one_sample <- function(x) {
  estimate <- suppressWarnings(vegan::estimateR(x))
  fetch <- function(name) {
    if (is.matrix(estimate)) {
      if (name %in% rownames(estimate)) return(as.numeric(estimate[name, 1L]))
      if (name %in% colnames(estimate)) return(as.numeric(estimate[1L, name]))
    }
    if (!is.null(names(estimate)) && name %in% names(estimate)) return(as.numeric(estimate[[name]]))
    NA_real_
  }
  c(
    Observed = sum(x > 0),
    Chao1 = fetch("S.chao1"),
    ACE = fetch("S.ACE"),
    Shannon = vegan::diversity(x, index = "shannon"),
    Simpson = vegan::diversity(x, index = "simpson"),
    InvSimpson = vegan::diversity(x, index = "invsimpson")
  )
}

compute_alpha_samples <- function(inputs) {
  pieces <- lapply(MARKER_ORDER, function(marker) {
    counts <- inputs$counts[[marker]]$rarefied_primary$matrix
    estimates <- t(vapply(seq_len(nrow(counts)), function(i) alpha_one_sample(counts[i, ]),
                          FUN.VALUE = setNames(numeric(length(METRIC_ORDER)), METRIC_ORDER)))
    assert_true(identical(colnames(estimates), METRIC_ORDER), sprintf("Alpha metric order drift for %s", marker))
    if (any(!is.finite(estimates))) {
      bad <- which(!is.finite(estimates), arr.ind = TRUE)
      stopf("Nonfinite alpha estimate for %s at %s/%s", marker,
            rownames(counts)[bad[1L, 1L]], colnames(estimates)[bad[1L, 2L]])
    }
    long <- data.frame(
      Marker = marker,
      SampleID = rep(rownames(counts), each = length(METRIC_ORDER)),
      Metric = rep(METRIC_ORDER, times = nrow(counts)),
      Value = as.vector(t(estimates)),
      CountSurface = "rarefied_primary",
      RarefactionDepth = if (marker == "16S") 20000L else 50000L,
      stringsAsFactors = FALSE
    )
    merge(long, inputs$master[, c("SampleID", "Site", "GeologicalHabitat", "SampleType",
                                  "Season", "Campaign", "FieldReplicate", "ReplicateStatus", "LocationID",
                                  "PairID_Substrate", "PairID_Season", "SiteSeasonBlock")],
          by = "SampleID", sort = FALSE, all.x = TRUE)
  })
  out <- do.call(rbind, pieces)
  out$Marker <- factor(out$Marker, levels = MARKER_ORDER)
  out$Metric <- factor(out$Metric, levels = METRIC_ORDER)
  out <- out[order(out$Marker, match(out$SampleID, inputs$master$SampleID), out$Metric), ]
  out$Marker <- as.character(out$Marker)
  out$Metric <- as.character(out$Metric)
  rownames(out) <- NULL
  out[, c("Marker", "SampleID", "Site", "GeologicalHabitat", "SampleType", "Season", "Campaign",
          "FieldReplicate", "ReplicateStatus", "LocationID", "PairID_Substrate", "PairID_Season", "SiteSeasonBlock",
          "Metric", "Value", "CountSurface", "RarefactionDepth")]
}

compute_cell_summaries <- function(alpha) {
  keys <- c("Marker", "Site", "GeologicalHabitat", "SampleType", "Season", "Campaign", "Metric")
  groups <- split(alpha, interaction(alpha[keys], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(x) {
    values <- x$Value
    data.frame(Marker = x$Marker[[1L]], Site = x$Site[[1L]],
               GeologicalHabitat = x$GeologicalHabitat[[1L]], SampleType = x$SampleType[[1L]],
               Season = x$Season[[1L]], Campaign = x$Campaign[[1L]], Metric = x$Metric[[1L]], N = length(values),
               Mean = mean(values), SD = stats::sd(values), Median = stats::median(values),
               IQR = stats::IQR(values), Min = min(values), Max = max(values),
               AggregationRule = "arithmetic summary of per-sample alpha values; never pooled counts",
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(match(out$Marker, MARKER_ORDER), match(out$Metric, METRIC_ORDER),
                   out$Site, match(out$Season, c("Spring", "Autumn")),
                   match(out$SampleType, c("Sea-sand", "Water"))), ]
  rownames(out) <- NULL
  out
}

cell_mean <- function(cells, marker, metric, site, campaign, sample_type) {
  hit <- cells$Marker == marker & cells$Metric == metric & cells$Site == site &
    cells$Campaign == campaign & cells$SampleType == sample_type
  if (sum(hit) != 1L) stopf("Expected one cell mean for %s/%s/%s/%s/%s; found %d",
                            marker, metric, site, campaign, sample_type, sum(hit))
  cells$Mean[hit]
}

alpha_contrast_specs <- function() {
  list(
    list(ModelID = "ALPHA_A72_BLOCK_CONTRAST", SetID = "A72_SUBSTRATE_PAIRED_CELLS",
         Contrast = "mean paired (Sea-sand - Water) within Site x Season",
         Unit = "SiteSeasonBlock", ExpectedN = 12L, PairIDsUsed = TRUE),
    list(ModelID = "ALPHA_A60_TYPE_BY_CAMPAIGN", SetID = "A60_COMPLETE_FACTORIAL",
         Contrast = "[mean paired (Sea-sand - Water)_Autumn - mean paired (Sea-sand - Water)_Spring]",
         Unit = "Site", ExpectedN = 5L, PairIDsUsed = TRUE),
    list(ModelID = "ALPHA_A36_SAND_CAMPAIGN", SetID = "A36_SAND_SEASON",
         Contrast = "mean fixed-location (Sea-sand Autumn - Spring)",
         Unit = "Site", ExpectedN = 6L, PairIDsUsed = TRUE),
    list(ModelID = "ALPHA_A30_WATER_CAMPAIGN", SetID = "A30_WATER_SEASON",
         Contrast = "mean fixed-location (Water Autumn - Spring)",
         Unit = "Site", ExpectedN = 5L, PairIDsUsed = TRUE),
    list(ModelID = "ALPHA_A42_FIRST_SAMPLETYPE", SetID = "A42_FIRST_SEASON_SITE",
         Contrast = "mean paired (Spring Sea-sand - Spring Water)",
         Unit = "Site", ExpectedN = 7L, PairIDsUsed = TRUE)
  )
}

compute_alpha_pair_contrasts <- function(alpha) {
  paired <- alpha[!is.na(alpha$PairID_Substrate) & nzchar(alpha$PairID_Substrate), , drop = FALSE]
  keys <- unique(paired[, c("Marker", "Metric", "PairID_Substrate")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- paired[paired$Marker == key$Marker & paired$Metric == key$Metric &
                  paired$PairID_Substrate == key$PairID_Substrate, , drop = FALSE]
    assert_true(nrow(x) == 2L && setequal(x$SampleType, c("Sea-sand", "Water")),
                sprintf("Invalid alpha pair %s/%s/%s", key$Marker, key$Metric, key$PairID_Substrate))
    sand <- x[x$SampleType == "Sea-sand", , drop = FALSE]
    water <- x[x$SampleType == "Water", , drop = FALSE]
    assert_true(sand$Site == water$Site && sand$Season == water$Season &&
                  sand$FieldReplicate == water$FieldReplicate,
                sprintf("Pair metadata mismatch: %s", key$PairID_Substrate))
    data.frame(Marker = key$Marker, Metric = key$Metric,
      PairID_Substrate = key$PairID_Substrate, LocationID = sand$LocationID,
      Site = sand$Site, Season = sand$Season,
      Campaign = sand$Campaign, FieldReplicate = sand$FieldReplicate,
      SeaSandSampleID = sand$SampleID, WaterSampleID = water$SampleID,
      SeaSandValue = sand$Value, WaterValue = water$Value,
      PairedDifferenceSeaSandMinusWater = sand$Value - water$Value,
      PairingEvidence = "confirmed same-number samples co-collected 1:1 at the same location and time",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  assert_true(nrow(out) == 2L * length(METRIC_ORDER) * 36L,
              sprintf("Expected 432 alpha pair rows, found %d", nrow(out)))
  assert_true(all(table(interaction(out$Marker, out$Metric, drop = TRUE)) == 36L),
              "Each marker x metric must have 36 substrate pairs")
  out <- out[order(match(out$Marker, MARKER_ORDER), match(out$Metric, METRIC_ORDER),
                   out$Site, match(out$Season, c("Spring", "Autumn")), out$FieldReplicate), ]
  rownames(out) <- NULL
  out
}

compute_alpha_season_contrasts <- function(alpha) {
  paired <- alpha[!is.na(alpha$PairID_Season) & nzchar(alpha$PairID_Season), , drop = FALSE]
  keys <- unique(paired[, c("Marker", "Metric", "PairID_Season")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    x <- paired[paired$Marker == key$Marker & paired$Metric == key$Metric &
                  paired$PairID_Season == key$PairID_Season, , drop = FALSE]
    assert_true(nrow(x) == 2L && setequal(x$Season, c("Spring", "Autumn")) &&
                  length(unique(x$SampleType)) == 1L && length(unique(x$LocationID)) == 1L,
                sprintf("Invalid alpha season pair %s/%s/%s", key$Marker, key$Metric, key$PairID_Season))
    spring <- x[x$Season == "Spring", , drop = FALSE]
    autumn <- x[x$Season == "Autumn", , drop = FALSE]
    data.frame(Marker = key$Marker, Metric = key$Metric,
      PairID_Season = key$PairID_Season, LocationID = spring$LocationID,
      Site = spring$Site, SampleType = spring$SampleType, FieldReplicate = spring$FieldReplicate,
      SpringSampleID = spring$SampleID, AutumnSampleID = autumn$SampleID,
      SpringValue = spring$Value, AutumnValue = autumn$Value,
      PairedDifferenceAutumnMinusSpring = autumn$Value - spring$Value,
      PairingEvidence = "confirmed fixed location revisited from Spring to Autumn",
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  assert_true(nrow(out) == 2L * length(METRIC_ORDER) * 33L,
              sprintf("Expected 396 alpha season-pair rows, found %d", nrow(out)))
  assert_true(all(table(interaction(out$Marker, out$Metric, drop = TRUE)) == 33L),
              "Each marker x metric must have 33 across-season pairs")
  out <- out[order(match(out$Marker, MARKER_ORDER), match(out$Metric, METRIC_ORDER),
                   out$Site, match(out$SampleType, c("Sea-sand", "Water")), out$FieldReplicate), ]
  rownames(out) <- NULL
  out
}

pair_block <- function(pairs, marker, metric, site, season) {
  x <- pairs[pairs$Marker == marker & pairs$Metric == metric & pairs$Site == site &
               pairs$Season == season, , drop = FALSE]
  if (nrow(x) != 3L) stopf("Expected three confirmed pairs for %s/%s/%s/%s; found %d",
                           marker, metric, site, season, nrow(x))
  list(delta = mean(x$PairedDifferenceSeaSandMinusWater), sand = mean(x$SeaSandValue),
       water = mean(x$WaterValue), n = nrow(x), ids = paste(x$PairID_Substrate, collapse = "|"))
}

season_block <- function(season_pairs, marker, metric, site, sample_type) {
  x <- season_pairs[season_pairs$Marker == marker & season_pairs$Metric == metric &
                      season_pairs$Site == site & season_pairs$SampleType == sample_type, , drop = FALSE]
  if (nrow(x) != 3L) stopf("Expected three fixed-location season pairs for %s/%s/%s/%s; found %d",
                           marker, metric, site, sample_type, nrow(x))
  list(delta = mean(x$PairedDifferenceAutumnMinusSpring), spring = mean(x$SpringValue),
       autumn = mean(x$AutumnValue), n = nrow(x), ids = paste(x$PairID_Season, collapse = "|"))
}

interaction_block <- function(substrate_pairs, marker, metric, site) {
  x <- substrate_pairs[substrate_pairs$Marker == marker & substrate_pairs$Metric == metric &
                         substrate_pairs$Site == site, , drop = FALSE]
  spring <- x[x$Season == "Spring", , drop = FALSE]
  autumn <- x[x$Season == "Autumn", , drop = FALSE]
  common <- intersect(spring$LocationID, autumn$LocationID)
  if (length(common) != 3L) stopf("Expected three repeated complete locations for %s/%s/%s; found %d",
                                  marker, metric, site, length(common))
  values <- vapply(common, function(location_id) {
    autumn$PairedDifferenceSeaSandMinusWater[autumn$LocationID == location_id] -
      spring$PairedDifferenceSeaSandMinusWater[spring$LocationID == location_id]
  }, numeric(1L))
  list(delta = mean(values), autumn = mean(autumn$PairedDifferenceSeaSandMinusWater[
    match(common, autumn$LocationID)]), spring = mean(spring$PairedDifferenceSeaSandMinusWater[
      match(common, spring$LocationID)]), n = length(common), ids = paste(common, collapse = "|"))
}

contrast_ledger_one <- function(cells, pairs, season_pairs, marker, metric, spec) {
  if (spec$ModelID == "ALPHA_A72_BLOCK_CONTRAST") {
    candidates <- unique(pairs[pairs$Marker == marker & pairs$Metric == metric,
                               c("Site", "Season", "Campaign"), drop = FALSE])
    rows <- lapply(seq_len(nrow(candidates)), function(i) {
      site <- candidates$Site[[i]]; season <- candidates$Season[[i]]; campaign <- candidates$Campaign[[i]]
      b <- pair_block(pairs, marker, metric, site, season)
      data.frame(IndependentUnitID = paste(site, season, sep = "|"), Site = site,
                 Season = season, Campaign = campaign, Component1 = b$sand, Component2 = b$water,
                 Component3 = NA_real_, Component4 = NA_real_, ContrastValue = b$delta,
                 ComponentDefinition = "C1=mean paired Sea-sand; C2=mean paired Water",
                 PairN = b$n, PairIDList = b$ids, stringsAsFactors = FALSE)
    })
  } else {
    eligible_sites <- switch(spec$ModelID,
      ALPHA_A60_TYPE_BY_CAMPAIGN = unique(cells$Site[
        cells$Marker == marker & cells$Metric == metric &
          ave(cells$SampleType, cells$Marker, cells$Metric, cells$Site,
              FUN = function(z) length(z)) >= 4L]),
      ALPHA_A36_SAND_CAMPAIGN = unique(cells$Site[cells$Marker == marker & cells$Metric == metric & cells$SampleType == "Sea-sand"]),
      ALPHA_A30_WATER_CAMPAIGN = unique(cells$Site[cells$Marker == marker & cells$Metric == metric & cells$SampleType == "Water"]),
      ALPHA_A42_FIRST_SAMPLETYPE = unique(cells$Site[cells$Marker == marker & cells$Metric == metric & cells$Campaign == "First"])
    )
    rows <- lapply(eligible_sites, function(site) {
      if (spec$ModelID == "ALPHA_A60_TYPE_BY_CAMPAIGN") {
        repeated <- interaction_block(pairs, marker, metric, site)
        data.frame(IndependentUnitID = site, Site = site, Season = "Spring|Autumn",
                   Campaign = "First|Second", Component1 = repeated$autumn, Component2 = repeated$spring,
                   Component3 = NA_real_, Component4 = NA_real_, ContrastValue = repeated$delta,
                   ComponentDefinition = "C1=mean within-location substrate difference_Autumn; C2=mean within-location substrate difference_Spring",
                   PairN = repeated$n, PairIDList = repeated$ids,
                   stringsAsFactors = FALSE)
      } else if (spec$ModelID == "ALPHA_A36_SAND_CAMPAIGN") {
        b <- tryCatch(season_block(season_pairs, marker, metric, site, "Sea-sand"), error = function(e) NULL)
        if (is.null(b)) return(NULL)
        data.frame(IndependentUnitID = site, Site = site, Season = "Spring|Autumn", Campaign = "First|Second",
                   Component1 = b$autumn, Component2 = b$spring, Component3 = NA_real_, Component4 = NA_real_,
                   ContrastValue = b$delta, ComponentDefinition = "C1=mean fixed-location Sea-sand_Autumn; C2=mean fixed-location Sea-sand_Spring",
                   PairN = b$n, PairIDList = b$ids,
                   stringsAsFactors = FALSE)
      } else if (spec$ModelID == "ALPHA_A30_WATER_CAMPAIGN") {
        b <- tryCatch(season_block(season_pairs, marker, metric, site, "Water"), error = function(e) NULL)
        if (is.null(b)) return(NULL)
        data.frame(IndependentUnitID = site, Site = site, Season = "Spring|Autumn", Campaign = "First|Second",
                   Component1 = b$autumn, Component2 = b$spring, Component3 = NA_real_, Component4 = NA_real_,
                   ContrastValue = b$delta, ComponentDefinition = "C1=mean fixed-location Water_Autumn; C2=mean fixed-location Water_Spring",
                   PairN = b$n, PairIDList = b$ids,
                   stringsAsFactors = FALSE)
      } else {
        spring <- pair_block(pairs, marker, metric, site, "Spring")
        data.frame(IndependentUnitID = site, Site = site, Season = "Spring", Campaign = "First",
                   Component1 = spring$sand, Component2 = spring$water,
                   Component3 = NA_real_, Component4 = NA_real_, ContrastValue = spring$delta,
                   ComponentDefinition = "C1=mean paired Spring Sea-sand; C2=mean paired Spring Water",
                   PairN = spring$n, PairIDList = spring$ids,
                   stringsAsFactors = FALSE)
      }
    })
  }
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out) != spec$ExpectedN) {
    stopf("%s %s %s contrast ledger expected n=%d, found %d",
          marker, metric, spec$ModelID, spec$ExpectedN, nrow(out))
  }
  cbind(data.frame(Marker = marker, Metric = metric, ModelID = spec$ModelID,
                   SetID = spec$SetID, ContrastDefinition = spec$Contrast,
                   IndependentUnit = spec$Unit, PairIDsUsed = spec$PairIDsUsed,
                   stringsAsFactors = FALSE), out)
}

compute_alpha_contrast_ledger <- function(cells, pairs, season_pairs) {
  rows <- list(); k <- 0L
  for (marker in MARKER_ORDER) for (metric in METRIC_ORDER) for (spec in alpha_contrast_specs()) {
    k <- k + 1L
    rows[[k]] <- contrast_ledger_one(cells, pairs, season_pairs, marker, metric, spec)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

exact_signflip_p <- function(values) {
  n <- length(values)
  assert_true(n > 0L && n <= 20L, sprintf("Unsupported exact sign-flip n=%d", n))
  pattern <- 0:(2^n - 1L)
  signed_means <- vapply(pattern, function(code) {
    signs <- ifelse(bitwAnd(code, bitwShiftL(1L, 0:(n - 1L))) == 0L, -1, 1)
    mean(signs * values)
  }, numeric(1L))
  observed <- abs(mean(values))
  mean(abs(signed_means) >= observed - TOLERANCE)
}

compute_alpha_tests <- function(ledger) {
  keys <- unique(ledger[, c("Marker", "Metric", "ModelID", "SetID", "ContrastDefinition",
                            "IndependentUnit", "PairIDsUsed")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    hit <- ledger$Marker == key$Marker & ledger$Metric == key$Metric & ledger$ModelID == key$ModelID
    values <- ledger$ContrastValue[hit]
    n <- length(values); effect <- mean(values); standard_error <- stats::sd(values) / sqrt(n)
    critical <- stats::qt(0.975, df = n - 1L)
    data.frame(Marker = key$Marker, Metric = key$Metric, ModelID = key$ModelID,
               SetID = key$SetID, ContrastDefinition = key$ContrastDefinition,
               IndependentUnit = key$IndependentUnit, NIndependentUnits = n,
               EffectMeanContrast = effect, SDContrasts = stats::sd(values),
               SEContrastMean = standard_error, CILevel = 0.95,
               CILower = effect - critical * standard_error,
               CIUpper = effect + critical * standard_error,
               CIMethod = "two-sided t interval for mean independent-unit contrast",
               PExactSignFlip = exact_signflip_p(values),
               Permutations = sprintf("%d_exact", 2^n), PairIDsUsed = key$PairIDsUsed,
               CountSurface = "rarefied_primary", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  split_key <- interaction(out$Marker, out$Metric, drop = TRUE)
  out$PAdjustedBH <- ave(out$PExactSignFlip, split_key, FUN = stats::p.adjust, method = "BH")
  out$PAdjustedBHWithinMarkerAllAlpha <- ave(out$PExactSignFlip, out$Marker,
                                             FUN = stats::p.adjust, method = "BH")
  out$MultiplicityFamily <- paste(out$Marker, out$Metric, "five_prespecified_contrasts", sep = "|")
  out <- out[order(match(out$Marker, MARKER_ORDER), match(out$Metric, METRIC_ORDER),
                   match(out$ModelID, vapply(alpha_contrast_specs(), `[[`, character(1L), "ModelID"))), ]
  rownames(out) <- NULL
  out
}

beta_model_specs <- function() {
  list(
    list(ModelID = "BETA_SAMPLETYPE", SetID = "A72_SUBSTRATE_PAIRED_CELLS",
         Formula = "DistanceObject ~ SiteSeasonBlock + SampleType", RHS = "~ SiteSeasonBlock + SampleType",
         Strata = "PairID_Substrate", DispersionFactor = "SampleType", PairIDsUsed = TRUE),
    list(ModelID = "BETA_TYPE_CAMPAIGN", SetID = "A60_COMPLETE_FACTORIAL",
         Formula = "DistanceObject ~ Site + SampleType * Season", RHS = "~ Site + SampleType * Season",
         Strata = "LocationID", DispersionFactor = "SampleType:Season", PairIDsUsed = TRUE),
    list(ModelID = "BETA_SAND_CAMPAIGN", SetID = "A36_SAND_SEASON",
         Formula = "DistanceObject ~ Site + Season", RHS = "~ Site + Season",
         Strata = "LocationID", DispersionFactor = "Season", PairIDsUsed = TRUE),
    list(ModelID = "BETA_WATER_CAMPAIGN", SetID = "A30_WATER_SEASON",
         Formula = "DistanceObject ~ Site + Season", RHS = "~ Site + Season",
         Strata = "LocationID", DispersionFactor = "Season", PairIDsUsed = TRUE),
    list(ModelID = "BETA_SITE", SetID = "A42_FIRST_SEASON_SITE",
         Formula = "DistanceObject ~ SampleType + Site", RHS = "~ SampleType + Site",
         Strata = "SampleType", DispersionFactor = "Site", PairIDsUsed = FALSE)
  )
}

surface_matrix <- function(inputs, marker, surface) {
  counts <- inputs$counts[[marker]][[surface]]$matrix
  if (surface == "rarefied_primary") {
    transformed <- counts
    normalization <- "fixed-depth rarefied counts"
  } else {
    transformed <- sweep(counts, 1L, rowSums(counts), "/")
    normalization <- "sample-wise total-sum scaling of unrarefied counts"
  }
  list(matrix = transformed, normalization = normalization)
}

distance_for_set <- function(inputs, marker, surface, set_id) {
  ids <- set_sample_ids(inputs, set_id)
  surf <- surface_matrix(inputs, marker, surface)
  matrix <- surf$matrix[ids, , drop = FALSE]
  assert_true(identical(rownames(matrix), ids), "Distance matrix sample order drift")
  list(distance = vegan::vegdist(matrix, method = "bray"), ids = ids,
       normalization = surf$normalization)
}

compute_pcoa <- function(inputs) {
  coordinates <- list(); eigenvalues <- list(); registry <- list(); k <- 0L
  for (marker in MARKER_ORDER) for (surface in SURFACE_ORDER) for (set_id in SET_ORDER) {
    k <- k + 1L
    item <- distance_for_set(inputs, marker, surface, set_id)
    raw <- vegan::wcmdscale(item$distance, eig = TRUE, add = FALSE)
    lingoes <- vegan::wcmdscale(item$distance, eig = TRUE, add = "lingoes")
    points <- as.matrix(lingoes$points)
    if (is.null(colnames(points))) colnames(points) <- paste0("Dim", seq_len(ncol(points)))
    coord <- data.frame(Marker = marker, Surface = surface, SetID = set_id,
                        SampleID = rep(rownames(points), times = ncol(points)),
                        Axis = rep(seq_len(ncol(points)), each = nrow(points)),
                        Coordinate = as.vector(points), Correction = "Lingoes",
                        AdditiveConstant = if (is.null(lingoes$ac)) 0 else as.numeric(lingoes$ac),
                        Normalization = item$normalization, stringsAsFactors = FALSE)
    coordinates[[k]] <- merge(coord,
      inputs$master[, c("SampleID", "Site", "GeologicalHabitat", "SampleType", "Season", "Campaign",
                        "FieldReplicate", "LocationID", "PairID_Substrate", "PairID_Season", "SiteSeasonBlock")],
      by = "SampleID", sort = FALSE, all.x = TRUE)
    eigen_one <- function(fit, correction) {
      eig <- as.numeric(fit$eig)
      positive_total <- sum(eig[eig > 0])
      data.frame(Marker = marker, Surface = surface, SetID = set_id,
                 Correction = correction, Axis = seq_along(eig), Eigenvalue = eig,
                 RelativePositiveVariance = ifelse(eig > 0 & positive_total > 0, eig / positive_total, NA_real_),
                 IsNegative = eig < -TOLERANCE,
                 AdditiveConstant = if (is.null(fit$ac)) 0 else as.numeric(fit$ac),
                 Normalization = item$normalization, stringsAsFactors = FALSE)
    }
    eigenvalues[[2L * k - 1L]] <- eigen_one(raw, "none")
    eigenvalues[[2L * k]] <- eigen_one(lingoes, "Lingoes")
    registry[[k]] <- data.frame(Marker = marker, Surface = surface, SetID = set_id,
      SampleN = length(item$ids), FeatureN = ncol(inputs$counts[[marker]][[surface]]$matrix),
      Distance = "Bray-Curtis", Normalization = item$normalization, Correction = "Lingoes",
      AdditiveConstant = if (is.null(lingoes$ac)) 0 else as.numeric(lingoes$ac),
      CoordinateAxes = ncol(points), stringsAsFactors = FALSE)
  }
  list(coordinates = do.call(rbind, coordinates), eigenvalues = do.call(rbind, eigenvalues),
       surface_registry = do.call(rbind, registry))
}

permutation_control <- function(metadata, strata, seed) {
  assert_true(strata %in% names(metadata), sprintf("Unknown strata field: %s", strata))
  set.seed(seed)
  permute::how(nperm = N_PERMUTATIONS, blocks = factor(metadata[[strata]]))
}

seed_for <- function(marker, surface, model_index, offset = 0L) {
  BASE_SEED + (match(marker, MARKER_ORDER) - 1L) * 1000L +
    (match(surface, SURFACE_ORDER) - 1L) * 100L + model_index * 10L + offset
}

adonis_rows <- function(fit, marker, surface, spec, seed, normalization) {
  tab <- as.data.frame(fit)
  tab$Term <- rownames(tab)
  rownames(tab) <- NULL
  get_col <- function(candidates) {
    hit <- intersect(candidates, names(tab))
    if (!length(hit)) return(rep(NA_real_, nrow(tab)))
    as.numeric(tab[[hit[[1L]]]])
  }
  out <- data.frame(Marker = marker, Surface = surface, ModelID = spec$ModelID,
    SetID = spec$SetID, Formula = spec$Formula, Strata = spec$Strata,
    PermutationRestriction = paste("within", spec$Strata), PermutationsRequested = N_PERMUTATIONS,
    Seed = seed, By = "margin", Term = tab$Term, Df = get_col(c("Df")),
    SumOfSquares = get_col(c("SumOfSqs", "SumsOfSqs")), R2 = get_col(c("R2")),
    F = get_col(c("F")), P = get_col(c("Pr(>F)")), Distance = "Bray-Curtis",
    Normalization = normalization, PairIDsUsed = spec$PairIDsUsed, stringsAsFactors = FALSE)
  inferential <- !(out$Term %in% c("Residual", "Total")) & is.finite(out$P)
  out$PAdjustedBHWithinModel <- NA_real_
  if (any(inferential)) out$PAdjustedBHWithinModel[inferential] <- stats::p.adjust(out$P[inferential], method = "BH")
  out
}

permdisp_rows <- function(distance, metadata, marker, surface, spec, seed, normalization) {
  group <- if (spec$DispersionFactor == "SampleType:Season") {
    interaction(metadata$SampleType, metadata$Season, sep = "|", drop = TRUE)
  } else {
    factor(metadata[[spec$DispersionFactor]])
  }
  dispersion <- vegan::betadisper(distance, group = group, type = "median", bias.adjust = TRUE)
  control <- permutation_control(metadata, spec$Strata, seed)
  set.seed(seed)
  test <- vegan::permutest(dispersion, permutations = control)
  tab <- as.data.frame(test$tab)
  result <- data.frame(Marker = marker, Surface = surface, ModelID = spec$ModelID,
    SetID = spec$SetID, Formula = spec$Formula, TargetFactor = spec$DispersionFactor, Strata = spec$Strata,
    PermutationRestriction = paste("within", spec$Strata), PermutationsRequested = N_PERMUTATIONS,
    Seed = seed, Method = "betadisper(type=median,bias.adjust=TRUE) + permutest",
    Df = as.numeric(tab[1L, intersect(c("Df"), names(tab))[[1L]]]),
    SumOfSquares = as.numeric(tab[1L, intersect(c("Sum Sq", "SumOfSqs"), names(tab))[[1L]]]),
    MeanSquare = as.numeric(tab[1L, intersect(c("Mean Sq", "MeanSqs"), names(tab))[[1L]]]),
    F = as.numeric(tab[1L, intersect(c("F", "F value"), names(tab))[[1L]]]),
    P = as.numeric(tab[1L, intersect(c("Pr(>F)"), names(tab))[[1L]]]),
    Distance = "Bray-Curtis", Normalization = normalization, PairIDsUsed = spec$PairIDsUsed,
    stringsAsFactors = FALSE)
  sample_rows <- data.frame(Marker = marker, Surface = surface, ModelID = spec$ModelID,
    SetID = spec$SetID, Formula = spec$Formula, SampleID = names(dispersion$distances),
    TargetFactor = spec$DispersionFactor, TargetLevel = as.character(group),
    DistanceToGroupMedian = as.numeric(dispersion$distances), Strata = spec$Strata,
    StrataLevel = as.character(metadata[[spec$Strata]]), Distance = "Bray-Curtis",
    Normalization = normalization, PairIDsUsed = spec$PairIDsUsed, stringsAsFactors = FALSE)
  list(result = result, samples = sample_rows)
}

compute_beta_tests <- function(inputs) {
  permanova <- list(); dispersion <- list(); dispersion_samples <- list(); k <- 0L
  specs <- beta_model_specs()
  for (marker in MARKER_ORDER) for (surface in SURFACE_ORDER) for (i in seq_along(specs)) {
    spec <- specs[[i]]; k <- k + 1L
    item <- distance_for_set(inputs, marker, surface, spec$SetID)
    metadata <- inputs$master[match(item$ids, inputs$master$SampleID), , drop = FALSE]
    assert_true(identical(metadata$SampleID, item$ids), "PERMANOVA metadata order drift")
    DistanceObject <- item$distance
    formula <- stats::as.formula(paste("DistanceObject", spec$RHS))
    environment(formula) <- environment()
    adonis_seed <- seed_for(marker, surface, i, 1L)
    control <- permutation_control(metadata, spec$Strata, adonis_seed)
    set.seed(adonis_seed)
    fit <- vegan::adonis2(formula, data = metadata, permutations = control,
                          method = "bray", by = "margin")
    permanova[[k]] <- adonis_rows(fit, marker, surface, spec, adonis_seed, item$normalization)
    dispersion_seed <- seed_for(marker, surface, i, 2L)
    disp <- permdisp_rows(item$distance, metadata, marker, surface, spec,
                          dispersion_seed, item$normalization)
    dispersion[[k]] <- disp$result
    dispersion_samples[[k]] <- disp$samples
  }
  perma <- do.call(rbind, permanova)
  target_terms <- c(BETA_SAMPLETYPE = "SampleType", BETA_TYPE_CAMPAIGN = "SampleType:Season",
                    BETA_SAND_CAMPAIGN = "Season", BETA_WATER_CAMPAIGN = "Season",
                    BETA_SITE = "Site")
  perma$PAlgorithmic <- perma$P
  perma$PAdjustedBHWithinModelAlgorithmic <- perma$PAdjustedBHWithinModel
  perma$TermRole <- ifelse(perma$Term == "Residual", "residual",
                    ifelse(perma$Term == "Total", "total",
                      ifelse(perma$Term == unname(target_terms[perma$ModelID]), "target",
                             "nuisance_fixed_adjustment")))
  perma$InferenceReportable <- perma$TermRole == "target"
  perma$PReportable <- ifelse(perma$InferenceReportable, perma$PAlgorithmic, NA_real_)
  perma$PAdjustedBHTargetWithinMarkerSurface <- NA_real_
  perma$PAdjustedBHWithinMarkerSurface <- NA_real_
  eligible <- is.finite(perma$P) & !(perma$Term %in% c("Residual", "Total"))
  groups <- interaction(perma$Marker, perma$Surface, drop = TRUE)
  for (g in unique(groups[eligible])) {
    hit <- eligible & groups == g
    perma$PAdjustedBHWithinMarkerSurface[hit] <- stats::p.adjust(perma$P[hit], method = "BH")
  }
  reportable <- perma$InferenceReportable & is.finite(perma$PReportable)
  for (g in unique(groups[reportable])) {
    hit <- reportable & groups == g
    assert_true(sum(hit) == 5L, sprintf("Expected five reportable beta targets in %s", g))
    perma$PAdjustedBHTargetWithinMarkerSurface[hit] <- stats::p.adjust(perma$PReportable[hit], method = "BH")
  }
  disp <- do.call(rbind, dispersion)
  groups <- interaction(disp$Marker, disp$Surface, drop = TRUE)
  disp$TermRole <- "target_dispersion_audit"
  disp$InferenceReportable <- TRUE
  disp$PAlgorithmic <- disp$P
  disp$PReportable <- disp$P
  disp$PAdjustedBHTargetWithinMarkerSurface <- ave(disp$PReportable, groups, FUN = stats::p.adjust, method = "BH")
  disp$PAdjustedBHWithinMarkerSurface <- disp$PAdjustedBHTargetWithinMarkerSurface
  list(permanova = perma, permdisp = disp, permdisp_samples = do.call(rbind, dispersion_samples))
}

alpha_readme <- function() c(
  "# Alpha-diversity package v3", "",
  "This versioned package contains figure-free, per-sample alpha-diversity estimates for the protected shared 75-sample 16S/18S design.", "",
  "- 16S uses the 20,000-read rarefied table; 18S uses the 50,000-read rarefied table.",
  "- Observed, Chao1, ACE, Shannon, Simpson, and inverse Simpson are computed for each sample before any grouping.",
  "- Cell summaries are arithmetic summaries of replicate-level metrics; counts are never pooled to create a cell metric.",
  "- The three suffixes are confirmed biological spatial field replicates collected by bulk harvest about 0.5-1 m apart.",
  "- Same-number Sea-sand and Water samples are confirmed 1:1 co-collections; substrate pair differences are computed first.",
  "- The same fixed LocationID was revisited in Spring and Autumn; seasonal and interaction differences are computed within LocationID before site aggregation.",
  "- Spring and Autumn are reader-facing seasons; First and Second are retained as source campaign provenance.",
  "- Incomplete cells remain unmatched; no missing visit or substrate is imputed.",
  "- P values are exhaustive two-sided sign-flip tests. The 95% interval is a t interval for the mean independent-unit contrast and is labeled separately from the randomization p value.",
  "- BH adjustment is reported both across the five prespecified contrasts within each Marker x Metric family (primary) and across all alpha tests within each marker (global sensitivity).", "",
  "Run `Rscript scripts/run_alpha_beta_hierarchical_v1.R --check` before using these tables in figures or manuscript text."
)

beta_readme <- function() c(
  "# Beta-diversity package v4", "",
  "This versioned package contains figure-free Bray-Curtis ordination and restricted permutation analyses for the protected shared 75-sample 16S/18S design.", "",
  "- Primary surface: Bray-Curtis on marker-specific fixed-depth rarefied counts.",
  "- Sensitivity surface: Bray-Curtis on sample-wise total-sum-scaled unrarefied counts.",
  "- PCoA coordinates use a Lingoes correction; both uncorrected and corrected eigenvalue ledgers are retained.",
  "- The A72 SampleType PERMANOVA fits Site x Season before SampleType and permutes labels within the 36 confirmed substrate PairIDs.",
  "- A60 interaction and A36/A30 season models restrict permutations within confirmed fixed LocationID; incomplete visits are excluded by their registered sets.",
  "- PERMDISP audits the target factor for each model using distances to spatial medians and the same restriction field.",
  "- Only prespecified target terms have reportable p-values; nuisance fixed-adjustment p-values remain algorithmic audit fields.",
  "- Habitat-level A7 results remain descriptive and are not tested here.", "",
  "Run `Rscript scripts/run_alpha_beta_hierarchical_v1.R --check` before plotting or interpretation."
)

build_qc <- function(inputs, alpha = NULL, beta = NULL, package_name) {
  rows <- data.frame(Check = character(), Observed = character(), Expected = character(),
                     Status = character(), stringsAsFactors = FALSE)
  add <- function(check, observed, expected, status) {
    rows <<- rbind(rows, data.frame(Check = check, Observed = as.character(observed),
      Expected = as.character(expected), Status = status, stringsAsFactors = FALSE))
  }
  add("master_sample_count", nrow(inputs$master), 75, ifelse(nrow(inputs$master) == 75, "PASS", "FAIL"))
  for (set_id in names(EXPECTED_SET_COUNTS)) {
    observed <- inputs$set_counts[[set_id]]
    add(paste0("set_count_", set_id), observed, EXPECTED_SET_COUNTS[[set_id]],
        ifelse(observed == EXPECTED_SET_COUNTS[[set_id]], "PASS", "FAIL"))
  }
  observed_pairs <- length(unique(inputs$master$PairID_Substrate[
    !is.na(inputs$master$PairID_Substrate) & nzchar(inputs$master$PairID_Substrate)]))
  add("confirmed_substrate_pair_ids", observed_pairs, 36L,
      ifelse(observed_pairs == 36L, "PASS", "FAIL"))
  observed_season_pairs <- length(unique(inputs$master$PairID_Season[
    !is.na(inputs$master$PairID_Season) & nzchar(inputs$master$PairID_Season)]))
  add("confirmed_season_pair_ids", observed_season_pairs, 33L,
      ifelse(observed_season_pairs == 33L, "PASS", "FAIL"))
  add("fixed_location_ids", length(unique(inputs$master$LocationID)), 21L,
      ifelse(length(unique(inputs$master$LocationID)) == 21L, "PASS", "FAIL"))
  if (package_name == "alpha") {
    add("sample_metric_rows", nrow(alpha$samples), 2L * 75L * 6L,
        ifelse(nrow(alpha$samples) == 900L, "PASS", "FAIL"))
    add("cell_summary_rows", nrow(alpha$cells), 2L * 25L * 6L,
        ifelse(nrow(alpha$cells) == 300L, "PASS", "FAIL"))
    add("confirmed_pair_contrast_rows", nrow(alpha$pairs), 2L * 6L * 36L,
        ifelse(nrow(alpha$pairs) == 432L, "PASS", "FAIL"))
    add("confirmed_season_pair_contrast_rows", nrow(alpha$season_pairs), 2L * 6L * 33L,
        ifelse(nrow(alpha$season_pairs) == 396L, "PASS", "FAIL"))
    add("contrast_ledger_rows", nrow(alpha$ledger), 2L * 6L * sum(c(12L, 5L, 6L, 5L, 7L)),
        ifelse(nrow(alpha$ledger) == 420L, "PASS", "FAIL"))
    add("contrast_result_rows", nrow(alpha$tests), 2L * 6L * 5L,
        ifelse(nrow(alpha$tests) == 60L, "PASS", "FAIL"))
  } else {
    expected_models <- 2L * 2L * 5L
    add("pcoa_surfaces", nrow(beta$pcoa$surface_registry), 2L * 2L * 6L,
        ifelse(nrow(beta$pcoa$surface_registry) == 24L, "PASS", "FAIL"))
    add("permdisp_model_rows", nrow(beta$tests$permdisp), expected_models,
        ifelse(nrow(beta$tests$permdisp) == expected_models, "PASS", "FAIL"))
    add("permanova_models_present", length(unique(paste(beta$tests$permanova$Marker,
      beta$tests$permanova$Surface, beta$tests$permanova$ModelID))), expected_models,
      ifelse(length(unique(paste(beta$tests$permanova$Marker, beta$tests$permanova$Surface,
                                 beta$tests$permanova$ModelID))) == expected_models, "PASS", "FAIL"))
  }
  rows
}

write_common_provenance <- function(stage, inputs, package_name, qc) {
  suffix <- if (package_name == "alpha") "v3" else "v4"
  write_tsv_new(input_manifest(inputs), file.path(stage, "provenance", paste0("input_sha256_manifest_", suffix, ".tsv")))
  write_tsv_new(analysis_config(inputs, package_name), file.path(stage, "provenance", paste0("analysis_config_", suffix, ".tsv")))
  write_text_new(capture.output(utils::sessionInfo()), file.path(stage, "provenance", paste0("r_session_info_", suffix, ".txt")))
  write_tsv_new(qc, file.path(stage, "qc", paste0(package_name, "_qc_summary_", suffix, ".tsv")))
}

write_alpha_package <- function(stage, inputs, alpha) {
  write_tsv_new(alpha$samples, file.path(stage, "results", "alpha_sample_metrics_v3.tsv"))
  write_tsv_new(alpha$cells, file.path(stage, "results", "alpha_cell_summaries_v3.tsv"))
  write_tsv_new(alpha$pairs, file.path(stage, "results", "alpha_substrate_pair_contrasts_v3.tsv"))
  write_tsv_new(alpha$season_pairs, file.path(stage, "results", "alpha_season_pair_contrasts_v3.tsv"))
  write_tsv_new(alpha$ledger, file.path(stage, "results", "alpha_independent_contrast_ledger_v3.tsv"))
  write_tsv_new(alpha$tests, file.path(stage, "results", "alpha_exact_signflip_results_v3.tsv"))
  write_text_new(alpha_readme(), file.path(stage, "README_v3.md"))
  qc <- build_qc(inputs, alpha = alpha, package_name = "alpha")
  assert_true(all(qc$Status == "PASS"), "Alpha QC did not pass before package finalization")
  write_common_provenance(stage, inputs, "alpha", qc)
  write_tsv_new(schema_rows(list(
    alpha_sample_metrics_v3 = alpha$samples,
    alpha_cell_summaries_v3 = alpha$cells,
    alpha_substrate_pair_contrasts_v3 = alpha$pairs,
    alpha_season_pair_contrasts_v3 = alpha$season_pairs,
    alpha_independent_contrast_ledger_v3 = alpha$ledger,
    alpha_exact_signflip_results_v3 = alpha$tests,
    alpha_qc_summary_v3 = qc
  )), file.path(stage, "provenance", "output_schema_v3.tsv"))
  write_tsv_new(make_output_manifest(stage), file.path(stage, "provenance", "output_sha256_manifest_v3.tsv"))
}

write_beta_package <- function(stage, inputs, beta) {
  write_tsv_new(beta$pcoa$coordinates, file.path(stage, "results", "beta_pcoa_coordinates_lingoes_v4.tsv"))
  write_tsv_new(beta$pcoa$eigenvalues, file.path(stage, "results", "beta_pcoa_eigenvalue_ledger_v4.tsv"))
  write_tsv_new(beta$pcoa$surface_registry, file.path(stage, "results", "beta_surface_registry_v4.tsv"))
  write_tsv_new(beta$tests$permanova, file.path(stage, "results", "beta_permanova_restricted_v4.tsv"))
  write_tsv_new(beta$tests$permdisp, file.path(stage, "results", "beta_permdisp_restricted_v4.tsv"))
  write_tsv_new(beta$tests$permdisp_samples, file.path(stage, "results", "beta_permdisp_sample_distances_v4.tsv"))
  inference_contract <- data.frame(
    RuleID = c("BETA-V4-001", "BETA-V4-002", "BETA-V4-003", "BETA-V4-004", "BETA-V4-005"),
    Rule = c("Only rows with InferenceReportable=TRUE are manuscript-facing inferential tests.",
             "A72 SampleType permutations are restricted within confirmed PairID_Substrate.",
             "A60 interaction and A36/A30 season permutations are restricted within confirmed fixed LocationID.",
             "Nuisance fixed-adjustment terms retain algorithmic statistics but have no reportable p-value.",
             "PERMANOVA targets are interpreted with matched PERMDISP and target-family BH correction."),
    stringsAsFactors = FALSE)
  write_tsv_new(inference_contract, file.path(stage, "results", "beta_inference_contract_v4.tsv"))
  write_text_new(beta_readme(), file.path(stage, "README_v4.md"))
  qc <- build_qc(inputs, beta = beta, package_name = "beta")
  assert_true(all(qc$Status == "PASS"), "Beta QC did not pass before package finalization")
  write_common_provenance(stage, inputs, "beta", qc)
  write_tsv_new(schema_rows(list(
    beta_pcoa_coordinates_lingoes_v4 = beta$pcoa$coordinates,
    beta_pcoa_eigenvalue_ledger_v4 = beta$pcoa$eigenvalues,
    beta_surface_registry_v4 = beta$pcoa$surface_registry,
    beta_permanova_restricted_v4 = beta$tests$permanova,
    beta_permdisp_restricted_v4 = beta$tests$permdisp,
    beta_permdisp_sample_distances_v4 = beta$tests$permdisp_samples,
    beta_inference_contract_v4 = inference_contract,
    beta_qc_summary_v4 = qc
  )), file.path(stage, "provenance", "output_schema_v4.tsv"))
  write_tsv_new(make_output_manifest(stage), file.path(stage, "provenance", "output_sha256_manifest_v4.tsv"))
}

safe_remove_stage <- function(path, base, target_basename) {
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  base_norm <- normalizePath(base, winslash = "/", mustWork = TRUE)
  allowed_prefix <- paste0(base_norm, "/.", target_basename, ".staging_")
  if (dir.exists(path) && startsWith(path_norm, allowed_prefix)) unlink(path, recursive = TRUE, force = TRUE)
}

build_all <- function(inputs, base) {
  finals <- output_paths(base)
  existing <- finals[dir.exists(finals) | file.exists(finals)]
  if (length(existing)) stopf("Refusing to overwrite existing output(s):\n%s", paste(existing, collapse = "\n"))
  tag <- paste0(Sys.getpid(), "_", format(Sys.time(), "%Y%m%dT%H%M%S"))
  stages <- c(alpha = file.path(base, paste0(".", basename(finals[["alpha"]]), ".staging_", tag)),
              beta = file.path(base, paste0(".", basename(finals[["beta"]]), ".staging_", tag)))
  if (any(dir.exists(stages) | file.exists(stages))) stopf("Staging path collision")
  dir.create(stages[["alpha"]], recursive = FALSE, showWarnings = FALSE)
  dir.create(stages[["beta"]], recursive = FALSE, showWarnings = FALSE)
  committed_alpha <- FALSE
  success <- FALSE
  on.exit({
    if (!success) {
      if (committed_alpha && dir.exists(finals[["alpha"]]) && !dir.exists(stages[["alpha"]])) {
        file.rename(finals[["alpha"]], stages[["alpha"]])
      }
      safe_remove_stage(stages[["alpha"]], base, basename(finals[["alpha"]]))
      safe_remove_stage(stages[["beta"]], base, basename(finals[["beta"]]))
    }
  }, add = TRUE)

  alpha_samples <- compute_alpha_samples(inputs)
  alpha_cells <- compute_cell_summaries(alpha_samples)
  alpha_pairs <- compute_alpha_pair_contrasts(alpha_samples)
  alpha_season_pairs <- compute_alpha_season_contrasts(alpha_samples)
  alpha_ledger <- compute_alpha_contrast_ledger(alpha_cells, alpha_pairs, alpha_season_pairs)
  alpha_tests <- compute_alpha_tests(alpha_ledger)
  alpha <- list(samples = alpha_samples, cells = alpha_cells, pairs = alpha_pairs,
                season_pairs = alpha_season_pairs, ledger = alpha_ledger, tests = alpha_tests)
  beta <- list(pcoa = compute_pcoa(inputs), tests = compute_beta_tests(inputs))
  write_alpha_package(stages[["alpha"]], inputs, alpha)
  write_beta_package(stages[["beta"]], inputs, beta)

  if (!file.rename(stages[["alpha"]], finals[["alpha"]])) stopf("Could not commit alpha staging directory")
  committed_alpha <- TRUE
  if (!file.rename(stages[["beta"]], finals[["beta"]])) stopf("Could not commit beta staging directory")
  success <- TRUE
  cat(sprintf("Created %s\nCreated %s\n", finals[["alpha"]], finals[["beta"]]))
}

verify_output_manifest <- function(root) {
  suffix <- if (grepl("alpha_diversity_v3", root, fixed = TRUE)) "v3" else "v4"
  path <- file.path(root, "provenance", paste0("output_sha256_manifest_", suffix, ".tsv"))
  manifest <- read_tsv(path)
  assert_schema(manifest, c("RelativePath", "SizeBytes", "SHA256"), paste(basename(root), "output manifest"))
  assert_unique(manifest$RelativePath, paste(basename(root), "manifest RelativePath"))
  for (i in seq_len(nrow(manifest))) {
    target <- file.path(root, manifest$RelativePath[[i]])
    assert_true(file.exists(target), sprintf("Manifest target missing: %s", target))
    assert_true(file_size(target) == as.numeric(manifest$SizeBytes[[i]]), sprintf("Size mismatch: %s", target))
    assert_true(identical(sha256_file(target), tolower(manifest$SHA256[[i]])), sprintf("SHA256 mismatch: %s", target))
  }
  actual <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  actual <- actual[file.info(actual)$isdir %in% FALSE]
  actual_rel <- substring(normalizePath(actual, winslash = "/", mustWork = TRUE),
                          nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L)
  expected_rel <- c(manifest$RelativePath, paste0("provenance/output_sha256_manifest_", suffix, ".tsv"))
  if (!setequal(actual_rel, expected_rel)) {
    stopf("Output manifest coverage differs in %s; unmanifested=%s missing=%s", root,
          paste(setdiff(actual_rel, expected_rel), collapse = "|"),
          paste(setdiff(expected_rel, actual_rel), collapse = "|"))
  }
  invisible(TRUE)
}

verify_input_manifest <- function(root) {
  suffix <- if (grepl("alpha_diversity_v3", root, fixed = TRUE)) "v3" else "v4"
  manifest <- read_tsv(file.path(root, "provenance", paste0("input_sha256_manifest_", suffix, ".tsv")))
  assert_schema(manifest, c("Role", "AbsolutePath", "SizeBytes", "SHA256"), "input manifest")
  for (i in seq_len(nrow(manifest))) {
    path <- manifest$AbsolutePath[[i]]
    assert_true(file.exists(path), sprintf("Current input missing: %s", path))
    assert_true(file_size(path) == as.numeric(manifest$SizeBytes[[i]]), sprintf("Current input size drift: %s", path))
    assert_true(identical(sha256_file(path), tolower(manifest$SHA256[[i]])), sprintf("Current input SHA256 drift: %s", path))
  }
}

check_existing <- function(base) {
  finals <- output_paths(base)
  if (!all(dir.exists(finals))) stopf("Both final output directories are required for --check")
  for (root in finals) {
    suffix <- if (root == finals[["alpha"]]) "v3" else "v4"
    verify_output_manifest(root)
    verify_input_manifest(root)
    qc_path <- file.path(root, "qc", paste0(ifelse(root == finals[["alpha"]], "alpha", "beta"), "_qc_summary_", suffix, ".tsv"))
    qc <- read_tsv(qc_path)
    assert_schema(qc, c("Check", "Observed", "Expected", "Status"), "QC summary")
    assert_true(all(qc$Status == "PASS"), sprintf("QC summary contains non-PASS row: %s", qc_path))
  }
  alpha_samples <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_sample_metrics_v3.tsv"))
  alpha_cells <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_cell_summaries_v3.tsv"))
  alpha_pairs <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_substrate_pair_contrasts_v3.tsv"))
  alpha_season_pairs <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_season_pair_contrasts_v3.tsv"))
  alpha_ledger <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_independent_contrast_ledger_v3.tsv"))
  alpha_tests <- read_tsv(file.path(finals[["alpha"]], "results", "alpha_exact_signflip_results_v3.tsv"))
  assert_true(nrow(alpha_samples) == 900L, "Alpha sample row count differs from 900")
  assert_true(nrow(alpha_cells) == 300L, "Alpha cell summary row count differs from 300")
  assert_true(nrow(alpha_pairs) == 432L, "Alpha pair contrast row count differs from 432")
  assert_true(nrow(alpha_season_pairs) == 396L, "Alpha season-pair contrast row count differs from 396")
  assert_true(nrow(alpha_ledger) == 420L, "Alpha contrast ledger row count differs from 420")
  assert_true(nrow(alpha_tests) == 60L, "Alpha test row count differs from 60")
  assert_true(all(alpha_ledger$PairIDsUsed[alpha_ledger$ModelID %in%
    c("ALPHA_A72_BLOCK_CONTRAST", "ALPHA_A60_TYPE_BY_CAMPAIGN",
      "ALPHA_A36_SAND_CAMPAIGN", "ALPHA_A30_WATER_CAMPAIGN",
      "ALPHA_A42_FIRST_SAMPLETYPE")] %in% c(TRUE, "TRUE")),
    "Confirmed alpha paired models do not record PairID use")

  beta_root <- finals[["beta"]]
  beta_registry <- read_tsv(file.path(beta_root, "results", "beta_surface_registry_v4.tsv"))
  permanova <- read_tsv(file.path(beta_root, "results", "beta_permanova_restricted_v4.tsv"))
  permdisp <- read_tsv(file.path(beta_root, "results", "beta_permdisp_restricted_v4.tsv"))
  assert_true(nrow(beta_registry) == 24L, "Beta surface registry row count differs from 24")
  assert_true(nrow(permdisp) == 20L, "PERMDISP row count differs from 20")
  assert_true(all(permanova$PermutationsRequested == N_PERMUTATIONS), "PERMANOVA permutation count drift")
  assert_true(all(permdisp$PermutationsRequested == N_PERMUTATIONS), "PERMDISP permutation count drift")
  assert_true(all(permanova$PairIDsUsed[permanova$ModelID == "BETA_SAMPLETYPE"] %in% c(TRUE, "TRUE")),
              "A72 PERMANOVA does not record PairID use")
  assert_true(all(permanova$Strata[permanova$ModelID == "BETA_SAMPLETYPE"] == "PairID_Substrate"),
              "A72 PERMANOVA strata drift")
  season_models <- c("BETA_TYPE_CAMPAIGN", "BETA_SAND_CAMPAIGN", "BETA_WATER_CAMPAIGN")
  assert_true(all(permanova$PairIDsUsed[permanova$ModelID %in% season_models] %in% c(TRUE, "TRUE")),
              "Season-sensitive PERMANOVA rows do not record fixed-location pairing")
  assert_true(all(permanova$Strata[permanova$ModelID %in% season_models] == "LocationID"),
              "Season-sensitive PERMANOVA strata must be LocationID")
  assert_true(sum(permanova$InferenceReportable %in% c(TRUE, "TRUE")) == 20L,
              "PERMANOVA must have exactly one reportable target per marker x surface x model")
  cat(sprintf("PASS: %s\nPASS: %s\n", finals[["alpha"]], finals[["beta"]]))
}

dry_run_report <- function(inputs, base) {
  finals <- output_paths(base)
  cat("DRY-RUN PASS: all inputs, schemas, sample identities, depths, sets, and model rows validated.\n")
  cat(sprintf("Samples: %d\n", nrow(inputs$master)))
  cat(sprintf("Set counts: %s\n", paste(sprintf("%s=%d", names(inputs$set_counts), inputs$set_counts), collapse = ", ")))
  for (marker in MARKER_ORDER) for (surface in SURFACE_ORDER) {
    item <- inputs$counts[[marker]][[surface]]
    cat(sprintf("%s %s: %d samples x %d features; depth range %.0f-%.0f\n",
                marker, surface, nrow(item$matrix), ncol(item$matrix), min(item$depths), max(item$depths)))
  }
  cat(sprintf("Planned alpha output (must not exist): %s\n", finals[["alpha"]]))
  cat(sprintf("Planned beta output (must not exist): %s\n", finals[["beta"]]))
  if (any(dir.exists(finals) | file.exists(finals))) cat("NOTE: --run would refuse because a final output already exists.\n")
}

main <- function() {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (isTRUE(cli$help)) { usage(); return(invisible(0L)) }
  check_packages()
  if (cli$mode == "check") {
    check_existing(cli$base)
  } else {
    inputs <- load_inputs(cli$base)
    if (cli$mode == "dry-run") dry_run_report(inputs, cli$base) else build_all(inputs, cli$base)
  }
  invisible(0L)
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L, runLast = FALSE)
})
