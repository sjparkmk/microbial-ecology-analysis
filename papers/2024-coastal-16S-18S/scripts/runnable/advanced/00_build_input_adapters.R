#!/usr/bin/env Rscript

# Canonical 75-sample inputs -> EasyAmplicon-compatible adapter tables.
# This script never edits source tables and refuses to overwrite generated files.

options(stringsAsFactors = FALSE, warn = 1)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert <- function(ok, ...) if (!isTRUE(ok)) stopf(...)

script_path <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", a, value = TRUE)
  assert(length(hit) == 1L, "Could not determine script path")
  normalizePath(sub("^--file=", "", hit), winslash = "/", mustWork = TRUE)
}

arg_value <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE)
  i <- match(flag, a)
  if (is.na(i)) return(default)
  assert(i < length(a), "Missing value after %s", flag)
  a[[i + 1L]]
}

has_flag <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)

read_tsv <- function(path) {
  assert(file.exists(path), "Missing source input: %s", path)
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", comment.char = "", na.strings = c("NA"))
}

write_tsv_new <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  assert(!file.exists(path), "Refusing to overwrite: %s", path)
  tmp <- paste0(path, ".tmp")
  assert(!file.exists(tmp), "Temporary output already exists: %s", tmp)
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  assert(file.rename(tmp, path), "Atomic promote failed: %s", path)
  invisible(path)
}

copy_new <- function(source, destination) {
  assert(file.exists(source), "Missing configuration source: %s", source)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  assert(!file.exists(destination), "Refusing to overwrite: %s", destination)
  assert(file.copy(source, destination, overwrite = FALSE, copy.date = TRUE),
         "Could not copy %s to %s", source, destination)
  invisible(destination)
}

sha256 <- function(path) {
  assert(requireNamespace("digest", quietly = TRUE),
         "R package 'digest' is required for SHA-256 provenance")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

clean_taxon <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | x == "" | tolower(x) %in%
    c("na", "nan", "unknown", "unassigned", "unclassified")
  x[bad] <- "Unassigned"
  x
}

standardize_counts <- function(path, sample_ids, expected_depth = NA_integer_) {
  x <- read_tsv(path)
  assert(ncol(x) >= 2L, "Count table has fewer than two columns: %s", path)
  ids <- as.character(x[[1L]])
  assert(!anyNA(ids) && !any(ids == "") && !anyDuplicated(ids),
         "Invalid or duplicated feature IDs: %s", path)
  sample_cols <- names(x)[-1L]
  assert(setequal(sample_cols, sample_ids),
         "Count/design SampleID sets differ: %s", path)
  x <- x[, c(names(x)[1L], sample_ids), drop = FALSE]
  m <- as.matrix(x[, -1L, drop = FALSE])
  storage.mode(m) <- "numeric"
  assert(all(is.finite(m)) && all(m >= 0) && all(m == floor(m)),
         "Counts must be finite, nonnegative integers: %s", path)
  assert(all(rowSums(m) > 0), "Zero-sum features remain in analysis table: %s", path)
  assert(all(colSums(m) > 0), "Zero-depth samples found: %s", path)
  if (is.finite(expected_depth)) {
    assert(all(colSums(m) == expected_depth),
           "Expected every sample depth to equal %s: %s",
           format(expected_depth, scientific = FALSE), path)
  }
  out <- data.frame(feature_id = ids, m, check.names = FALSE)
  names(out)[1L] <- "#OTUID"
  list(table = out, matrix = m, feature_ids = ids,
       min_depth = min(colSums(m)), max_depth = max(colSums(m)))
}

align_rows <- function(x, id_col, ids, label) {
  assert(id_col %in% names(x), "%s ID column missing: %s", label, id_col)
  key <- as.character(x[[id_col]])
  assert(!anyNA(key) && !anyDuplicated(key), "%s contains invalid/duplicate IDs", label)
  idx <- match(ids, key)
  assert(!anyNA(idx), "%s is missing %d count-table feature IDs", label, sum(is.na(idx)))
  x[idx, , drop = FALSE]
}

easy_taxonomy_16s <- function(path, ids) {
  x <- align_rows(read_tsv(path), "OTUID", ids, "16S taxonomy")
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  assert(all(ranks %in% names(x)), "16S taxonomy is missing canonical rank columns")
  out <- x[, c("OTUID", ranks), drop = FALSE]
  out[ranks] <- lapply(out[ranks], clean_taxon)
  out
}

easy_taxonomy_18s <- function(path, ids) {
  x <- align_rows(read_tsv(path), "ASV_ID", ids, "18S PR2-primary taxonomy")
  source_cols <- c("Kingdom", "Division_or_Phylum", "Class", "Order",
                   "Family", "Genus", "Species")
  assert(all(c("Domain", source_cols) %in% names(x)),
         "18S PR2-primary taxonomy is missing canonical columns")
  out <- x[, c("ASV_ID", source_cols), drop = FALSE]
  names(out) <- c("OTUID", "Kingdom", "Phylum", "Class", "Order",
                  "Family", "Genus", "Species")
  out[-1L] <- lapply(out[-1L], clean_taxon)
  out
}

preserve_primary_18s <- function(path, ids) {
  align_rows(read_tsv(path), "ASV_ID", ids, "18S PR2-primary taxonomy")
}

preserve_parallel_18s <- function(path, ids) {
  x <- align_rows(read_tsv(path), "ASV_ID", ids,
                  "18S PR2/SILVA parallel taxonomy")
  required <- c(
    "pr2_source_database", "pr2_raw_assignment", "pr2_canonical_domain",
    "pr2_canonical_kingdom", "pr2_canonical_division_or_phylum",
    "pr2_canonical_class", "pr2_canonical_order", "pr2_canonical_family",
    "pr2_canonical_genus", "pr2_canonical_species", "silva_source_database",
    "silva_raw_assignment", "silva_canonical_domain", "silva_canonical_kingdom",
    "silva_canonical_division_or_phylum", "silva_canonical_class",
    "silva_canonical_order", "silva_canonical_family", "silva_canonical_genus",
    "silva_canonical_species", "agreement_level", "agreement_status",
    "analysis_keep", "exclusion_reason"
  )
  miss <- setdiff(required, names(x))
  assert(length(miss) == 0L,
         "Parallel taxonomy lost required columns: %s", paste(miss, collapse = ", "))
  x
}

adapter_metadata <- function(design) {
  required <- c(
    "SampleID", "Site", "GeologicalHabitat", "SampleType", "Season", "Campaign",
    "CampaignOriginal", "FieldReplicate", "LocationID", "PairID_Substrate",
    "PairID_Season", "SiteSeasonBlock", "EnvironmentRecordID",
    "EnvironmentEligible", "ReplicateStatus", "ReplicateSpacingM", "SamplingMethod"
  )
  miss <- setdiff(required, names(design))
  assert(length(miss) == 0L, "Design columns missing: %s", paste(miss, collapse = ", "))
  assert(nrow(design) == 75L && !anyDuplicated(design$SampleID),
         "Design must contain exactly 75 unique samples")
  assert(setequal(unique(design$Season), c("Spring", "Autumn")),
         "Season must be Spring/Autumn")
  assert(setequal(unique(design$SampleType), c("Sea-sand", "Water")),
         "SampleType must be Sea-sand/Water")

  data.frame(
    SampleID = design$SampleID,
    Group = design$GeologicalHabitat,
    group = design$GeologicalHabitat,
    site = design$Site,
    type = design$SampleType,
    time = design$Season,
    timepoint = design$Season,
    subID = design$EnvironmentRecordID,
    Site = design$Site,
    GeologicalHabitat = design$GeologicalHabitat,
    SampleType = design$SampleType,
    Season = design$Season,
    Campaign = design$Campaign,
    CampaignOriginal = design$CampaignOriginal,
    FieldReplicate = design$FieldReplicate,
    LocationID = design$LocationID,
    PairID_Substrate = design$PairID_Substrate,
    PairID_Season = design$PairID_Season,
    SiteSeasonBlock = design$SiteSeasonBlock,
    EnvironmentRecordID = design$EnvironmentRecordID,
    EnvironmentEligible = design$EnvironmentEligible,
    ReplicateStatus = design$ReplicateStatus,
    ReplicateSpacingM = design$ReplicateSpacingM,
    SamplingMethod = design$SamplingMethod,
    check.names = FALSE
  )
}

metadata_dictionary <- function() {
  data.frame(
    adapter_column = c(
      "SampleID", "Group", "group", "site", "type", "time", "timepoint", "subID",
      "Site", "GeologicalHabitat", "SampleType", "Season", "Campaign",
      "CampaignOriginal", "FieldReplicate", "LocationID", "PairID_Substrate",
      "PairID_Season", "SiteSeasonBlock", "EnvironmentRecordID",
      "EnvironmentEligible", "ReplicateStatus", "ReplicateSpacingM", "SamplingMethod"
    ),
    canonical_source = c(
      "SampleID", "GeologicalHabitat", "GeologicalHabitat", "Site", "SampleType",
      "Season", "Season", "EnvironmentRecordID", "Site", "GeologicalHabitat",
      "SampleType", "Season", "Campaign", "CampaignOriginal", "FieldReplicate",
      "LocationID", "PairID_Substrate", "PairID_Season", "SiteSeasonBlock",
      "EnvironmentRecordID", "EnvironmentEligible", "ReplicateStatus",
      "ReplicateSpacingM", "SamplingMethod"
    ),
    role = c(
      "sample_key", "legacy_alias", "legacy_alias", "legacy_alias", "legacy_alias",
      "legacy_alias", "legacy_alias", "environment_join_alias", rep("canonical", 16)
    ),
    semantics = c(
      "Unique sample identifier", "Ocean/Tidal/Island", "Ocean/Tidal/Island",
      "Seven named sites", "Sea-sand/Water; never soil/sediment in revised outputs",
      "Spring/Autumn; never First/Second in revised outputs",
      "Spring/Autumn; never F/S in revised outputs",
      "Explicit environment-record key; only eligible Spring Sea-sand rows enter RDA/Mantel",
      "Seven named sites", "Ocean/Tidal/Island", "Sea-sand/Water", "Spring/Autumn",
      "Spring/Autumn reader-facing campaign", "Raw First/Second provenance only",
      "Spatial bulk-harvest replicate 1-3", "Fixed field location across seasons",
      "Contemporaneous Sea-sand/Water pair", "Spring/Autumn repeated-location pair",
      "Site x season aggregation block", "Environment table join key",
      "Eligibility flag; analysis still enforces Spring Sea-sand restriction",
      "Biological spatial field replicate", "Approximately 0.5-1 m", "Bulk harvest"
    ),
    stringsAsFactors = FALSE
  )
}

this_script <- script_path()
script_dir <- dirname(this_script)
default_root <- normalizePath(file.path(script_dir, "..", "..", ".."),
                              winslash = "/", mustWork = TRUE)
project_root <- normalizePath(arg_value("--project-root", default_root),
                              winslash = "/", mustWork = TRUE)
output_root <- arg_value(
  "--output-root",
  file.path(project_root, "revision", "re_analysis",
            "12_original_pipeline_reanalysis_v1")
)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
check_only <- has_flag("--check-only")

expected_output_root <- normalizePath(
  file.path(project_root, "revision", "re_analysis",
            "12_original_pipeline_reanalysis_v1"),
  winslash = "/", mustWork = FALSE
)
assert(identical(tolower(output_root), tolower(expected_output_root)),
       "Output root must be the frozen v1 result directory: %s", expected_output_root)

source_paths <- list(
  design = file.path(project_root, "revision/re_analysis/00_design_v4/master_sample_metadata_v4.tsv"),
  counts16_raw = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_asv_counts_analysis_75_v1.tsv"),
  tax16_raw = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_taxonomy_analysis_75_v1.tsv"),
  counts16_rare = file.path(project_root, "revision/re_analysis/04_input_16s_rarefied_20000_v1/results/16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  counts18_raw = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_asv_counts_pr2_analysis_v1.tsv"),
  tax18_primary_raw = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_taxonomy_pr2_primary_analysis_v1.tsv"),
  tax18_parallel_raw = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_taxonomy_parallel_pr2_silva123_all_v2.tsv"),
  counts18_rare = file.path(project_root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
  tax18_primary_rare = file.path(project_root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_taxonomy_pr2_primary_rarefied_50000_analysis_v1.tsv"),
  tax18_parallel_rare = file.path(project_root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_taxonomy_parallel_pr2_silva123_rarefied_50000_analysis_v2.tsv")
)
invisible(lapply(source_paths, function(p) assert(file.exists(p), "Missing input: %s", p)))

design <- read_tsv(source_paths$design)
metadata <- adapter_metadata(design)
sample_ids <- metadata$SampleID

c16_raw <- standardize_counts(source_paths$counts16_raw, sample_ids)
c16_rare <- standardize_counts(source_paths$counts16_rare, sample_ids, 20000L)
c18_raw <- standardize_counts(source_paths$counts18_raw, sample_ids)
c18_rare <- standardize_counts(source_paths$counts18_rare, sample_ids, 50000L)

assert(identical(c16_raw$feature_ids, easy_taxonomy_16s(source_paths$tax16_raw,
                                                       c16_raw$feature_ids)$OTUID),
       "16S raw taxonomy order/linkage check failed")
tax16 <- easy_taxonomy_16s(source_paths$tax16_raw, c16_raw$feature_ids)
tax18_easy <- easy_taxonomy_18s(source_paths$tax18_primary_raw, c18_raw$feature_ids)
tax18_primary_raw <- preserve_primary_18s(source_paths$tax18_primary_raw,
                                          c18_raw$feature_ids)
tax18_primary_rare <- preserve_primary_18s(source_paths$tax18_primary_rare,
                                           c18_rare$feature_ids)
tax18_parallel_raw <- preserve_parallel_18s(source_paths$tax18_parallel_raw,
                                            c18_raw$feature_ids)
tax18_parallel_rare <- preserve_parallel_18s(source_paths$tax18_parallel_rare,
                                             c18_rare$feature_ids)

assert(all(c16_rare$feature_ids %in% c16_raw$feature_ids),
       "16S rarefied feature IDs are not a subset of raw 16S IDs")
assert(all(c18_rare$feature_ids %in% c18_raw$feature_ids),
       "18S rarefied feature IDs are not a subset of raw PR2-primary IDs")

cat(sprintf("Validated 75 samples; 16S raw/rarefied features=%d/%d, 18S raw/rarefied features=%d/%d\n",
            length(c16_raw$feature_ids), length(c16_rare$feature_ids),
            length(c18_raw$feature_ids), length(c18_rare$feature_ids)))
cat(sprintf("Rarefied depths: 16S=%d; 18S=%d reads/sample\n",
            c16_rare$min_depth, c18_rare$min_depth))

if (check_only) {
  cat("Check-only validation completed; no files written.\n")
  quit(save = "no", status = 0L)
}

input_dir <- file.path(output_root, "inputs")
registry_dir <- file.path(output_root, "registry")

outputs <- list(
  metadata = list(
    path = file.path(input_dir, "metadata_easyamplicon_75_v1.tsv"),
    table = metadata, marker = "both", dataset = "metadata_adapter",
    source_key = "design"
  ),
  metadata_dictionary = list(
    path = file.path(input_dir, "metadata_easyamplicon_dictionary_v1.tsv"),
    table = metadata_dictionary(), marker = "both", dataset = "metadata_dictionary",
    source_key = "design"
  ),
  counts16_raw = list(
    path = file.path(input_dir, "16s_otutab_unrarefied_75_v1.tsv"),
    table = c16_raw$table, marker = "16S", dataset = "unrarefied_counts",
    source_key = "counts16_raw"
  ),
  counts16_rare = list(
    path = file.path(input_dir, "16s_otutab_rarefied_20000_75_v1.tsv"),
    table = c16_rare$table, marker = "16S", dataset = "rarefied_counts_20000",
    source_key = "counts16_rare"
  ),
  tax16 = list(
    path = file.path(input_dir, "16s_taxonomy_easyamplicon_v1.tsv"),
    table = tax16, marker = "16S", dataset = "taxonomy_easyamplicon",
    source_key = "tax16_raw"
  ),
  counts18_raw = list(
    path = file.path(input_dir, "18s_otutab_unrarefied_pr2_75_v1.tsv"),
    table = c18_raw$table, marker = "18S", dataset = "unrarefied_counts_pr2",
    source_key = "counts18_raw"
  ),
  counts18_rare = list(
    path = file.path(input_dir, "18s_otutab_rarefied_50000_pr2_75_v1.tsv"),
    table = c18_rare$table, marker = "18S", dataset = "rarefied_counts_pr2_50000",
    source_key = "counts18_rare"
  ),
  tax18_easy = list(
    path = file.path(input_dir, "18s_taxonomy_pr2_easyamplicon_v1.tsv"),
    table = tax18_easy, marker = "18S", dataset = "taxonomy_pr2_easyamplicon",
    source_key = "tax18_primary_raw"
  ),
  tax18_primary_raw = list(
    path = file.path(input_dir, "18s_taxonomy_pr2_primary_unrarefied_preserved_v1.tsv"),
    table = tax18_primary_raw, marker = "18S", dataset = "taxonomy_pr2_primary_preserved",
    source_key = "tax18_primary_raw"
  ),
  tax18_primary_rare = list(
    path = file.path(input_dir, "18s_taxonomy_pr2_primary_rarefied_50000_preserved_v1.tsv"),
    table = tax18_primary_rare, marker = "18S", dataset = "taxonomy_pr2_primary_rarefied_preserved",
    source_key = "tax18_primary_rare"
  ),
  tax18_parallel_raw = list(
    path = file.path(input_dir, "18s_taxonomy_parallel_pr2_silva123_unrarefied_full_v1.tsv"),
    table = tax18_parallel_raw, marker = "18S", dataset = "taxonomy_parallel_full",
    source_key = "tax18_parallel_raw"
  ),
  tax18_parallel_rare = list(
    path = file.path(input_dir, "18s_taxonomy_parallel_pr2_silva123_rarefied_50000_full_v1.tsv"),
    table = tax18_parallel_rare, marker = "18S", dataset = "taxonomy_parallel_rarefied_full",
    source_key = "tax18_parallel_rare"
  )
)

for (z in outputs) write_tsv_new(z$table, z$path)

config_files <- c(
  "original_script_crosswalk_v1.tsv",
  "module_decision_registry_v1.tsv",
  "module_status_categories_v1.tsv",
  "module_execution_manifest_v1.tsv"
)
for (f in config_files) {
  copy_new(file.path(script_dir, "config", f), file.path(registry_dir, f))
}

manifest <- do.call(rbind, lapply(names(outputs), function(nm) {
  z <- outputs[[nm]]
  source <- source_paths[[z$source_key]]
  data.frame(
    adapter_id = nm,
    marker = z$marker,
    dataset = z$dataset,
    source_path = normalizePath(source, winslash = "/", mustWork = TRUE),
    source_sha256 = sha256(source),
    output_path = normalizePath(z$path, winslash = "/", mustWork = TRUE),
    output_sha256 = sha256(z$path),
    rows = nrow(z$table),
    columns = ncol(z$table),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv_new(manifest, file.path(registry_dir, "input_adapter_manifest_v1.tsv"))

depth_qc <- data.frame(
  marker = c("16S", "16S", "18S", "18S"),
  dataset = c("unrarefied", "rarefied_20000", "unrarefied_pr2", "rarefied_50000_pr2"),
  samples = 75L,
  features = c(length(c16_raw$feature_ids), length(c16_rare$feature_ids),
               length(c18_raw$feature_ids), length(c18_rare$feature_ids)),
  minimum_depth = c(c16_raw$min_depth, c16_rare$min_depth,
                    c18_raw$min_depth, c18_rare$min_depth),
  maximum_depth = c(c16_raw$max_depth, c16_rare$max_depth,
                    c18_raw$max_depth, c18_rare$max_depth),
  expected_rarefied_depth = c(NA_integer_, 20000L, NA_integer_, 50000L),
  depth_contract_pass = c(TRUE, c16_rare$min_depth == 20000L && c16_rare$max_depth == 20000L,
                          TRUE, c18_rare$min_depth == 50000L && c18_rare$max_depth == 50000L),
  stringsAsFactors = FALSE
)
write_tsv_new(depth_qc, file.path(registry_dir, "input_depth_qc_v1.tsv"))

cat(sprintf("Wrote %d immutable adapter tables to %s\n", length(outputs), input_dir))
cat(sprintf("Wrote registry/configuration artifacts to %s\n", registry_dir))
