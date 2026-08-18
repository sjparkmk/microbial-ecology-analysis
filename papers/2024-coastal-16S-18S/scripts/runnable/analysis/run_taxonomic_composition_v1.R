#!/usr/bin/env Rscript

# Build taxonomic-composition source tables for the frozen 75-sample revision.
#
# Design and input packages are immutable.  --run writes only to a uniquely
# named sibling staging directory and publishes it by rename after every
# invariant and schema check passes.  This script intentionally makes no
# figures and performs no differential-abundance testing.

SCRIPT_VERSION <- "1.0.0"
PACKAGE_ID <- "07_taxonomic_composition_v1"
EXPECTED_SAMPLES <- 75L
SUM_TOLERANCE <- 1e-10
TOP_N <- 12L

usage <- function(status = 0L) {
  message(paste(
    "Usage:",
    "  Rscript run_taxonomic_composition_v1.R --dry-run [--root PATH]",
    "  Rscript run_taxonomic_composition_v1.R --run [--root PATH]",
    "  Rscript run_taxonomic_composition_v1.R --check [--root PATH]",
    "",
    "Exactly one mode is required.  PATH defaults to revision/re_analysis.",
    "The v1 plotting backend is separate; this script writes tables only.",
    sep = "\n"
  ))
  quit(save = "no", status = status, runLast = FALSE)
}

fail <- function(message) stop(message, call. = FALSE)
require_true <- function(condition, message) {
  if (!isTRUE(condition)) fail(message)
}

script_path <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  require_true(length(arg) == 1L, "Cannot resolve the script path")
  normalizePath(sub("^--file=", "", arg), winslash = "/", mustWork = TRUE)
}

parse_args <- function(args) {
  if (length(args) == 0L) usage(2L)
  if ("--help" %in% args || "-h" %in% args) usage(0L)
  out <- list(run = FALSE, dry_run = FALSE, check = FALSE, root = NULL)
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (token %in% c("--run", "--dry-run", "--check")) {
      if (token == "--run") out$run <- TRUE
      if (token == "--dry-run") out$dry_run <- TRUE
      if (token == "--check") out$check <- TRUE
      i <- i + 1L
      next
    }
    require_true(token == "--root", paste0("Unknown argument: ", token))
    require_true(i < length(args), "Missing value for --root")
    out$root <- args[[i + 1L]]
    i <- i + 2L
  }
  require_true(sum(c(out$run, out$dry_run, out$check)) == 1L,
               "Select exactly one of --run, --dry-run, or --check")
  if (is.null(out$root)) {
    out$root <- file.path(dirname(script_path()), "..")
  }
  out$root <- normalizePath(out$root, winslash = "/", mustWork = TRUE)
  out
}

input_paths <- function(root) {
  list(
    metadata = file.path(root, "00_design_v2", "master_sample_metadata_v2.tsv"),
    membership = file.path(root, "00_design_v2", "analysis_set_membership_v2.tsv"),
    counts16 = file.path(root, "04_input_16s_rarefied_20000_v1", "results",
                         "16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
    taxonomy16 = file.path(root, "04_input_16s_rarefied_20000_v1", "results",
                           "16s_taxonomy_rarefied_20000_analysis_v1.tsv"),
    counts18 = file.path(root, "03_18s_rarefied_50000_v1", "results",
                         "18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
    taxonomy18 = file.path(root, "03_18s_rarefied_50000_v1", "results",
                           "18s_taxonomy_pr2_primary_rarefied_50000_analysis_v1.tsv"),
    parallel18 = file.path(root, "03_18s_rarefied_50000_v1", "results",
                           "18s_taxonomy_parallel_pr2_silva123_rarefied_50000_analysis_v2.tsv")
  )
}

read_tsv <- function(path) {
  require_true(file.exists(path) && !dir.exists(path), paste0("Missing input: ", path))
  read.delim(path, sep = "\t", quote = "", comment.char = "", header = TRUE,
             check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = character(0), fileEncoding = "UTF-8")
}

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  require_true(length(missing) == 0L,
               paste0(label, " is missing columns: ", paste(missing, collapse = ", ")))
}

as_flag <- function(x, label) {
  y <- tolower(trimws(as.character(x)))
  require_true(all(y %in% c("true", "false")),
               paste0(label, " must contain only true/false"))
  y == "true"
}

read_counts <- function(path, marker) {
  x <- read_tsv(path)
  require_true(ncol(x) == EXPECTED_SAMPLES + 1L,
               paste0(marker, " count table must have 75 sample columns"))
  id_header <- names(x)[[1L]]
  require_true(id_header %in% c("#OTUID", "OTUID", "ASV_ID"),
               paste0(marker, " has unexpected ID header: ", id_header))
  ids <- as.character(x[[1L]])
  sample_ids <- names(x)[-1L]
  require_true(all(nzchar(ids)) && !anyDuplicated(ids),
               paste0(marker, " feature IDs are empty or duplicated"))
  require_true(all(nzchar(sample_ids)) && !anyDuplicated(sample_ids),
               paste0(marker, " sample IDs are empty or duplicated"))
  require_true(!any(grepl("^SS_BS(?:_|$)", sample_ids, perl = TRUE)),
               paste0(marker, " contains forbidden SS_BS samples"))
  raw <- unlist(x[-1L], use.names = FALSE)
  values <- suppressWarnings(as.numeric(raw))
  require_true(length(values) == length(raw) && !anyNA(values) && all(is.finite(values)),
               paste0(marker, " count table has non-numeric or non-finite cells"))
  require_true(all(values >= 0) && all(values == floor(values)),
               paste0(marker, " count table has negative or non-integer cells"))
  mat <- matrix(values, nrow = nrow(x), ncol = ncol(x) - 1L,
                dimnames = list(ids, sample_ids))
  list(id_header = id_header, counts = mat)
}

normalize_taxon <- function(x) {
  y <- trimws(as.character(x))
  bad <- is.na(y) | !nzchar(y) |
    tolower(y) %in% c("unassigned", "unknown", "na", "n/a", "not_assigned")
  y[bad] <- "Unassigned"
  y
}

sample_relative_abundance <- function(counts, label) {
  totals <- colSums(counts)
  require_true(all(totals > 0), paste0(label, " contains a zero-depth sample"))
  out <- sweep(counts, 2L, totals, "/")
  require_true(max(abs(colSums(out) - 1)) <= SUM_TOLERANCE,
               paste0(label, " sample-wise relative abundances do not sum to one"))
  out
}

aggregate_taxon <- function(sample_ra, taxonomy, id_column, rank, label) {
  require_columns(taxonomy, c(id_column, rank), label)
  idx <- match(rownames(sample_ra), taxonomy[[id_column]])
  require_true(!anyNA(idx), paste0(label, " does not cover every count-table feature"))
  taxa <- normalize_taxon(taxonomy[[rank]][idx])
  out <- rowsum(sample_ra, group = taxa, reorder = FALSE)
  require_true(max(abs(colSums(out) - 1)) <= SUM_TOLERANCE,
               paste0(label, " ", rank, " aggregation does not sum to one"))
  out
}

metadata_for_samples <- function(metadata, sample_ids) {
  idx <- match(sample_ids, metadata$SampleID)
  require_true(!anyNA(idx), "Metadata does not cover every sample")
  metadata[idx, , drop = FALSE]
}

make_cell_id <- function(metadata) {
  paste(metadata$SiteCodeOriginal, metadata$SampleType, metadata$Campaign, sep = "__")
}

sample_long <- function(taxon_mat, metadata, marker, dataset, rank) {
  md <- metadata_for_samples(metadata, colnames(taxon_mat))
  cell_id <- make_cell_id(md)
  n_taxa <- nrow(taxon_mat)
  out <- data.frame(
    Marker = marker,
    Dataset = dataset,
    Rank = rank,
    SampleID = rep(colnames(taxon_mat), each = n_taxa),
    Site = rep(md$Site, each = n_taxa),
    SiteCodeOriginal = rep(md$SiteCodeOriginal, each = n_taxa),
    GeologicalHabitat = rep(md$GeologicalHabitat, each = n_taxa),
    SampleType = rep(md$SampleType, each = n_taxa),
    Campaign = rep(md$Campaign, each = n_taxa),
    CellID = rep(cell_id, each = n_taxa),
    Taxon = rep(rownames(taxon_mat), times = ncol(taxon_mat)),
    RelativeAbundance = as.vector(taxon_mat),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out
}

cell_summary <- function(long) {
  keys <- unique(long[c("CellID", "Site", "SiteCodeOriginal", "GeologicalHabitat",
                        "SampleType", "Campaign")])
  stat <- aggregate(
    long$RelativeAbundance,
    by = list(CellID = long$CellID, Taxon = long$Taxon),
    FUN = function(z) c(mean = mean(z), sd = stats::sd(z), n = length(z))
  )
  vals <- as.data.frame(stat$x, stringsAsFactors = FALSE)
  names(vals) <- c("MeanRelativeAbundance", "SDRelativeAbundance", "N")
  stat$x <- NULL
  out <- cbind(stat, keys[match(stat$CellID, keys$CellID), setdiff(names(keys), "CellID"),
                          drop = FALSE], vals, stringsAsFactors = FALSE)
  out <- out[c("CellID", "Site", "SiteCodeOriginal", "GeologicalHabitat",
               "SampleType", "Campaign", "Taxon", "MeanRelativeAbundance",
               "SDRelativeAbundance", "N")]
  out[order(out$CellID, out$Taxon), , drop = FALSE]
}

decorate_cell_summary <- function(cell, marker, dataset, rank) {
  data.frame(Marker = marker, Dataset = dataset, Rank = rank, cell,
             stringsAsFactors = FALSE, check.names = FALSE)
}

safe_slug <- function(rank) tolower(rank)

select_top_taxa <- function(taxon_mat, n = TOP_N) {
  assigned <- setdiff(rownames(taxon_mat), "Unassigned")
  means <- rowMeans(taxon_mat[assigned, , drop = FALSE])
  ordered <- assigned[order(-means, assigned)]
  head(ordered, n)
}

collapse_top <- function(taxon_mat, top_taxa) {
  group <- ifelse(rownames(taxon_mat) == "Unassigned", "Unassigned",
                  ifelse(rownames(taxon_mat) %in% top_taxa,
                         rownames(taxon_mat), "Other"))
  out <- rowsum(taxon_mat, group = group, reorder = FALSE)
  for (required in c("Other", "Unassigned")) {
    if (!required %in% rownames(out)) {
      zero <- matrix(0, nrow = 1L, ncol = ncol(out),
                     dimnames = list(required, colnames(out)))
      out <- rbind(out, zero)
    }
  }
  order_names <- c(top_taxa[top_taxa %in% rownames(out)], "Other", "Unassigned")
  out <- out[order_names, , drop = FALSE]
  require_true(max(abs(colSums(out) - 1)) <= SUM_TOLERANCE,
               "Top-taxon collapse does not sum to one")
  out
}

top_sample_table <- function(top_mat, metadata, marker, dataset, rank) {
  x <- sample_long(top_mat, metadata, marker, dataset, rank)
  names(x)[names(x) == "Taxon"] <- "TaxonDisplay"
  order_map <- setNames(seq_len(nrow(top_mat)), rownames(top_mat))
  x$CategoryOrder <- unname(order_map[x$TaxonDisplay])
  x[c("Marker", "Dataset", "Rank", "SampleID", "Site", "SiteCodeOriginal",
      "GeologicalHabitat", "SampleType", "Campaign", "CellID", "TaxonDisplay",
      "CategoryOrder", "RelativeAbundance")]
}

top_cell_table <- function(top_sample) {
  tmp <- top_sample
  names(tmp)[names(tmp) == "TaxonDisplay"] <- "Taxon"
  names(tmp)[names(tmp) == "RelativeAbundance"] <- "RelativeAbundance"
  cell <- cell_summary(tmp)
  names(cell)[names(cell) == "Taxon"] <- "TaxonDisplay"
  order_map <- unique(top_sample[c("TaxonDisplay", "CategoryOrder")])
  cell$CategoryOrder <- order_map$CategoryOrder[match(cell$TaxonDisplay,
                                                     order_map$TaxonDisplay)]
  data.frame(
    Marker = unique(top_sample$Marker),
    Dataset = unique(top_sample$Dataset),
    Rank = unique(top_sample$Rank),
    cell[c("CellID", "Site", "SiteCodeOriginal", "GeologicalHabitat",
           "SampleType", "Campaign", "TaxonDisplay")],
    CategoryOrder = cell$CategoryOrder,
    cell[c("MeanRelativeAbundance", "SDRelativeAbundance", "N")],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_tsv <- function(x, path) {
  require_true(!file.exists(path), paste0("Refusing to overwrite staged file: ", path))
  require_true(dir.exists(dirname(path)), paste0("Missing output directory: ", dirname(path)))
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "", eol = "\n", fileEncoding = "UTF-8")
}

write_json <- function(x, path) {
  require_true(!file.exists(path), paste0("Refusing to overwrite staged file: ", path))
  jsonlite::write_json(x, path = path, pretty = TRUE, auto_unbox = TRUE,
                       null = "null", na = "null")
  cat("\n", file = path, append = TRUE)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

relative_path <- function(path, root) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  r <- normalizePath(root, winslash = "/", mustWork = TRUE)
  require_true(startsWith(p, paste0(r, "/")), "Manifest path escaped package root")
  substring(p, nchar(r) + 2L)
}

write_manifest <- function(stage) {
  manifest_path <- file.path(stage, "provenance", "sha256_manifest_taxonomic_composition_v1.tsv")
  files <- list.files(stage, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- setdiff(normalizePath(files, winslash = "/", mustWork = TRUE),
                   normalizePath(manifest_path, winslash = "/", mustWork = FALSE))
  files <- sort(files)
  manifest <- data.frame(
    package = PACKAGE_ID,
    path = vapply(files, relative_path, character(1), root = stage),
    bytes = as.numeric(file.info(files)$size),
    sha256 = vapply(files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, manifest_path)
}

validate_sources <- function(paths) {
  invisible(lapply(paths, function(p) {
    require_true(file.exists(p) && !dir.exists(p), paste0("Missing input: ", p))
  }))
  metadata <- read_tsv(paths$metadata)
  require_columns(metadata,
                  c("SampleID", "Site", "SiteCodeOriginal", "GeologicalHabitat",
                    "SampleType", "Campaign", "Included16S", "Included18S"),
                  "design v2 metadata")
  require_true(nrow(metadata) == EXPECTED_SAMPLES && !anyDuplicated(metadata$SampleID),
               "Design v2 metadata must contain 75 unique SampleIDs")
  require_true(all(as_flag(metadata$Included16S, "Included16S")) &&
                 all(as_flag(metadata$Included18S, "Included18S")),
               "All design-v2 samples must be shared 16S/18S inclusions")
  require_true(all(metadata$SampleType %in% c("Sea-sand", "Water")),
               "Unexpected SampleType in design v2")
  require_true(all(metadata$Campaign %in% c("First", "Second")),
               "Unexpected Campaign in design v2")
  require_true(!any(grepl("^SS_BS(?:_|$)", metadata$SampleID, perl = TRUE)),
               "Design v2 contains forbidden SS_BS samples")

  membership <- read_tsv(paths$membership)
  require_columns(membership, c("SetID", "SampleID", "Included"), "design membership v2")
  a75_rows <- membership$SetID == "A75_FULL"
  a75 <- membership$SampleID[a75_rows][
    as_flag(membership$Included[a75_rows], "A75 Included")
  ]
  require_true(length(a75) == EXPECTED_SAMPLES && setequal(a75, metadata$SampleID),
               "A75_FULL membership does not exactly match design metadata")

  c16 <- read_counts(paths$counts16, "16S")
  c18 <- read_counts(paths$counts18, "18S")
  require_true(setequal(colnames(c16$counts), metadata$SampleID),
               "16S samples do not exactly match design v2")
  require_true(setequal(colnames(c18$counts), metadata$SampleID),
               "18S samples do not exactly match design v2")
  require_true(setequal(colnames(c16$counts), colnames(c18$counts)),
               "16S and 18S sample sets differ")
  c16$counts <- c16$counts[, metadata$SampleID, drop = FALSE]
  c18$counts <- c18$counts[, metadata$SampleID, drop = FALSE]
  require_true(all(colSums(c16$counts) == 20000),
               "16S analysis counts are not exactly 20,000 reads per sample")
  require_true(all(colSums(c18$counts) == 50000),
               "18S analysis counts are not exactly 50,000 reads per sample")

  t16 <- read_tsv(paths$taxonomy16)
  require_columns(t16, c("OTUID", "Phylum", "Class", "Genus"), "16S taxonomy")
  require_true(!anyDuplicated(t16$OTUID) && setequal(t16$OTUID, rownames(c16$counts)),
               "16S taxonomy IDs do not exactly match counts")

  t18 <- read_tsv(paths$taxonomy18)
  require_columns(t18, c("ASV_ID", "Kingdom", "Division_or_Phylum", "Class", "Genus"),
                  "18S PR2-primary taxonomy")
  require_true(!anyDuplicated(t18$ASV_ID) && setequal(t18$ASV_ID, rownames(c18$counts)),
               "18S PR2-primary taxonomy IDs do not exactly match counts")

  p18 <- read_tsv(paths$parallel18)
  require_columns(p18,
                  c("ASV_ID", "pr2_canonical_kingdom", "pr2_canonical_division_or_phylum",
                    "pr2_canonical_class", "pr2_canonical_genus", "agreement_level",
                    "agreement_status", "analysis_keep", "flag_pr2_metazoa"),
                  "18S PR2/SILVA parallel taxonomy v2")
  require_true(!anyDuplicated(p18$ASV_ID) && setequal(p18$ASV_ID, rownames(c18$counts)),
               "18S parallel-v2 IDs do not exactly match counts")
  pidx <- match(t18$ASV_ID, p18$ASV_ID)
  require_true(all(as_flag(p18$analysis_keep[pidx], "parallel analysis_keep")),
               "Parallel-v2 taxonomy contains a non-analysis feature")
  comparisons <- list(
    Kingdom = "pr2_canonical_kingdom",
    Division_or_Phylum = "pr2_canonical_division_or_phylum",
    Class = "pr2_canonical_class",
    Genus = "pr2_canonical_genus"
  )
  for (rank in names(comparisons)) {
    a <- normalize_taxon(t18[[rank]])
    b <- normalize_taxon(p18[[comparisons[[rank]]]][pidx])
    require_true(identical(a, b), paste0("PR2-primary and parallel-v2 disagree at ", rank))
  }
  flags <- as_flag(p18$flag_pr2_metazoa, "flag_pr2_metazoa")
  require_true(any(flags) && any(!flags),
               "Exact PR2 Metazoa flag must contain both flagged and retained features")

  list(metadata = metadata, membership = membership, counts16 = c16$counts,
       counts18 = c18$counts, taxonomy16 = t16, taxonomy18 = t18,
       parallel18 = p18, metazoa_flags = setNames(flags, p18$ASV_ID))
}

composition_invariant <- function(marker, dataset, rank, sample_long_table, cell_table,
                                  feature_count, taxon_count) {
  sample_sums <- aggregate(RelativeAbundance ~ SampleID, sample_long_table, sum)
  cell_sums <- aggregate(MeanRelativeAbundance ~ CellID, cell_table, sum)
  ok <- max(abs(sample_sums$RelativeAbundance - 1)) <= SUM_TOLERANCE &&
    max(abs(cell_sums$MeanRelativeAbundance - 1)) <= SUM_TOLERANCE
  data.frame(
    Marker = marker, Dataset = dataset, Rank = rank,
    SourceFeatureCount = feature_count, TaxonCount = taxon_count,
    SampleCount = nrow(sample_sums), CellCount = nrow(cell_sums),
    MinimumSampleSum = min(sample_sums$RelativeAbundance),
    MaximumSampleSum = max(sample_sums$RelativeAbundance),
    MinimumCellMeanSum = min(cell_sums$MeanRelativeAbundance),
    MaximumCellMeanSum = max(cell_sums$MeanRelativeAbundance),
    Tolerance = SUM_TOLERANCE,
    Status = if (ok) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}

top_invariant <- function(sample_table, cell_table) {
  ss <- aggregate(RelativeAbundance ~ Marker + Dataset + Rank + SampleID,
                  sample_table, sum)
  cs <- aggregate(MeanRelativeAbundance ~ Marker + Dataset + Rank + CellID,
                  cell_table, sum)
  data.frame(
    Marker = unique(sample_table$Marker), Dataset = unique(sample_table$Dataset),
    Rank = unique(sample_table$Rank),
    MinimumSampleSum = min(ss$RelativeAbundance),
    MaximumSampleSum = max(ss$RelativeAbundance),
    MinimumCellMeanSum = min(cs$MeanRelativeAbundance),
    MaximumCellMeanSum = max(cs$MeanRelativeAbundance),
    HasOther = "Other" %in% sample_table$TaxonDisplay,
    HasUnassigned = "Unassigned" %in% sample_table$TaxonDisplay,
    Status = if (max(abs(ss$RelativeAbundance - 1)) <= SUM_TOLERANCE &&
                  max(abs(cs$MeanRelativeAbundance - 1)) <= SUM_TOLERANCE &&
                  all(c("Other", "Unassigned") %in% sample_table$TaxonDisplay)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
}

write_rank_outputs <- function(stage, sample_ra, taxonomy, id_column, ranks,
                               metadata, marker, dataset) {
  invariants <- list()
  matrices <- list()
  for (rank in ranks) {
    taxon_mat <- aggregate_taxon(sample_ra, taxonomy, id_column, rank,
                                 paste(marker, dataset))
    matrices[[rank]] <- taxon_mat
    long <- sample_long(taxon_mat, metadata, marker, dataset, rank)
    cell <- decorate_cell_summary(cell_summary(long), marker, dataset, rank)
    dataset_slug <- if (dataset == "primary_rarefied") "" else if (
      dataset == "full_pr2_primary"
    ) "_full" else paste0("_", dataset)
    stem <- paste0(tolower(marker), dataset_slug, "_", safe_slug(rank))
    write_tsv(long, file.path(stage, "results",
                              paste0(stem, "_sample_relative_abundance_long_v1.tsv")))
    write_tsv(cell, file.path(stage, "results",
                              paste0(stem, "_cell_summary_v1.tsv")))
    invariants[[rank]] <- composition_invariant(
      marker, dataset, rank, long, cell, nrow(sample_ra), nrow(taxon_mat)
    )
  }
  list(matrices = matrices, invariants = do.call(rbind, invariants))
}

write_top_outputs <- function(stage, taxon_mat, metadata, marker, dataset, rank,
                              stem) {
  top_taxa <- select_top_taxa(taxon_mat, TOP_N)
  collapsed <- collapse_top(taxon_mat, top_taxa)
  sample <- top_sample_table(collapsed, metadata, marker, dataset, rank)
  cell <- top_cell_table(sample)
  write_tsv(sample, file.path(stage, "results", paste0(stem, "_figure_sample_v1.tsv")))
  write_tsv(cell, file.path(stage, "results", paste0(stem, "_figure_cell_v1.tsv")))
  selection <- data.frame(
    Marker = marker, Dataset = dataset, Rank = rank,
    TaxonDisplay = rownames(collapsed),
    CategoryOrder = seq_len(nrow(collapsed)),
    CategoryRole = ifelse(rownames(collapsed) == "Other", "other_assigned_taxa",
                          ifelse(rownames(collapsed) == "Unassigned", "unassigned",
                                 "top_assigned_taxon")),
    GlobalMeanRelativeAbundance = rowMeans(collapsed),
    stringsAsFactors = FALSE
  )
  list(selection = selection, invariant = top_invariant(sample, cell))
}

parallel_summary <- function(parallel, counts18) {
  idx <- match(parallel$ASV_ID, rownames(counts18))
  parallel$RarefiedReads <- rowSums(counts18[idx, , drop = FALSE])
  parallel$PositiveSampleCount <- rowSums(counts18[idx, , drop = FALSE] > 0)
  parallel$analysis_keep <- tolower(as.character(parallel$analysis_keep))
  parallel$flag_pr2_metazoa <- tolower(as.character(parallel$flag_pr2_metazoa))
  groups <- split(seq_len(nrow(parallel)), interaction(
    parallel$agreement_level, parallel$agreement_status,
    parallel$analysis_keep, parallel$flag_pr2_metazoa,
    drop = TRUE, lex.order = TRUE
  ))
  rows <- lapply(groups, function(i) data.frame(
    agreement_level = parallel$agreement_level[i[[1L]]],
    agreement_status = parallel$agreement_status[i[[1L]]],
    analysis_keep = parallel$analysis_keep[i[[1L]]],
    flag_pr2_metazoa = parallel$flag_pr2_metazoa[i[[1L]]],
    ASVCount = length(i),
    RarefiedReadCount = sum(parallel$RarefiedReads[i]),
    MeanPositiveSampleCount = mean(parallel$PositiveSampleCount[i]),
    stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, rows)
  out[order(out$agreement_level, out$agreement_status, out$flag_pr2_metazoa), , drop = FALSE]
}

validate_composition_in_memory <- function(data) {
  metadata <- data$metadata
  ra16 <- sample_relative_abundance(data$counts16, "16S rarefied counts")
  ra18 <- sample_relative_abundance(data$counts18, "18S rarefied counts")
  flags <- data$metazoa_flags[rownames(data$counts18)]
  nonmet_counts <- data$counts18[!flags, , drop = FALSE]
  ra18_nonmet <- sample_relative_abundance(nonmet_counts,
                                           "18S non-metazoan sensitivity counts")
  t18_nonmet <- data$taxonomy18[match(rownames(nonmet_counts), data$taxonomy18$ASV_ID),
                                , drop = FALSE]
  specs <- list(
    list(marker = "16S", dataset = "primary_rarefied", ra = ra16,
         taxonomy = data$taxonomy16, id = "OTUID",
         ranks = c("Phylum", "Class", "Genus"), top_rank = "Phylum"),
    list(marker = "18S", dataset = "full_pr2_primary", ra = ra18,
         taxonomy = data$taxonomy18, id = "ASV_ID",
         ranks = c("Kingdom", "Division_or_Phylum", "Class", "Genus"),
         top_rank = "Class"),
    list(marker = "18S", dataset = "non_metazoan_sensitivity", ra = ra18_nonmet,
         taxonomy = t18_nonmet, id = "ASV_ID",
         ranks = c("Kingdom", "Division_or_Phylum", "Class", "Genus"),
         top_rank = "Class")
  )
  invariant_rows <- list()
  top_rows <- list()
  k <- 0L
  for (spec in specs) {
    matrices <- list()
    for (rank in spec$ranks) {
      taxon_mat <- aggregate_taxon(spec$ra, spec$taxonomy, spec$id, rank,
                                   paste(spec$marker, spec$dataset))
      matrices[[rank]] <- taxon_mat
      long <- sample_long(taxon_mat, metadata, spec$marker, spec$dataset, rank)
      cell <- decorate_cell_summary(cell_summary(long), spec$marker, spec$dataset, rank)
      k <- k + 1L
      invariant_rows[[k]] <- composition_invariant(
        spec$marker, spec$dataset, rank, long, cell, nrow(spec$ra), nrow(taxon_mat)
      )
    }
    top_taxa <- select_top_taxa(matrices[[spec$top_rank]], TOP_N)
    collapsed <- collapse_top(matrices[[spec$top_rank]], top_taxa)
    top_sample <- top_sample_table(collapsed, metadata, spec$marker,
                                   spec$dataset, spec$top_rank)
    top_cell <- top_cell_table(top_sample)
    top_rows[[length(top_rows) + 1L]] <- top_invariant(top_sample, top_cell)
  }
  invariants <- do.call(rbind, invariant_rows)
  top_invariants <- do.call(rbind, top_rows)
  require_true(nrow(invariants) == 11L && all(invariants$Status == "PASS"),
               "In-memory rank-table validation failed")
  require_true(nrow(top_invariants) == 3L && all(top_invariants$Status == "PASS"),
               "In-memory plot-table validation failed")
  list(rank_invariants = nrow(invariants), top_invariants = nrow(top_invariants),
       minimum_non_metazoan_reads = min(colSums(nonmet_counts)))
}

build_package <- function(stage, data, paths) {
  for (subdir in c("results", "qc", "provenance")) {
    require_true(dir.create(file.path(stage, subdir), recursive = TRUE, showWarnings = FALSE),
                 paste0("Cannot create staging subdirectory: ", subdir))
  }
  metadata <- data$metadata
  ra16 <- sample_relative_abundance(data$counts16, "16S rarefied counts")
  ra18 <- sample_relative_abundance(data$counts18, "18S rarefied counts")

  flags <- data$metazoa_flags[rownames(data$counts18)]
  require_true(!anyNA(flags), "PR2 Metazoa flags do not cover 18S counts")
  nonmet_counts <- data$counts18[!flags, , drop = FALSE]
  require_true(nrow(nonmet_counts) > 0L, "Metazoa filtering removed every 18S feature")
  ra18_nonmet <- sample_relative_abundance(nonmet_counts,
                                           "18S non-metazoan sensitivity counts")
  t18_nonmet <- data$taxonomy18[match(rownames(nonmet_counts), data$taxonomy18$ASV_ID),
                                , drop = FALSE]

  out16 <- write_rank_outputs(stage, ra16, data$taxonomy16, "OTUID",
                              c("Phylum", "Class", "Genus"), metadata,
                              "16S", "primary_rarefied")
  out18 <- write_rank_outputs(stage, ra18, data$taxonomy18, "ASV_ID",
                              c("Kingdom", "Division_or_Phylum", "Class", "Genus"),
                              metadata, "18S", "full_pr2_primary")
  out18_nonmet <- write_rank_outputs(
    stage, ra18_nonmet, t18_nonmet, "ASV_ID",
    c("Kingdom", "Division_or_Phylum", "Class", "Genus"), metadata,
    "18S", "non_metazoan_sensitivity"
  )

  top16 <- write_top_outputs(stage, out16$matrices$Phylum, metadata, "16S",
                             "primary_rarefied", "Phylum", "16s_phylum_top12")
  top18 <- write_top_outputs(stage, out18$matrices$Class, metadata, "18S",
                             "full_pr2_primary", "Class", "18s_full_class_top12")
  top18_nonmet <- write_top_outputs(
    stage, out18_nonmet$matrices$Class, metadata, "18S",
    "non_metazoan_sensitivity", "Class",
    "18s_non_metazoan_sensitivity_class_top12"
  )

  invariants <- rbind(out16$invariants, out18$invariants, out18_nonmet$invariants)
  require_true(all(invariants$Status == "PASS"), "A composition sum invariant failed")
  write_tsv(invariants, file.path(stage, "qc", "composition_sum_invariants_v1.tsv"))
  top_invariants <- rbind(top16$invariant, top18$invariant, top18_nonmet$invariant)
  require_true(all(top_invariants$Status == "PASS"), "A top-category sum invariant failed")
  write_tsv(top_invariants, file.path(stage, "qc", "top_category_sum_invariants_v1.tsv"))
  write_tsv(rbind(top16$selection, top18$selection, top18_nonmet$selection),
            file.path(stage, "results", "top_taxa_selection_ledger_v1.tsv"))

  source_depth <- colSums(data$counts18)
  retained_depth <- colSums(nonmet_counts)
  metazoa_depth <- source_depth - retained_depth
  filter_ledger <- data.frame(
    SampleID = colnames(data$counts18),
    SourceRarefiedReads = source_depth,
    ExactPR2MetazoaReadsRemoved = metazoa_depth,
    NonMetazoanReadsRetained = retained_depth,
    NonMetazoanRetentionFraction = retained_depth / source_depth,
    RetainedFeatureCountGlobal = nrow(nonmet_counts),
    RemovedFeatureCountGlobal = sum(flags),
    RenormalizedWithinSample = TRUE,
    DatasetLabel = "non_metazoan_sensitivity",
    InterpretationRestriction = "not_microeukaryote_only",
    stringsAsFactors = FALSE
  )
  require_true(all(filter_ledger$NonMetazoanReadsRetained > 0),
               "At least one sample has zero non-metazoan reads")
  write_tsv(filter_ledger, file.path(stage, "qc", "non_metazoan_filter_ledger_v1.tsv"))

  write_tsv(parallel_summary(data$parallel18, data$counts18),
            file.path(stage, "results",
                      "18s_pr2_silva_taxonomy_concordance_summary_v1.tsv"))

  input_index <- data.frame(
    InputKey = names(paths),
    Path = unname(unlist(paths)),
    Bytes = as.numeric(file.info(unname(unlist(paths)))$size),
    SHA256 = vapply(unname(unlist(paths)), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv(input_index, file.path(stage, "provenance", "input_source_index_v1.tsv"))

  run_config <- list(
    package = PACKAGE_ID,
    script_version = SCRIPT_VERSION,
    design_package = "00_design_v2",
    samples = EXPECTED_SAMPLES,
    rarefaction_depths = list(`16S` = 20000L, `18S` = 50000L),
    composition_rule = "sample-wise relative abundance first; arithmetic cell mean second",
    cell_definition = "SiteCodeOriginal x SampleType x Campaign",
    top_n_assigned_taxa = TOP_N,
    explicit_categories = c("Other", "Unassigned"),
    taxonomy_18s_primary = "PR2 v5.1.1",
    non_metazoan_rule = "remove exact parallel-v2 flag_pr2_metazoa features; renormalize each sample",
    interpretation_restriction = "non-metazoan sensitivity; never microeukaryote-only",
    plotting_performed = FALSE,
    r_version = R.version.string,
    package_versions = list(
      digest = as.character(utils::packageVersion("digest")),
      jsonlite = as.character(utils::packageVersion("jsonlite"))
    )
  )
  write_json(run_config, file.path(stage, "provenance",
                                   "run_config_taxonomic_composition_v1.json"))
  summary <- list(
    package = PACKAGE_ID,
    status = "complete",
    samples = EXPECTED_SAMPLES,
    design_cells = length(unique(make_cell_id(metadata))),
    features_16s = nrow(data$counts16),
    features_18s_full = nrow(data$counts18),
    features_18s_exact_pr2_metazoa_removed = sum(flags),
    features_18s_non_metazoan_retained = nrow(nonmet_counts),
    composition_invariants_passed = nrow(invariants),
    top_category_invariants_passed = nrow(top_invariants),
    downstream_figures_created = FALSE
  )
  write_json(summary, file.path(stage, "provenance",
                                "package_summary_taxonomic_composition_v1.json"))

  readme <- c(
    "# Taxonomic composition package v1",
    "",
    "This package contains table-only composition outputs for the protected shared 75-sample design.",
    "16S uses 20,000 reads per sample and 18S uses 50,000 reads per sample.",
    "For every dataset, relative abundance was calculated within each sample before arithmetic means, standard deviations, and n were calculated for Site x SampleType x Campaign cells.",
    "",
    "The 18S primary tables use the full PR2 v5.1.1 analysis set. The sensitivity tables remove features marked by the exact flag_pr2_metazoa field in the protected PR2/SILVA parallel-v2 taxonomy and then renormalize each sample. They are labeled non-metazoan sensitivity and must not be described as microeukaryote-only.",
    "",
    "Plot-facing top-12 tables always contain explicit Other and Unassigned categories. CategoryOrder fixes the stack/legend order. GeologicalHabitat, SampleType, and Campaign labels come from design package 00_design_v2.",
    "",
    "No figures, differential-abundance tests, environmental models, or marker-concordance tests are generated here."
  )
  writeLines(readme, con = file.path(stage, "provenance",
                                     "README_taxonomic_composition_v1.md"),
             sep = "\n", useBytes = TRUE)
  write_manifest(stage)
}

expected_key_files <- c(
  "results/16s_phylum_sample_relative_abundance_long_v1.tsv",
  "results/16s_phylum_cell_summary_v1.tsv",
  "results/16s_phylum_top12_figure_sample_v1.tsv",
  "results/16s_phylum_top12_figure_cell_v1.tsv",
  "results/18s_full_class_sample_relative_abundance_long_v1.tsv",
  "results/18s_full_class_cell_summary_v1.tsv",
  "results/18s_full_class_top12_figure_sample_v1.tsv",
  "results/18s_full_class_top12_figure_cell_v1.tsv",
  "results/18s_non_metazoan_sensitivity_class_sample_relative_abundance_long_v1.tsv",
  "results/18s_non_metazoan_sensitivity_class_cell_summary_v1.tsv",
  "results/18s_non_metazoan_sensitivity_class_top12_figure_sample_v1.tsv",
  "results/18s_non_metazoan_sensitivity_class_top12_figure_cell_v1.tsv",
  "qc/composition_sum_invariants_v1.tsv",
  "qc/top_category_sum_invariants_v1.tsv",
  "qc/non_metazoan_filter_ledger_v1.tsv",
  "provenance/README_taxonomic_composition_v1.md",
  "provenance/run_config_taxonomic_composition_v1.json",
  "provenance/package_summary_taxonomic_composition_v1.json",
  "provenance/sha256_manifest_taxonomic_composition_v1.tsv"
)

check_columns_file <- function(path, columns) {
  x <- read_tsv(path)
  require_columns(x, columns, path)
  x
}

check_package <- function(target) {
  require_true(dir.exists(target), paste0("Missing output package: ", target))
  missing <- expected_key_files[!file.exists(file.path(target, expected_key_files))]
  require_true(length(missing) == 0L,
               paste0("Output package is missing key files: ", paste(missing, collapse = ", ")))
  manifest_path <- file.path(target, "provenance",
                             "sha256_manifest_taxonomic_composition_v1.tsv")
  manifest <- check_columns_file(manifest_path, c("package", "path", "bytes", "sha256"))
  require_true(all(manifest$package == PACKAGE_ID), "Manifest package ID mismatch")
  actual_files <- list.files(target, recursive = TRUE, full.names = FALSE)
  actual_files <- gsub("\\\\", "/", actual_files)
  actual_files <- setdiff(actual_files,
                          "provenance/sha256_manifest_taxonomic_composition_v1.tsv")
  require_true(setequal(actual_files, manifest$path),
               "Manifest file set does not match package file set")
  full <- file.path(target, manifest$path)
  require_true(all(file.exists(full)), "Manifest references a missing file")
  require_true(all(as.numeric(file.info(full)$size) == as.numeric(manifest$bytes)),
               "Manifest byte-size mismatch")
  hashes <- vapply(full, sha256_file, character(1))
  require_true(identical(unname(hashes), as.character(manifest$sha256)),
               "Manifest SHA-256 mismatch")

  sample_schema <- c("Marker", "Dataset", "Rank", "SampleID", "Site",
                     "SiteCodeOriginal", "GeologicalHabitat", "SampleType",
                     "Campaign", "CellID", "TaxonDisplay", "CategoryOrder",
                     "RelativeAbundance")
  cell_schema <- c("Marker", "Dataset", "Rank", "CellID", "Site",
                   "SiteCodeOriginal", "GeologicalHabitat", "SampleType",
                   "Campaign", "TaxonDisplay", "CategoryOrder",
                   "MeanRelativeAbundance", "SDRelativeAbundance", "N")
  for (stem in c("16s_phylum_top12", "18s_full_class_top12",
                 "18s_non_metazoan_sensitivity_class_top12")) {
    s <- check_columns_file(file.path(target, "results",
                                      paste0(stem, "_figure_sample_v1.tsv")), sample_schema)
    c <- check_columns_file(file.path(target, "results",
                                      paste0(stem, "_figure_cell_v1.tsv")), cell_schema)
    sums <- aggregate(RelativeAbundance ~ SampleID, s, sum)$RelativeAbundance
    cell_sums <- aggregate(MeanRelativeAbundance ~ CellID, c, sum)$MeanRelativeAbundance
    require_true(max(abs(sums - 1)) <= SUM_TOLERANCE &&
                   max(abs(cell_sums - 1)) <= SUM_TOLERANCE,
                 paste0("Top-table sum check failed: ", stem))
    require_true(all(c("Other", "Unassigned") %in% s$TaxonDisplay),
                 paste0("Explicit category check failed: ", stem))
  }
  inv <- read_tsv(file.path(target, "qc", "composition_sum_invariants_v1.tsv"))
  top_inv <- read_tsv(file.path(target, "qc", "top_category_sum_invariants_v1.tsv"))
  require_true(all(inv$Status == "PASS") && all(top_inv$Status == "PASS"),
               "Stored invariant ledger contains a failure")
  invisible(TRUE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  for (pkg in c("digest", "jsonlite")) {
    require_true(requireNamespace(pkg, quietly = TRUE),
                 paste0("Required R package is not installed: ", pkg))
  }
  paths <- input_paths(opts$root)
  data <- validate_sources(paths)
  target <- file.path(opts$root, PACKAGE_ID)
  mode <- if (opts$run) "run" else if (opts$dry_run) "dry-run" else "check"
  cat(paste0("script_version=", SCRIPT_VERSION, "\n"))
  cat(paste0("mode=", mode, "\n"))
  cat(paste0("analysis_root=", opts$root, "\n"))
  cat(paste0("design_package=00_design_v2\n"))
  cat(paste0("samples=", nrow(data$metadata), "\n"))
  cat(paste0("design_cells=", length(unique(make_cell_id(data$metadata))), "\n"))
  cat(paste0("features_16s=", nrow(data$counts16), "\n"))
  cat(paste0("features_18s=", nrow(data$counts18), "\n"))
  cat(paste0("features_flag_pr2_metazoa=", sum(data$metazoa_flags), "\n"))
  cat("composition_rule=sample-wise_RA_then_cell_arithmetic_mean_sd_n\n")
  cat("non_metazoan_label_restriction=never_microeukaryote_only\n")

  if (opts$dry_run) {
    smoke <- validate_composition_in_memory(data)
    cat(paste0("in_memory_rank_invariants=", smoke$rank_invariants, "\n"))
    cat(paste0("in_memory_top_invariants=", smoke$top_invariants, "\n"))
    cat(paste0("minimum_non_metazoan_reads=", smoke$minimum_non_metazoan_reads, "\n"))
    cat(paste0("target_exists=", dir.exists(target), "\n"))
    cat("status=dry_run_validated_no_outputs_written\n")
    return(invisible(NULL))
  }
  if (opts$check) {
    check_package(target)
    cat("status=check_passed\n")
    return(invisible(NULL))
  }

  require_true(!file.exists(target), paste0("Refusing to overwrite existing target: ", target))
  stage <- paste0(target, ".staging_", Sys.getpid())
  require_true(!file.exists(stage), paste0("Staging path already exists: ", stage))
  require_true(dir.create(stage, recursive = FALSE, showWarnings = FALSE),
               paste0("Cannot create sibling staging directory: ", stage))
  published <- FALSE
  on.exit({
    if (!published && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  build_package(stage, data, paths)
  check_package(stage)
  require_true(file.rename(stage, target),
               paste0("Atomic staging rename failed: ", stage, " -> ", target))
  published <- TRUE
  check_package(target)
  cat(paste0("output_package=", normalizePath(target, winslash = "/"), "\n"))
  cat("status=run_complete\n")
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1L, runLast = FALSE)
})
