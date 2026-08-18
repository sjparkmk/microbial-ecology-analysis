#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath("/your/data/", winslash = "/", mustWork = TRUE) # insert your path
pkg <- file.path(root, "revision/re_analysis/54_revision_niche_breadth_common75_v1")
dirs <- file.path(pkg, c("source_data", "results", "provenance", "qc"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

inputs <- list(
  metadata = file.path(root, "revision/re_analysis/00_design_v4/master_sample_metadata_v4.tsv"),
  membership = file.path(root, "revision/re_analysis/00_design_v4/analysis_set_membership_v4.tsv"),
  count16 = file.path(root, "revision/re_analysis/04_input_16s_rarefied_20000_v1/results/16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  tax16 = file.path(root, "revision/re_analysis/04_input_16s_rarefied_20000_v1/results/16s_taxonomy_rarefied_20000_analysis_v1.tsv"),
  count18 = file.path(root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
  tax18 = file.path(root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_taxonomy_pr2_primary_rarefied_50000_analysis_v1.tsv")
)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
assert <- function(ok, ...) if (!isTRUE(ok)) stopf(...)
read_tsv <- function(path) {
  assert(file.exists(path), "Missing input: %s", path)
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", comment.char = "", na.strings = character())
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "")
}
is_assigned <- function(x) {
  z <- trimws(as.character(x))
  nzchar(z) & !(tolower(z) %in% c("na", "n/a", "unassigned", "unknown", "unclassified"))
}
levins <- function(x) {
  n <- ncol(x)
  s <- rowSums(x)
  ss <- rowSums(x * x)
  b <- rep(NA_real_, nrow(x))
  ok <- s > 0 & ss > 0 & n > 1
  b[ok] <- s[ok]^2 / ss[ok]
  ba <- (b - 1) / (n - 1)
  list(B = b, BA = ba)
}
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
qfun <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(quantile(x, p, type = 8, names = FALSE))
}

metadata <- read_tsv(inputs$metadata)
membership <- read_tsv(inputs$membership)
assert(nrow(metadata) == 75L, "Expected 75 metadata rows, found %d", nrow(metadata))
assert(!anyDuplicated(metadata$SampleID), "Duplicate SampleID in metadata")
assert(all(metadata$Included16S == "TRUE") && all(metadata$Included18S == "TRUE"),
       "Metadata is not the common-75 marker set")
metadata$CellID <- paste(metadata$Site, metadata$Season, metadata$SampleType, sep = "|")
cell_n <- table(metadata$CellID)
assert(length(cell_n) == 25L && all(cell_n == 3L),
       "Expected 25 observed cells with three samples each")

a60_ids <- membership$SampleID[membership$SetID == "A60_COMPLETE_FACTORIAL" &
                                 membership$Included == "TRUE"]
assert(length(a60_ids) == 60L && !anyDuplicated(a60_ids), "Expected 60 unique A60 samples")
assert(all(a60_ids %in% metadata$SampleID), "A60 includes samples outside common 75")

set.seed(76109)
n_boot <- 2000L
site_levels <- unique(metadata$Site)

run_marker <- function(marker, rank, count_path, tax_path, feature_id,
                       depth_expected, higher_col) {
  counts_df <- read_tsv(count_path)
  taxonomy <- read_tsv(tax_path)
  assert(names(counts_df)[1] %in% c("#OTUID", "OTUID", "ASV_ID"),
         "%s count ID column unexpected", marker)
  names(counts_df)[1] <- "FeatureID"
  names(taxonomy)[names(taxonomy) == feature_id] <- "FeatureID"
  assert(all(c("FeatureID", rank, higher_col) %in% names(taxonomy)),
         "%s taxonomy columns missing", marker)
  assert(!anyDuplicated(counts_df$FeatureID), "%s duplicate count feature IDs", marker)
  assert(!anyDuplicated(taxonomy$FeatureID), "%s duplicate taxonomy feature IDs", marker)
  assert(setequal(names(counts_df)[-1], metadata$SampleID),
         "%s sample columns do not equal common-75 metadata", marker)
  sample_ids <- metadata$SampleID
  count_mat <- as.matrix(counts_df[, sample_ids, drop = FALSE])
  storage.mode(count_mat) <- "numeric"
  rownames(count_mat) <- counts_df$FeatureID
  assert(all(colSums(count_mat) == depth_expected),
         "%s does not have expected fixed depth %d", marker, depth_expected)
  ti <- match(rownames(count_mat), taxonomy$FeatureID)
  assert(!anyNA(ti), "%s count features missing taxonomy", marker)
  taxonomy <- taxonomy[ti, , drop = FALSE]
  keep_asv <- is_assigned(taxonomy[[rank]])
  assert(any(keep_asv), "%s has no assigned %s features", marker, rank)
  assigned_counts <- count_mat[keep_asv, , drop = FALSE]
  assigned_tax <- trimws(taxonomy[[rank]][keep_asv])
  agg_counts <- rowsum(assigned_counts, group = assigned_tax, reorder = TRUE)
  agg_ra <- sweep(agg_counts, 2L, colSums(count_mat), "/")

  prevalence_n <- rowSums(agg_counts > 0)
  prevalence <- prevalence_n / ncol(agg_counts)
  sites_detected <- vapply(seq_len(nrow(agg_counts)), function(i) {
    sum(vapply(site_levels, function(s) {
      any(agg_counts[i, metadata$Site == s, drop = TRUE] > 0)
    }, logical(1L)))
  }, integer(1L))
  primary <- prevalence >= 0.10 & sites_detected >= 3L
  assert(any(primary), "%s has no taxa at prevalence >=10%% and >=3 Sites", marker)
  agg_counts <- agg_counts[primary, , drop = FALSE]
  agg_ra <- agg_ra[primary, , drop = FALSE]
  prevalence_n <- prevalence_n[primary]
  prevalence <- prevalence[primary]
  sites_detected <- sites_detected[primary]
  taxa <- rownames(agg_ra)

  split_cells <- split(seq_len(nrow(metadata)), metadata$CellID)
  cell_ra <- vapply(split_cells, function(i) rowMeans(agg_ra[, i, drop = FALSE]),
                    numeric(nrow(agg_ra)))
  rownames(cell_ra) <- taxa
  cell_meta <- metadata[match(colnames(cell_ra), metadata$CellID),
                        c("CellID", "Site", "Season", "SampleType"), drop = FALSE]
  assert(!anyNA(cell_meta$Site), "%s cell metadata match failed", marker)
  full <- levins(cell_ra)

  a60_meta <- metadata[match(a60_ids, metadata$SampleID), , drop = FALSE]
  a60_cell_ids <- unique(a60_meta$CellID)
  assert(length(a60_cell_ids) == 20L, "%s A60 must have 20 cells", marker)
  a60_cell <- cell_ra[, a60_cell_ids, drop = FALSE]
  assert(all(table(a60_meta$CellID) == 3L), "%s A60 cells are not triplicates", marker)
  a60 <- levins(a60_cell)
  a60_samples_detected <- rowSums(agg_counts[, a60_ids, drop = FALSE] > 0)
  a60_cells_detected <- rowSums(a60_cell > 0)
  a60_evaluable <- a60_samples_detected >= 6L & a60_cells_detected >= 3L

  loso_rows <- list()
  loso_mat <- matrix(NA_real_, nrow = nrow(cell_ra), ncol = length(site_levels),
                     dimnames = list(taxa, site_levels))
  for (s in site_levels) {
    idx <- which(cell_meta$Site != s)
    z <- levins(cell_ra[, idx, drop = FALSE])$BA
    loso_mat[, s] <- z
    loso_rows[[s]] <- data.frame(
      Marker = marker, Rank = rank, Taxon = taxa, OmittedSite = s,
      RemainingCells = length(idx), StandardizedLevinsBA = z,
      DeltaFromFull = z - full$BA, stringsAsFactors = FALSE
    )
  }

  # Preserve the same 25 sampled context states in every bootstrap. Within each
  # Site, resample its three fixed spatial positions with replacement and carry
  # each selected position jointly across all available Sample type x Season
  # records. This preserves paired/repeated covariance without treating Sites as
  # a superpopulation; Site omission is evaluated separately below.
  boot <- matrix(NA_real_, nrow = nrow(cell_ra), ncol = n_boot)
  site_cells <- split(cell_meta$CellID, cell_meta$Site)
  site_locations <- split(metadata$LocationID, metadata$Site)
  site_locations <- lapply(site_locations, unique)
  assert(all(vapply(site_locations, length, integer(1L)) == 3L),
         "%s expected three fixed positions per Site", marker)
  for (b in seq_len(n_boot)) {
    boot_cell <- matrix(NA_real_, nrow = nrow(agg_ra), ncol = length(split_cells),
                        dimnames = list(taxa, names(split_cells)))
    for (s in site_levels) {
      draw <- sample(site_locations[[s]], 3L, replace = TRUE)
      for (cid in site_cells[[s]]) {
        i <- split_cells[[cid]]
        chosen <- i[match(draw, metadata$LocationID[i])]
        assert(!anyNA(chosen), "%s position-to-cell bootstrap match failed", marker)
        boot_cell[, cid] <- rowMeans(agg_ra[, chosen, drop = FALSE])
      }
    }
    boot[, b] <- levins(boot_cell)$BA
  }
  valid_boot <- rowSums(is.finite(boot))
  boot_low <- apply(boot, 1L, qfun, p = 0.025)
  boot_med <- apply(boot, 1L, qfun, p = 0.5)
  boot_high <- apply(boot, 1L, qfun, p = 0.975)

  lineage <- aggregate(taxonomy[keep_asv, higher_col, drop = FALSE],
                       by = list(Taxon = assigned_tax),
                       FUN = function(x) paste(sort(unique(x[is_assigned(x)])), collapse = "; "))
  names(lineage)[2] <- "HigherTaxonomy"
  lineage <- lineage[match(taxa, lineage$Taxon), , drop = FALSE]
  n_asv <- table(assigned_tax)[taxa]

  out <- data.frame(
    Marker = marker, Rank = rank, Taxon = taxa,
    HigherTaxonomy = lineage$HigherTaxonomy,
    AggregatedASVs = as.integer(n_asv),
    SamplesDetected = as.integer(prevalence_n),
    SitesDetected = as.integer(sites_detected),
    SamplePrevalence = prevalence,
    PassPrevalence10Pct = TRUE,
    PassPrevalence20Pct = prevalence >= 0.20 & sites_detected >= 3L,
    GrandMeanRelativeAbundancePct = rowMeans(agg_ra) * 100,
    ObservedCells = ncol(cell_ra),
    LevinsB = full$B,
    StandardizedLevinsBA = full$BA,
    SiteStratifiedPositionBootstrapValid = valid_boot,
    SiteStratifiedPositionBootstrapMedianBA = boot_med,
    SiteStratifiedPositionBootstrapLower95BA = boot_low,
    SiteStratifiedPositionBootstrapUpper95BA = boot_high,
    LOSOMinBA = apply(loso_mat, 1L, min, na.rm = TRUE),
    LOSOMaxBA = apply(loso_mat, 1L, max, na.rm = TRUE),
    LOSOMaxAbsoluteDelta = apply(abs(sweep(loso_mat, 1L, full$BA, "-")), 1L, max,
                                 na.rm = TRUE),
    A60ObservedCells = ncol(a60_cell),
    A60SamplesDetected = a60_samples_detected,
    A60CellsDetected = a60_cells_detected,
    A60Evaluable = a60_evaluable,
    A60StandardizedLevinsBA = ifelse(a60_evaluable, a60$BA, NA_real_),
    A60DeltaFromA75 = ifelse(a60_evaluable, a60$BA - full$BA, NA_real_),
    stringsAsFactors = FALSE
  )
  out$BreadthRankDescending <- rank(-out$StandardizedLevinsBA, ties.method = "min")
  full_q <- quantile(full$BA, c(0.1, 0.5, 0.9), type = 8, na.rm = TRUE)
  a60_q <- quantile(a60$BA[a60_evaluable], c(0.1, 0.5, 0.9), type = 8, na.rm = TRUE)
  full_tail <- ifelse(full$BA >= full_q[[3]], "Upper_decile",
                      ifelse(full$BA <= full_q[[1]], "Lower_decile", "Middle_80pct"))
  a60_tail <- ifelse(a60$BA >= a60_q[[3]], "Upper_decile",
                     ifelse(a60$BA <= a60_q[[1]], "Lower_decile", "Middle_80pct"))
  a60_tail[!a60_evaluable] <- "Not_evaluable"
  loso_same <- integer(length(taxa))
  opposite_half <- logical(length(taxa))
  for (s in site_levels) {
    z <- loso_mat[, s]
    zq <- quantile(z, c(0.1, 0.5, 0.9), type = 8, na.rm = TRUE)
    ztail <- ifelse(z >= zq[[3]], "Upper_decile",
                    ifelse(z <= zq[[1]], "Lower_decile", "Middle_80pct"))
    loso_same <- loso_same + (ztail == full_tail & full_tail != "Middle_80pct")
    opposite_half <- opposite_half |
      (full_tail == "Upper_decile" & z < zq[[2]]) |
      (full_tail == "Lower_decile" & z > zq[[2]])
  }
  opposite_half <- opposite_half |
    (a60_evaluable & full_tail == "Upper_decile" & a60$BA < a60_q[[2]]) |
    (a60_evaluable & full_tail == "Lower_decile" & a60$BA > a60_q[[2]])
  out$FullDataBreadthDecile <- full_tail[match(out$Taxon, taxa)]
  out$A60SameBreadthDecile <- (a60_tail == full_tail)[match(out$Taxon, taxa)] &
    out$FullDataBreadthDecile != "Middle_80pct"
  out$LOSOSameBreadthDecileCount <- loso_same[match(out$Taxon, taxa)]
  out$OppositeHalfReversal <- opposite_half[match(out$Taxon, taxa)]
  out$EligibleAsDescriptiveTailExample <- out$FullDataBreadthDecile != "Middle_80pct" &
    out$A60SameBreadthDecile & out$LOSOSameBreadthDecileCount >= 6L &
    !out$OppositeHalfReversal
  out <- out[order(-out$StandardizedLevinsBA, out$Taxon), ]

  upper_examples <- head(out[out$FullDataBreadthDecile == "Upper_decile" &
                               out$EligibleAsDescriptiveTailExample &
                               out$PassPrevalence20Pct, ], 5L)
  lower_pool <- out[out$FullDataBreadthDecile == "Lower_decile" &
                      out$EligibleAsDescriptiveTailExample &
                      out$PassPrevalence20Pct, ]
  lower_examples <- tail(lower_pool, 5L)
  examples <- rbind(upper_examples, lower_examples)
  examples$ExampleScope <- "descriptive sampled-context breadth; not an ecological niche category"

  cell_long <- do.call(rbind, lapply(seq_len(nrow(cell_ra)), function(i) {
    data.frame(
      Marker = marker, Rank = rank, Taxon = taxa[i],
      CellID = cell_meta$CellID, Site = cell_meta$Site,
      Season = cell_meta$Season, SampleType = cell_meta$SampleType,
      SpatialReplicates = 3L,
      CellMeanSampleWiseRelativeAbundance = as.numeric(cell_ra[i, ]),
      stringsAsFactors = FALSE
    )
  }))
  loso <- do.call(rbind, loso_rows)

  loso_cor <- vapply(site_levels, function(s) safe_cor(full$BA, loso_mat[, s]), numeric(1L))
  stability <- data.frame(
    Marker = marker, Rank = rank,
    PrimaryTaxa = nrow(out),
    Prevalence20PctTaxa = sum(out$PassPrevalence20Pct),
    A60EvaluableTaxa = sum(a60_evaluable),
    A75A60SpearmanAll = safe_cor(full$BA[a60_evaluable], a60$BA[a60_evaluable]),
    A75A60SpearmanPrevalence20Pct = safe_cor(full$BA[out$PassPrevalence20Pct & out$A60Evaluable],
                                              a60$BA[out$PassPrevalence20Pct & out$A60Evaluable]),
    MedianAbsoluteA60A75Delta = median(abs(a60$BA[a60_evaluable] - full$BA[a60_evaluable]),
                                       na.rm = TRUE),
    MinimumLOSOSpearmanAll = min(loso_cor, na.rm = TRUE),
    MedianLOSOSpearmanAll = median(loso_cor, na.rm = TRUE),
    MinimumLOSOSpearmanPrevalence20Pct = min(vapply(site_levels, function(s) {
      k <- out$PassPrevalence20Pct
      safe_cor(full$BA[k], loso_mat[k, s])
    }, numeric(1L)), na.rm = TRUE),
    MedianBootstrapCIWidth = median(out$SiteStratifiedPositionBootstrapUpper95BA -
                                      out$SiteStratifiedPositionBootstrapLower95BA, na.rm = TRUE),
    TaxaWithAllValidBootstraps = sum(out$SiteStratifiedPositionBootstrapValid == n_boot),
    LOSOSpearmanAtLeast080 = sum(loso_cor >= 0.80),
    SupportGate = safe_cor(full$BA[a60_evaluable], a60$BA[a60_evaluable]) >= 0.80 &&
      median(abs(a60$BA[a60_evaluable] - full$BA[a60_evaluable]), na.rm = TRUE) <= 0.10 &&
      sum(loso_cor >= 0.80) >= 6L && min(loso_cor, na.rm = TRUE) >= 0.70,
    stringsAsFactors = FALSE
  )
  omission_cor <- data.frame(
    Marker = marker, Rank = rank, OmittedSite = site_levels,
    SpearmanVsFull = loso_cor, stringsAsFactors = FALSE
  )

  list(summary = out, examples = examples, cell = cell_long, loso = loso,
       stability = stability, omission_cor = omission_cor,
       n_features = nrow(count_mat), n_assigned_features = sum(keep_asv),
       n_assigned_taxa = nrow(agg_counts), depth = depth_expected)
}

r16 <- run_marker("16S", "Genus", inputs$count16, inputs$tax16, "OTUID", 20000L,
                  "Phylum")
r18 <- run_marker("18S", "Class", inputs$count18, inputs$tax18, "ASV_ID", 50000L,
                  "Division_or_Phylum")

summary_all <- rbind(r16$summary, r18$summary)
examples_all <- rbind(r16$examples, r18$examples)
cell_all <- rbind(r16$cell, r18$cell)
loso_all <- rbind(r16$loso, r18$loso)
stability_all <- rbind(r16$stability, r18$stability)
omission_all <- rbind(r16$omission_cor, r18$omission_cor)

write_tsv(summary_all, file.path(pkg, "results/distributional_breadth_taxon_summary_v1.tsv"))
write_tsv(examples_all, file.path(pkg, "results/distributional_breadth_descriptive_examples_v1.tsv"))
write_tsv(cell_all, file.path(pkg, "source_data/distributional_breadth_cell_means_v1.tsv"))
write_tsv(loso_all, file.path(pkg, "results/distributional_breadth_leave_one_site_out_v1.tsv"))
write_tsv(stability_all, file.path(pkg, "results/distributional_breadth_marker_stability_v1.tsv"))
write_tsv(omission_all, file.path(pkg, "results/distributional_breadth_loso_correlations_v1.tsv"))

run_summary <- data.frame(
  Item = c("Common samples per marker", "Observed cells", "Replicates per observed cell",
           "A60 samples", "A60 cells", "Primary prevalence threshold",
           "Sensitivity prevalence threshold", "Site-stratified position bootstrap replicates", "Bootstrap seed",
           "16S fixed depth", "18S fixed depth", "16S source ASVs",
           "18S source ASVs", "16S primary taxa", "18S primary taxa"),
  Value = c(75, 25, 3, 60, 20, "10% at count > 0", "20% at count > 0", n_boot,
            76109, r16$depth, r18$depth, r16$n_features, r18$n_features,
            nrow(r16$summary), nrow(r18$summary)),
  stringsAsFactors = FALSE
)
write_tsv(run_summary, file.path(pkg, "provenance/run_summary_v1.tsv"))

qa <- data.frame(
  Check = c("metadata_rows_75", "observed_cells_25", "cell_triplicates",
            "a60_samples_60", "a60_cells_20", "fixed_depths_exact",
            "ba_range", "no_forced_labels", "both_marker_support_gate"),
  Pass = c(nrow(metadata) == 75L, length(cell_n) == 25L, all(cell_n == 3L),
           length(a60_ids) == 60L, length(unique(metadata$CellID[metadata$SampleID %in% a60_ids])) == 20L,
           TRUE, all(summary_all$StandardizedLevinsBA >= -1e-12 &
                       summary_all$StandardizedLevinsBA <= 1 + 1e-12),
           !any(grepl("generalist|specialist", names(summary_all), ignore.case = TRUE)),
           all(stability_all$SupportGate)),
  Detail = c("common marker set", "available Site x Season x Sample type cells",
             "three separate spatial field replicates per cell", "complete factorial sensitivity",
             "five Sites x two Seasons x two Sample types", "16S=20000; 18S=50000",
             "standardized Levins B_A constrained to [0,1]", "continuous breadth only",
             "gate: A75-A60 Spearman >=0.80; >=6/7 Site omissions >=0.80; none <0.70"),
  stringsAsFactors = FALSE
)
write_tsv(qa, file.path(pkg, "qc/distributional_breadth_validation_v1.tsv"))
assert(all(qa$Pass), "One or more QA checks failed; inspect validation table")

sink(file.path(pkg, "provenance/sessionInfo_v1.txt"))
print(sessionInfo())
sink()

message("Distributional-breadth analysis complete: ", pkg)
