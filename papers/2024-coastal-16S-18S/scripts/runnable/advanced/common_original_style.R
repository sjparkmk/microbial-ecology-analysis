suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
  invisible(TRUE)
}

read_tsv_strict <- function(path) {
  assert_true(file.exists(path), paste0("Missing input: ", path))
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv_atomic <- function(x, path, row.names = FALSE, col.names = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(tmp)) stop("Refusing overwrite: ", path, call. = FALSE)
  write.table(x, tmp, sep = "\t", row.names = row.names, col.names = col.names,
              quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic promote failed: ", path, call. = FALSE)
  invisible(path)
}

load_count_matrix <- function(path) {
  x <- read_tsv_strict(path)
  assert_true(ncol(x) >= 2L, paste0("Count table has fewer than two columns: ", path))
  ids <- as.character(x[[1]])
  assert_true(!anyNA(ids) && !anyDuplicated(ids), paste0("Invalid feature IDs: ", path))
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  assert_true(all(is.finite(mat)) && all(mat >= 0) && all(mat == floor(mat)),
              paste0("Counts must be finite nonnegative integers: ", path))
  list(id_col = names(x)[1], matrix = mat)
}

load_design <- function(path) {
  x <- read_tsv_strict(path)
  need <- c("SampleID", "Site", "GeologicalHabitat", "SampleType", "Season",
            "Campaign", "FieldReplicate", "LocationID", "PairID_Substrate",
            "PairID_Season", "SiteSeasonBlock")
  miss <- setdiff(need, names(x))
  assert_true(length(miss) == 0L, paste0("Design columns missing: ", paste(miss, collapse = ", ")))
  assert_true(nrow(x) == 75L && !anyDuplicated(x$SampleID), "Design must contain 75 unique samples")
  x$Site <- factor(x$Site, levels = c("Dokdo", "Goseong", "Pohang", "Pyoseon",
                                             "Sehwa", "Seongsan", "Ulleungdo"))
  x$GeologicalHabitat <- factor(x$GeologicalHabitat, levels = c("Ocean", "Tidal", "Island"))
  x$SampleType <- factor(x$SampleType, levels = c("Sea-sand", "Water"))
  x$Season <- factor(x$Season, levels = c("Spring", "Autumn"))
  x
}

align_counts_design <- function(mat, design) {
  assert_true(setequal(colnames(mat), design$SampleID), "Count/design SampleID sets differ")
  mat <- mat[, design$SampleID, drop = FALSE]
  list(matrix = mat, design = design)
}

taxonomy_id_column <- function(x) {
  hits <- intersect(c("OTUID", "ASV_ID", "#OTUID", "FeatureID"), names(x))
  if (length(hits)) hits[[1]] else names(x)[[1]]
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
  assert_true(!anyNA(idx), paste0("Taxonomy linkage failure for rank ", rank))
  lab <- normalize_tax_label(taxonomy[[rank]][idx])
  out <- rowsum(mat, group = lab, reorder = FALSE)
  out[rowSums(out) > 0, , drop = FALSE]
}

relative_abundance <- function(mat) {
  den <- colSums(mat)
  assert_true(all(den > 0), "Zero-depth sample in relative-abundance input")
  sweep(mat, 2, den, "/")
}

clr_transform <- function(mat, pseudocount = 0.5) {
  x <- log(t(mat) + pseudocount)
  sweep(x, 1, rowMeans(x), "-")
}

cell_id <- function(design) paste(design$Site, design$SampleType, design$Season, sep = "|")

sample_long_relative <- function(mat, design, marker, rank, dataset) {
  rel <- relative_abundance(mat)
  x <- as.data.frame(as.table(rel), stringsAsFactors = FALSE)
  names(x) <- c("Taxon", "SampleID", "RelativeAbundance")
  x$RelativeAbundancePct <- 100 * x$RelativeAbundance
  x <- left_join(x, design, by = "SampleID")
  x$Marker <- marker
  x$Rank <- rank
  x$Dataset <- dataset
  x$CellID <- cell_id(x)
  x
}

cell_summary_relative <- function(long) {
  long %>%
    group_by(Marker, Dataset, Rank, Site, GeologicalHabitat, SampleType, Season, CellID, Taxon) %>%
    summarise(MeanRelativeAbundancePct = mean(RelativeAbundancePct),
              SDRelativeAbundancePct = sd(RelativeAbundancePct),
              N = n(), .groups = "drop")
}

select_top_taxa <- function(long, n = 12L) {
  totals <- long %>% group_by(Taxon) %>% summarise(m = mean(RelativeAbundancePct), .groups = "drop") %>%
    arrange(desc(m), Taxon)
  keep <- head(totals$Taxon[totals$Taxon != "Unassigned"], n)
  long$TaxonDisplay <- ifelse(long$Taxon %in% keep, long$Taxon,
                              ifelse(long$Taxon == "Unassigned", "Unassigned", "Other assigned"))
  list(data = long, selected = keep)
}

habitat_palette <- c(Ocean = "#0072B2", Tidal = "#E69F00", Island = "#009E73")
sampletype_palette <- c(`Sea-sand` = "#B97845", Water = "#3B82B8")
season_palette <- c(Spring = "#67A9CF", Autumn = "#EF8A62")

taxon_palette <- function(labels) {
  labels <- unique(as.character(labels))
  fixed <- c(`Other assigned` = "#BDBDBD", Unassigned = "#666666")
  main <- setdiff(labels, names(fixed))
  hues <- grDevices::hcl(h = seq(15, 375, length.out = length(main) + 1L)[seq_along(main)],
                         c = 55, l = 70)
  c(stats::setNames(hues, main), fixed[intersect(names(fixed), labels)])
}

theme_project_original_style <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
          strip.text = element_text(color = "black", face = "bold"),
          legend.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          plot.subtitle = element_text(color = "#333333"),
          plot.caption = element_text(color = "#555555", hjust = 0))
}

save_plot_bundle <- function(plot, base_path, width_mm, height_mm, dpi = 600) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  paths <- paste0(base_path, c(".pdf", ".svg", ".png"))
  if (any(file.exists(paths))) stop("Refusing figure overwrite: ", base_path, call. = FALSE)
  ggsave(paths[[1]], plot, width = width_mm, height = height_mm, units = "mm",
         device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(paths[[2]], plot, width = width_mm, height = height_mm, units = "mm",
         device = "svg", bg = "white", limitsize = FALSE)
  ggsave(paths[[3]], plot, width = width_mm, height = height_mm, units = "mm",
         dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(paths)
}

safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

classification_status <- function(valid, informative, manuscript = FALSE, supplement = FALSE) {
  if (!valid) return("failed_qc")
  if (!informative) return("not_informative")
  if (manuscript) return("main_candidate")
  if (supplement) return("supplement_candidate")
  "descriptive_only"
}
