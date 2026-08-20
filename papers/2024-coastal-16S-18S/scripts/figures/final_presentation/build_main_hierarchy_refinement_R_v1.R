#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(svglite)
  library(ragg)
  library(jsonlite)
  library(digest)
  library(pdftools)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 5)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script path.")
script_file <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)
package_root <- normalizePath(file.path(script_dir, "..", "generated", "hierarchy"), winslash = "/", mustWork = FALSE)
project_root <- normalizePath(file.path(script_dir, "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)

paths <- list(
  package = package_root,
  source = file.path(package_root, "source_data"),
  figures = file.path(package_root, "figures"),
  previews = file.path(package_root, "previews"),
  manifests = file.path(package_root, "manifests"),
  qc = file.path(package_root, "qc"),
  captions = file.path(package_root, "caption_facts")
)
invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE)

project_relative <- function(path) {
  p <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(project_root, "/")
  if (startsWith(tolower(p), tolower(prefix))) substring(p, nchar(prefix) + 1L) else p
}

write_package_manifest <- function() {
  manifest_path <- file.path(paths$manifests, "final_hierarchy_figures_v1_sha256_manifest.tsv")
  files <- list.files(package_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[normalizePath(files, winslash = "/", mustWork = TRUE) !=
                   normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
  manifest <- tibble(
    PackagePath = vapply(files, function(x) {
      p <- normalizePath(x, winslash = "/", mustWork = TRUE)
      substring(p, nchar(normalizePath(package_root, winslash = "/", mustWork = TRUE)) + 2L)
    }, character(1)),
    Bytes = as.numeric(file.info(files)$size),
    SHA256 = vapply(files, sha256, character(1))
  ) %>% arrange(PackagePath)
  write_tsv(manifest, manifest_path)
  invisible(manifest)
}

if ("--manifest-only" %in% args) {
  write_package_manifest()
  quit(save = "no", status = 0L)
}

accepted_inputs <- tibble::tribble(
  ~LogicalName, ~OriginalRelativePath, ~PackageName,
  "16S cPCoA scores", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S3_16S_beta_diversity_v1_cpcoa_scores.tsv", "accepted_16S_cpcoa_scores.tsv",
  "16S Bray-Curtis matrix", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S3_16S_beta_diversity_v1_bray_matrix.tsv", "accepted_16S_bray_matrix.tsv",
  "16S heatmap order", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S3_16S_beta_diversity_v1_sample_order.tsv", "accepted_16S_sample_order.tsv",
  "18S cPCoA scores", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S4_18S_beta_diversity_v1_cpcoa_scores.tsv", "accepted_18S_cpcoa_scores.tsv",
  "18S Bray-Curtis matrix", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S4_18S_beta_diversity_v1_bray_matrix.tsv", "accepted_18S_bray_matrix.tsv",
  "18S heatmap order", "revision/re_analysis/22_revision_supp_beta_figures_v1/source_data/Supplementary_Figure_S4_18S_beta_diversity_v1_sample_order.tsv", "accepted_18S_sample_order.tsv",
  "Site-synchronized alpha effects", "revision/rev_ms/rev_v1/01_hierarchical_exchangeability_v2/results/alpha_A72_site_synchronized_results_v2.tsv", "accepted_alpha_A72_site_synchronized_results_v2.tsv",
  "Alpha multiplicity sensitivity", "revision/rev_ms/rev_v1/01_hierarchical_exchangeability_v2/results/alpha_five_contrast_family_sensitivity_v2.tsv", "accepted_alpha_five_contrast_family_sensitivity_v2.tsv",
  "Hierarchical beta PERMANOVA", "revision/rev_ms/rev_v1/01_hierarchical_exchangeability_v2/results/beta_hierarchical_permanova_v2.tsv", "accepted_beta_hierarchical_permanova_v2.tsv",
  "Hierarchical beta PERMDISP", "revision/rev_ms/rev_v1/01_hierarchical_exchangeability_v2/results/beta_hierarchical_permdisp_v2.tsv", "accepted_beta_hierarchical_permdisp_v2.tsv"
) %>%
  mutate(
    OriginalPath = file.path(project_root, OriginalRelativePath),
    PackagePath = file.path(paths$source, PackageName)
  )

if (!all(file.exists(accepted_inputs$OriginalPath))) {
  stop("One or more accepted input TSVs are missing: ",
       paste(accepted_inputs$OriginalPath[!file.exists(accepted_inputs$OriginalPath)], collapse = "; "))
}

for (i in seq_len(nrow(accepted_inputs))) {
  ok <- file.copy(accepted_inputs$OriginalPath[i], accepted_inputs$PackagePath[i], overwrite = TRUE)
  if (!isTRUE(ok)) stop("Failed to copy accepted input: ", accepted_inputs$OriginalPath[i])
}

source_manifest <- accepted_inputs %>%
  transmute(
    LogicalName,
    AcceptedOriginalPath = OriginalRelativePath,
    PackagedSourcePath = paste0("source_data/", PackageName),
    OriginalSHA256 = vapply(OriginalPath, sha256, character(1)),
    PackagedSHA256 = vapply(PackagePath, sha256, character(1)),
    HashMatch = OriginalSHA256 == PackagedSHA256
  )
stopifnot(all(source_manifest$HashMatch))
write_tsv(source_manifest, file.path(paths$manifests, "accepted_source_traceability_v1.tsv"))

read_matrix_tsv <- function(path) {
  dat <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  ids <- dat[[1]]
  mat <- as.matrix(dat[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids
  mat
}

read_beta_bundle <- function(marker) {
  list(
    scores = read_tsv(file.path(paths$source, paste0("accepted_", marker, "_cpcoa_scores.tsv")), show_col_types = FALSE),
    matrix = read_matrix_tsv(file.path(paths$source, paste0("accepted_", marker, "_bray_matrix.tsv"))),
    order = read_tsv(file.path(paths$source, paste0("accepted_", marker, "_sample_order.tsv")), show_col_types = FALSE)
  )
}

site_levels <- c("Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
sample_type_levels <- c("Sea-sand", "Water")
season_levels <- c("Spring", "Autumn")
block_levels <- c("Sea-sand / Spring", "Sea-sand / Autumn", "Water / Spring", "Water / Autumn")

sample_type_colors <- c("Sea-sand" = "#D98C2B", "Water" = "#2C7FB8")
season_colors <- c("Spring" = "#3D3D3D", "Autumn" = "#B9B9B9")
site_colors <- c(
  "Goseong" = "#7B6FD0", "Pohang" = "#4FB3BF", "Pyoseon" = "#B79F31",
  "Sehwa" = "#D87891", "Seongsan" = "#6AA84F"
)

validate_beta_bundle <- function(bundle, marker) {
  s <- bundle$scores
  m <- bundle$matrix
  o <- bundle$order
  required_score <- c(
    "SampleID", "Marker", "Site", "SampleType", "Season", "FieldReplicate", "LocationID",
    "cPCoA1", "cPCoA2", "cPCoA1PercentConstrained", "cPCoA2PercentConstrained",
    "SetID", "Formula", "Correction"
  )
  stopifnot(all(required_score %in% names(s)))
  stopifnot(nrow(s) == 60L, dplyr::n_distinct(s$SampleID) == 60L)
  marker_value <- as.character(s$Marker)
  stopifnot(all(marker_value %in% c(marker, sub("S$", "", marker))))
  stopifnot(setequal(unique(s$Site), site_levels))
  stopifnot(setequal(unique(s$SampleType), sample_type_levels))
  stopifnot(setequal(unique(s$Season), season_levels))
  stopifnot(all(count(s, SampleType)$n == 30L), all(count(s, Season)$n == 30L))
  stopifnot(dplyr::n_distinct(s$LocationID) == 15L)
  pair_count <- s %>% count(LocationID, Season)
  stopifnot(nrow(pair_count) == 30L, all(pair_count$n == 2L))
  stopifnot(all(count(s, LocationID, Season, SampleType)$n == 1L))
  stopifnot(nrow(m) == 60L, ncol(m) == 60L)
  stopifnot(identical(rownames(m), colnames(m)))
  stopifnot(setequal(rownames(m), s$SampleID))
  stopifnot(max(abs(m - t(m))) < 1e-12, max(abs(diag(m))) < 1e-12)
  stopifnot(min(m) >= 0, max(m) <= 1 + 1e-12)
  stopifnot(nrow(o) == 60L, dplyr::n_distinct(o$SampleID) == 60L)
  stopifnot(setequal(o$SampleID, s$SampleID), setequal(o$SampleID, rownames(m)))
  stopifnot(setequal(as.integer(o$Order), seq_len(60L)))
  stopifnot(all(count(o, Block)$n == 15L))
  invisible(TRUE)
}

b16 <- read_beta_bundle("16S")
b18 <- read_beta_bundle("18S")
validate_beta_bundle(b16, "16S")
validate_beta_bundle(b18, "18S")

# Use one coordinate window for both marker panels.  This preserves equal x:y
# units while giving the A and B plotting boxes the same physical height.
cpcoa_xlim <- range(c(b16$scores$cPCoA1, b18$scores$cPCoA1), finite = TRUE)
cpcoa_ylim <- range(c(b16$scores$cPCoA2, b18$scores$cPCoA2), finite = TRUE)

prepare_scores <- function(scores) {
  scores %>% mutate(
    Site = factor(Site, levels = site_levels),
    SampleType = factor(SampleType, levels = sample_type_levels),
    Season = factor(Season, levels = season_levels)
  )
}

make_cpcoa_panel <- function(scores, marker, panel_label) {
  scores <- prepare_scores(scores)
  pair_lines <- scores %>%
    select(LocationID, Season, SampleType, cPCoA1, cPCoA2) %>%
    pivot_wider(names_from = SampleType, values_from = c(cPCoA1, cPCoA2)) %>%
    filter(
      !is.na(`cPCoA1_Sea-sand`), !is.na(cPCoA1_Water),
      !is.na(`cPCoA2_Sea-sand`), !is.na(cPCoA2_Water)
    )
  stopifnot(nrow(pair_lines) == 30L)
  axis1 <- unique(scores$cPCoA1PercentConstrained)
  axis2 <- unique(scores$cPCoA2PercentConstrained)
  stopifnot(length(axis1) == 1L, length(axis2) == 1L)

  ggplot(scores, aes(cPCoA1, cPCoA2)) +
    geom_hline(yintercept = 0, linewidth = 0.28, linetype = "dotted", colour = "#858585") +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "dotted", colour = "#858585") +
    geom_segment(
      data = pair_lines,
      aes(
        x = `cPCoA1_Sea-sand`, y = `cPCoA2_Sea-sand`,
        xend = cPCoA1_Water, yend = cPCoA2_Water
      ),
      inherit.aes = FALSE, linewidth = 0.28, colour = "#BEBEBE", alpha = 0.72
    ) +
    geom_point(
      aes(fill = SampleType, shape = SampleType),
      size = 2.2, stroke = 0.45, colour = "#222222", alpha = 0.92
    ) +
    scale_fill_manual(values = sample_type_colors, drop = FALSE) +
    scale_shape_manual(values = c("Sea-sand" = 21, "Water" = 24), drop = FALSE) +
    facet_wrap(~Season, nrow = 1) +
    labs(
      title = paste0(panel_label, " ", marker),
      x = sprintf("cPCoA1 (%.1f%% constrained)", axis1),
      y = sprintf("cPCoA2 (%.1f%% constrained)", axis2),
      fill = "Sample type", shape = "Sample type"
    ) +
    coord_equal(xlim = cpcoa_xlim, ylim = cpcoa_ylim, clip = "off") +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      plot.title = element_text(size = 9.2, face = "bold", hjust = 0, margin = margin(b = 2)),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      axis.title = element_text(size = 7.5, colour = "black", face = "bold"),
      axis.text = element_text(size = 6.6, colour = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 7.1, face = "bold"),
      legend.text = element_text(size = 6.9),
      legend.key.height = unit(3.2, "mm"),
      strip.background = element_rect(fill = "#F7F7F7", colour = "#777777", linewidth = 0.3),
      strip.text = element_text(size = 7.0, face = "bold", colour = "black"),
      panel.spacing = unit(2.6, "mm"),
      plot.margin = margin(1.5, 2, 1, 1.5, unit = "mm")
    )
}

figure2 <- (make_cpcoa_panel(b16$scores, "16S", "A") |
              make_cpcoa_panel(b18$scores, "18S", "B")) +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom")

alpha_effects <- read_tsv(
  file.path(paths$source, "accepted_alpha_A72_site_synchronized_results_v2.tsv"),
  show_col_types = FALSE
) %>%
  filter(Estimand == "synchronized_SiteSeason_block_mean")

alpha_sensitivity <- read_tsv(
  file.path(paths$source, "accepted_alpha_five_contrast_family_sensitivity_v2.tsv"),
  show_col_types = FALSE
) %>%
  filter(
    Estimand == "synchronized_SiteSeason_block_mean",
    ModelID == "ALPHA_A72_BLOCK_CONTRAST"
  )

alpha_plot_data <- alpha_effects %>%
  inner_join(
    alpha_sensitivity %>%
      select(
        Marker, Metric, EffectUsed, PUsed,
        QWithinMetricFive = PAdjustedBHFiveContrasts,
        QWithinMarkerAllAlpha = PAdjustedBHWithinMarkerAllAlpha
      ),
    by = c("Marker", "Metric")
  )

stopifnot(nrow(alpha_plot_data) == 12L)
stopifnot(all(alpha_plot_data$NDisplayedBlocks == 12L))
stopifnot(all(alpha_plot_data$NIndependentSites == 7L))
stopifnot(max(abs(alpha_plot_data$EffectMeanContrast - alpha_plot_data$EffectUsed)) < 1e-9)
stopifnot(max(abs(alpha_plot_data$PExactSiteSynchronized - alpha_plot_data$PUsed)) < 1e-12)
stopifnot(all(alpha_plot_data$CILower <= alpha_plot_data$EffectMeanContrast))
stopifnot(all(alpha_plot_data$CIUpper >= alpha_plot_data$EffectMeanContrast))

metric_levels <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
metric_labels <- c(
  "Observed" = "Observed ASVs", "Chao1" = "Chao1", "ACE" = "ACE",
  "Shannon" = "Shannon", "Simpson" = "Simpson", "InvSimpson" = "Inverse Simpson"
)
marker_colors <- c("A 16S" = "#31688E", "B 18S" = "#C65D3A")

format_q <- function(x) {
  ifelse(x < 0.001, "q < 0.001", paste0("q = ", sprintf("%.3f", x)))
}

alpha_plot_data <- alpha_plot_data %>%
  mutate(
    Marker = as.character(Marker),
    MarkerPanel = factor(if_else(Marker == "16", "A 16S", "B 18S"), levels = c("A 16S", "B 18S")),
    Metric = factor(Metric, levels = metric_levels),
    MetricLabel = factor(metric_labels[as.character(Metric)], levels = unname(metric_labels[metric_levels])),
    QLabel = format_q(QWithinMetricFive)
  )

figure3 <- ggplot(alpha_plot_data, aes(colour = MarkerPanel)) +
  geom_vline(xintercept = 0, linewidth = 0.40, colour = "#777777", linetype = "dashed") +
  geom_segment(
    aes(x = CILower, xend = CIUpper, y = 1, yend = 1),
    linewidth = 0.72, lineend = "round"
  ) +
  geom_segment(aes(x = CILower, xend = CILower, y = 0.92, yend = 1.08), linewidth = 0.48) +
  geom_segment(aes(x = CIUpper, xend = CIUpper, y = 0.92, yend = 1.08), linewidth = 0.48) +
  geom_point(aes(x = EffectMeanContrast, y = 1), size = 2.7, shape = 16) +
  geom_text(
    aes(x = Inf, y = 1.30, label = QLabel),
    hjust = 1.05, vjust = 0.5, size = 2.25, colour = "#333333", family = "sans"
  ) +
  facet_grid(MarkerPanel ~ MetricLabel, scales = "free_x", switch = "y") +
  scale_colour_manual(values = marker_colors, guide = "none") +
  scale_y_continuous(limits = c(0.78, 1.42), breaks = NULL) +
  scale_x_continuous(labels = label_number(accuracy = NULL, scale_cut = cut_short_scale())) +
  labs(x = "Mean Site x Season block contrast (Sea-sand - Water)", y = NULL) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    panel.border = element_rect(fill = NA, colour = "#888888", linewidth = 0.30),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 8.2, face = "bold", margin = margin(t = 5)),
    axis.text.x = element_text(size = 6.4, colour = "black"),
    strip.background = element_rect(fill = "#F7F7F7", colour = "#777777", linewidth = 0.3),
    strip.text.x = element_text(size = 6.9, face = "bold", colour = "black"),
    strip.text.y.left = element_text(size = 8.0, face = "bold", angle = 90, colour = "black"),
    strip.placement = "outside",
    panel.spacing.x = unit(1.7, "mm"),
    panel.spacing.y = unit(2.0, "mm"),
    plot.margin = margin(2, 2, 2, 2, unit = "mm")
  )

prepare_heatmap <- function(bundle, marker, panel_label) {
  ord <- bundle$order %>% arrange(Order)
  stopifnot(identical(as.integer(ord$Order), seq_len(60L)))
  ids <- ord$SampleID
  mat <- bundle$matrix[ids, ids, drop = FALSE]
  ord <- ord %>% mutate(
    SampleType = factor(SampleType, levels = sample_type_levels),
    Season = factor(Season, levels = season_levels),
    Site = factor(Site, levels = site_levels),
    Block = factor(Block, levels = block_levels)
  )
  ann_df <- data.frame(
    `Sample type` = ord$SampleType,
    Season = ord$Season,
    Site = ord$Site,
    check.names = FALSE,
    row.names = ids
  )
  ann_colors <- list(`Sample type` = sample_type_colors, Season = season_colors, Site = site_colors)
  top_ha <- HeatmapAnnotation(
    df = ann_df, col = ann_colors, show_legend = FALSE,
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontfamily = "sans", fontsize = 6.0, fontface = "bold"),
    simple_anno_size = unit(2.35, "mm"), gap = unit(0.45, "mm")
  )
  left_ha <- rowAnnotation(
    df = ann_df, col = ann_colors, show_legend = FALSE,
    show_annotation_name = FALSE,
    simple_anno_size = unit(2.35, "mm"), gap = unit(0.45, "mm")
  )
  heat_fun <- colorRamp2(
    c(0, 0.25, 0.50, 0.75, 1),
    c("#F7FCFD", "#D9F0E8", "#A6BDDB", "#756BB1", "#2D004B")
  )
  ht <- Heatmap(
    mat,
    name = "Bray-Curtis dissimilarity",
    col = heat_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = ids, column_order = ids,
    row_split = ord$Block, column_split = ord$Block,
    row_gap = unit(0.72, "mm"), column_gap = unit(0.72, "mm"),
    row_title = NULL,
    column_title = paste0(panel_label, " ", marker),
    column_title_gp = gpar(fontfamily = "sans", fontsize = 9, fontface = "bold"),
    column_title_rot = 0,
    column_names_gp = gpar(fontfamily = "sans", fontsize = 5.7, fontface = "bold"),
    show_row_names = FALSE, show_column_names = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    rect_gp = gpar(col = NA), border = TRUE,
    top_annotation = top_ha, left_annotation = left_ha
  )
  grob <- grid.grabExpr(
    ComplexHeatmap::draw(ht, padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm")),
    wrap = TRUE
  )
  list(grob = grob, matrix = mat, order = ord, heat_fun = heat_fun)
}

h16 <- prepare_heatmap(b16, "16S", "A")
h18 <- prepare_heatmap(b18, "18S", "B")

legend_pack <- packLegend(
  Legend(
    title = "Bray-Curtis\ndissimilarity", col_fun = h16$heat_fun,
    at = c(0, 0.25, 0.50, 0.75, 1),
    direction = "horizontal",
    title_gp = gpar(fontfamily = "sans", fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6.5),
    legend_width = unit(29, "mm"), grid_height = unit(3.2, "mm")
  ),
  Legend(
    title = "Sample type", labels = names(sample_type_colors),
    legend_gp = gpar(fill = unname(sample_type_colors), col = NA),
    title_gp = gpar(fontfamily = "sans", fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6.5), grid_width = unit(3.2, "mm"), grid_height = unit(3.2, "mm"),
    nrow = 1
  ),
  Legend(
    title = "Season", labels = names(season_colors),
    legend_gp = gpar(fill = unname(season_colors), col = NA),
    title_gp = gpar(fontfamily = "sans", fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6.5), grid_width = unit(3.2, "mm"), grid_height = unit(3.2, "mm"),
    nrow = 1
  ),
  Legend(
    title = "Site", labels = names(site_colors),
    legend_gp = gpar(fill = unname(site_colors), col = NA),
    title_gp = gpar(fontfamily = "sans", fontsize = 7.0, fontface = "bold"),
    labels_gp = gpar(fontfamily = "sans", fontsize = 6.5), grid_width = unit(3.2, "mm"), grid_height = unit(3.2, "mm"),
    nrow = 1
  ),
  direction = "horizontal", gap = unit(4.0, "mm"), max_width = unit(176, "mm")
)

legend_grob <- grid.grabExpr(
  ComplexHeatmap::draw(
    legend_pack,
    x = unit(0.50, "npc"), y = unit(0.50, "npc"), just = c("center", "center")
  ),
  wrap = TRUE
)

heat_panels <- wrap_plots(
  list(
    wrap_elements(full = h16$grob),
    wrap_elements(full = h18$grob)
  ),
  nrow = 1,
  widths = c(1, 1)
)
figure_s4 <- wrap_plots(
  list(heat_panels, wrap_elements(full = legend_grob)),
  ncol = 1,
  heights = c(1, 0.23)
)

figure2_source <- bind_rows(
  b16$scores %>% mutate(Marker = "16S"),
  b18$scores %>% mutate(Marker = "18S")
) %>% arrange(Marker, Season, Site, SampleType, FieldReplicate)
write_tsv(figure2_source, file.path(paths$source, "Figure_2_16S_18S_beta_diversity_v1_plot_data.tsv"))

write_tsv(
  alpha_plot_data %>%
    transmute(
      Marker, Metric = as.character(Metric), Estimand,
      EffectSeaSandMinusWater = EffectMeanContrast,
      DescriptiveT95CILower = CILower,
      DescriptiveT95CIUpper = CIUpper,
      ExactSiteSynchronizedP = PExactSiteSynchronized,
      BHQWithinFivePrespecifiedContrasts = QWithinMetricFive,
      BHQWithinAll30AlphaTests = QWithinMarkerAllAlpha,
      DisplayedSiteSeasonBlocks = NDisplayedBlocks,
      IndependentSites = NIndependentSites,
      ExactPatterns, Pairing
    ),
  file.path(paths$source, "Figure_3_16S_18S_alpha_effect_sizes_v1_plot_data.tsv")
)

save_bundle <- function(plot_object, stem, width_mm, height_mm) {
  out <- c(
    pdf = file.path(paths$figures, paste0(stem, ".pdf")),
    tiff = file.path(paths$figures, paste0(stem, ".tiff")),
    svg = file.path(paths$figures, paste0(stem, ".svg")),
    preview = file.path(paths$previews, paste0(stem, "_direct_preview.png")),
    pdf_render = file.path(paths$previews, paste0(stem, "_pdf_render_page_1.png"))
  )
  ggsave(out[["pdf"]], plot_object, device = cairo_pdf, width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(out[["svg"]], plot_object, device = svglite::svglite, width = width_mm, height = height_mm, units = "mm", bg = "white")
  ragg::agg_tiff(
    out[["tiff"]], width = width_mm, height = height_mm, units = "mm",
    res = 600, compression = "lzw", background = "white", scaling = 1
  )
  print(plot_object)
  dev.off()
  ragg::agg_png(
    out[["preview"]], width = width_mm, height = height_mm, units = "mm",
    res = 240, background = "white", scaling = 1
  )
  print(plot_object)
  dev.off()
  render_template <- file.path(paths$previews, paste0(stem, "_pdf_render_page_%d.%s"))
  pdftools::pdf_convert(
    out[["pdf"]], format = "png", pages = 1L, dpi = 180,
    filenames = render_template, verbose = FALSE
  )
  info <- pdftools::pdf_info(out[["pdf"]])
  if (info$pages != 1L) stop("Expected a one-page PDF for ", stem)
  tibble(
    Stem = stem,
    WidthMM = width_mm,
    HeightMM = height_mm,
    ExpectedTIFFWidthPixels = floor(width_mm / 25.4 * 600),
    ExpectedTIFFHeightPixels = floor(height_mm / 25.4 * 600),
    TIFFResolutionDPI = 600L,
    TIFFCompression = "LZW (ragg export argument)",
    PDFPages = info$pages,
    PDFVersion = info$version,
    PDFBytes = file.info(out[["pdf"]])$size,
    TIFFBytes = file.info(out[["tiff"]])$size,
    SVGBytes = file.info(out[["svg"]])$size,
    PDFSHA256 = sha256(out[["pdf"]]),
    TIFFSHA256 = sha256(out[["tiff"]]),
    SVGSHA256 = sha256(out[["svg"]])
  )
}

file_qa <- bind_rows(
  save_bundle(figure2, "Figure_2_16S_18S_beta_diversity_v1", 180, 100),
  save_bundle(figure3, "Figure_3_16S_18S_alpha_effect_sizes_v1", 180, 105),
  save_bundle(figure_s4, "Supplementary_Figure_S4_16S_18S_beta_heatmaps_v1", 180, 105)
)
write_tsv(file_qa, file.path(paths$qc, "technical_file_qa_v1.tsv"))

beta_perm <- read_tsv(
  file.path(paths$source, "accepted_beta_hierarchical_permanova_v2.tsv"),
  show_col_types = FALSE
)
beta_disp <- read_tsv(
  file.path(paths$source, "accepted_beta_hierarchical_permdisp_v2.tsv"),
  show_col_types = FALSE
)
beta_perm_primary <- beta_perm %>% filter(Surface == "rarefied_primary")
beta_disp_primary <- beta_disp %>% filter(Surface == "rarefied_primary")
stopifnot(nrow(beta_perm_primary) == 8L, nrow(beta_disp_primary) == 8L)
stopifnot(all(beta_perm_primary$PAdjustedBHFourHierarchicalTargets >= 0.05))
stopifnot(all(beta_disp_primary$PAdjustedBHFourHierarchicalTargets >= 0.05))

qa_table <- bind_rows(
  tibble(Test = "16S cPCoA sample count", Observed = as.character(nrow(b16$scores)), Expected = "60", Status = "PASS"),
  tibble(Test = "18S cPCoA sample count", Observed = as.character(nrow(b18$scores)), Expected = "60", Status = "PASS"),
  tibble(Test = "16S matched contemporaneous segment count", Observed = as.character(nrow(count(b16$scores, LocationID, Season))), Expected = "30", Status = "PASS"),
  tibble(Test = "18S matched contemporaneous segment count", Observed = as.character(nrow(count(b18$scores, LocationID, Season))), Expected = "30", Status = "PASS"),
  tibble(Test = "16S Bray matrix dimension", Observed = paste(dim(b16$matrix), collapse = " x "), Expected = "60 x 60", Status = "PASS"),
  tibble(Test = "18S Bray matrix dimension", Observed = paste(dim(b18$matrix), collapse = " x "), Expected = "60 x 60", Status = "PASS"),
  tibble(Test = "Alpha plotted rows", Observed = as.character(nrow(alpha_plot_data)), Expected = "12 (2 markers x 6 metrics)", Status = "PASS"),
  tibble(Test = "Alpha displayed block count", Observed = paste(sort(unique(alpha_plot_data$NDisplayedBlocks)), collapse = ","), Expected = "12", Status = "PASS"),
  tibble(Test = "Alpha independent Site count", Observed = paste(sort(unique(alpha_plot_data$NIndependentSites)), collapse = ","), Expected = "7", Status = "PASS"),
  tibble(Test = "Minimum rarefied-primary hierarchical PERMANOVA BH q", Observed = format(min(beta_perm_primary$PAdjustedBHFourHierarchicalTargets), digits = 8), Expected = ">= 0.05", Status = "PASS"),
  tibble(Test = "Minimum rarefied-primary matched PERMDISP BH q", Observed = format(min(beta_disp_primary$PAdjustedBHFourHierarchicalTargets), digits = 8), Expected = ">= 0.05", Status = "PASS"),
  tibble(Test = "Accepted source hash identity", Observed = as.character(sum(source_manifest$HashMatch)), Expected = as.character(nrow(source_manifest)), Status = "PASS")
)
write_tsv(qa_table, file.path(paths$qc, "numeric_source_validation_v1.tsv"))

qa_json <- list(
  status = "PASS",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  R_version = R.version.string,
  numerical_reanalysis = FALSE,
  figure2 = list(
    panel_map = c("A 16S", "B 18S"),
    sample_n_per_marker = 60L,
    locations = 15L,
    matched_segments_per_marker = 30L,
    axis_percent = list(
      `16S` = c(unique(b16$scores$cPCoA1PercentConstrained), unique(b16$scores$cPCoA2PercentConstrained)),
      `18S` = c(unique(b18$scores$cPCoA1PercentConstrained), unique(b18$scores$cPCoA2PercentConstrained))
    )
  ),
  figure3 = list(
    plotted_rows = nrow(alpha_plot_data),
    displayed_samples = 72L,
    contemporaneous_pairs = 36L,
    site_season_blocks = 12L,
    independent_sites = 7L,
    contrast = "Sea-sand minus Water",
    q_labels = "BH within the five prespecified contrasts for each marker and metric",
    stars = FALSE
  ),
  supplementary_figure_s4 = list(
    panel_map = c("A 16S", "B 18S"),
    matrix_dimension_per_marker = c(60L, 60L),
    accepted_order_preserved = TRUE
  ),
  inference_guardrail = list(
    minimum_primary_permanova_q = min(beta_perm_primary$PAdjustedBHFourHierarchicalTargets),
    minimum_primary_permdisp_q = min(beta_disp_primary$PAdjustedBHFourHierarchicalTargets),
    statement = "No rarefied-primary Site-synchronized hierarchical beta target or matched PERMDISP target passed BH correction."
  )
)
write_json(qa_json, file.path(paths$qc, "final_hierarchy_figures_validation_v1.json"), pretty = TRUE, auto_unbox = TRUE)

write_package_manifest()
message("PASS: three presentation-only hierarchy figures were generated from accepted TSVs.")
