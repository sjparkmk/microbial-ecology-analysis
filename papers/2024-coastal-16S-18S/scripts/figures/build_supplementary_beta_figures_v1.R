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
})

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/")
script_dir <- dirname(script_path)
package_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/")
project_root <- normalizePath(file.path(package_dir, "..", "..", ".."), winslash = "/")

paths <- list(
  scores = file.path(project_root, "revision", "re_analysis", "12_original_pipeline_reanalysis_v1", "core_v1", "source_data", "core_C04_partial_cpcoa_scores_source_v1.tsv"),
  models = file.path(project_root, "revision", "re_analysis", "12_original_pipeline_reanalysis_v1", "core_v1", "source_data", "core_C04_partial_cpcoa_models_source_v1.tsv"),
  membership = file.path(project_root, "revision", "re_analysis", "00_design_v4", "analysis_set_membership_v4.tsv"),
  pairwise = file.path(project_root, "revision", "re_analysis", "10_figure_source_data_v4", "source_data", "DS5_bray_pairwise_all_subsets_source_data_v1.tsv"),
  permanova = file.path(project_root, "revision", "re_analysis", "06_beta_diversity_v4", "results", "beta_permanova_restricted_v4.tsv"),
  permdisp = file.path(project_root, "revision", "re_analysis", "06_beta_diversity_v4", "results", "beta_permdisp_restricted_v4.tsv"),
  source_dir = file.path(package_dir, "source_data"),
  figure_dir = file.path(package_dir, "figures"),
  qc_dir = file.path(package_dir, "qc"),
  delivery_dir = file.path(project_root, "revision", "rev_ms", "figures")
)

invisible(lapply(paths[c("source_dir", "figure_dir", "qc_dir", "delivery_dir")], dir.create,
                 recursive = TRUE, showWarnings = FALSE))

required_inputs <- unlist(paths[c("scores", "models", "membership", "pairwise", "permanova", "permdisp")])
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required inputs: ", paste(missing_inputs, collapse = "; "))

read_tsv_clean <- function(path) {
  suppressMessages(readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ))
}

scores <- read_tsv_clean(paths$scores)
models <- read_tsv_clean(paths$models)
membership <- read_tsv_clean(paths$membership)
pairwise <- read_tsv_clean(paths$pairwise)
permanova <- read_tsv_clean(paths$permanova)
permdisp <- read_tsv_clean(paths$permdisp)

scores <- scores %>%
  mutate(
    cPCoA1 = as.numeric(cPCoA1),
    cPCoA2 = as.numeric(cPCoA2),
    cPCoA1PercentConstrained = as.numeric(cPCoA1PercentConstrained),
    cPCoA2PercentConstrained = as.numeric(cPCoA2PercentConstrained),
    FieldReplicate = as.integer(FieldReplicate)
  )

a60 <- membership %>%
  filter(SetID == "A60_COMPLETE_FACTORIAL", Included %in% c(TRUE, "TRUE")) %>%
  transmute(
    SampleID, Site, GeologicalHabitat, SampleType,
    Season, FieldReplicate = as.integer(FieldReplicate), LocationID
  )

stopifnot(nrow(a60) == 60L, n_distinct(a60$SampleID) == 60L)
stopifnot(setequal(unique(scores$SampleID), a60$SampleID))

site_levels <- c("Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
sample_type_levels <- c("Sea-sand", "Water")
season_levels <- c("Spring", "Autumn")
block_levels <- c(
  "Sea-sand / Spring", "Sea-sand / Autumn",
  "Water / Spring", "Water / Autumn"
)

sample_type_colors <- c("Sea-sand" = "#D98C2B", "Water" = "#2C7FB8")
season_colors <- c("Spring" = "#3D3D3D", "Autumn" = "#B9B9B9")
site_colors <- c(
  "Goseong" = "#7B6FD0", "Pohang" = "#4FB3BF", "Pyoseon" = "#B79F31",
  "Sehwa" = "#D87891", "Seongsan" = "#6AA84F"
)

a60 <- a60 %>%
  mutate(
    Site = factor(Site, levels = site_levels),
    SampleType = factor(SampleType, levels = sample_type_levels),
    Season = factor(Season, levels = season_levels),
    Block = factor(paste(SampleType, Season, sep = " / "), levels = block_levels)
  )

block_counts <- a60 %>% count(Block, .drop = FALSE)
stopifnot(all(block_counts$n == 15L))
stopifnot(n_distinct(a60$LocationID) == 15L)

reconstruct_bray <- function(marker) {
  value_col <- if (marker == "16S") "Bray16S" else "Bray18S"
  pw <- pairwise %>%
    filter(
      AnalysisSubset == "A75",
      SampleID1 %in% a60$SampleID,
      SampleID2 %in% a60$SampleID
    ) %>%
    select(SampleID1, SampleID2, Distance = all_of(value_col)) %>%
    mutate(Distance = as.numeric(Distance))
  stopifnot(nrow(pw) == choose(60L, 2L), all(is.finite(pw$Distance)))
  ids <- a60$SampleID
  mat <- matrix(0, nrow = length(ids), ncol = length(ids), dimnames = list(ids, ids))
  mat[cbind(match(pw$SampleID1, ids), match(pw$SampleID2, ids))] <- pw$Distance
  mat[cbind(match(pw$SampleID2, ids), match(pw$SampleID1, ids))] <- pw$Distance
  stopifnot(max(abs(mat - t(mat))) < 1e-12, all(diag(mat) == 0))
  stopifnot(min(mat) >= 0, max(mat) <= 1 + 1e-12)
  list(matrix = mat, pairs = pw)
}

order_within_blocks <- function(mat, metadata) {
  ordered_ids <- character()
  for (block in block_levels) {
    ids <- metadata$SampleID[metadata$Block == block]
    stopifnot(length(ids) == 15L)
    hc <- hclust(as.dist(mat[ids, ids, drop = FALSE]), method = "average")
    ordered_ids <- c(ordered_ids, ids[hc$order])
  }
  stopifnot(length(ordered_ids) == 60L, n_distinct(ordered_ids) == 60L)
  ordered_ids
}

make_cpcoa <- function(marker, score_df, model_row, stat_row, disp_row) {
  score_df <- score_df %>%
    mutate(
      SampleType = factor(SampleType, levels = sample_type_levels),
      Season = factor(Season, levels = season_levels)
    )

  pair_lines <- score_df %>%
    select(LocationID, Season, SampleType, cPCoA1, cPCoA2) %>%
    tidyr::pivot_wider(names_from = SampleType, values_from = c(cPCoA1, cPCoA2)) %>%
    filter(
      !is.na(`cPCoA1_Sea-sand`), !is.na(cPCoA1_Water),
      !is.na(`cPCoA2_Sea-sand`), !is.na(cPCoA2_Water)
    )
  stopifnot(nrow(pair_lines) == 30L)

  axis1 <- unique(score_df$cPCoA1PercentConstrained)
  axis2 <- unique(score_df$cPCoA2PercentConstrained)
  stopifnot(length(axis1) == 1L, length(axis2) == 1L)

  ggplot(score_df, aes(cPCoA1, cPCoA2)) +
    geom_hline(yintercept = 0, linewidth = 0.28, linetype = "dotted", colour = "#858585") +
    geom_vline(xintercept = 0, linewidth = 0.28, linetype = "dotted", colour = "#858585") +
    geom_segment(
      data = pair_lines,
      aes(
        x = `cPCoA1_Sea-sand`, y = `cPCoA2_Sea-sand`,
        xend = cPCoA1_Water, yend = cPCoA2_Water
      ),
      inherit.aes = FALSE, linewidth = 0.32, colour = "#C7C7C7", alpha = 0.7
    ) +
    geom_point(
      aes(fill = SampleType, shape = SampleType),
      size = 2.65, stroke = 0.5, colour = "#262626", alpha = 0.9
    ) +
    scale_fill_manual(values = sample_type_colors, drop = FALSE) +
    scale_shape_manual(values = c("Sea-sand" = 21, "Water" = 24), drop = FALSE) +
    facet_wrap(~Season, nrow = 1) +
    labs(
      x = sprintf("cPCoA1 (%.1f%% of constrained variation)", axis1),
      y = sprintf("cPCoA2 (%.1f%% of constrained variation)", axis2),
      fill = "Sample type", shape = "Sample type"
    ) +
    coord_equal(clip = "off") +
    theme_classic(base_size = 8, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.32, colour = "black"),
      axis.title = element_text(size = 8.2, colour = "black"),
      axis.text = element_text(size = 7.2, colour = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 7.3, face = "bold"),
      legend.text = element_text(size = 7.1),
      legend.key.height = unit(3.5, "mm"),
      legend.spacing.x = unit(2.2, "mm"),
      strip.background = element_rect(fill = "white", colour = "#777777", linewidth = 0.3),
      strip.text = element_text(size = 7.5, face = "bold", colour = "black"),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(2, 4, 1, 2, unit = "mm")
    )
}

make_heatmap_grob <- function(marker, mat, metadata, ordered_ids) {
  meta_ord <- metadata[match(ordered_ids, metadata$SampleID), , drop = FALSE]
  mat_ord <- mat[ordered_ids, ordered_ids, drop = FALSE]
  stopifnot(identical(rownames(mat_ord), colnames(mat_ord)))
  stopifnot(identical(meta_ord$SampleID, rownames(mat_ord)))

  ann_df <- data.frame(
    `Sample type` = factor(meta_ord$SampleType, levels = sample_type_levels),
    Season = factor(meta_ord$Season, levels = season_levels),
    Site = factor(meta_ord$Site, levels = site_levels),
    check.names = FALSE,
    row.names = meta_ord$SampleID
  )

  ann_colors <- list(
    `Sample type` = sample_type_colors,
    Season = season_colors,
    Site = site_colors
  )

  top_ha <- HeatmapAnnotation(
    df = ann_df,
    col = ann_colors,
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontfamily = "sans", fontsize = 6.3, fontface = "bold"),
    annotation_legend_param = list(
      `Sample type` = list(title = "Sample type", title_gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"), labels_gp = gpar(fontfamily = "sans", fontsize = 6.5)),
      Season = list(title = "Season", title_gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"), labels_gp = gpar(fontfamily = "sans", fontsize = 6.5)),
      Site = list(title = "Site", title_gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"), labels_gp = gpar(fontfamily = "sans", fontsize = 6.5), nrow = 5)
    ),
    simple_anno_size = unit(2.6, "mm"),
    gap = unit(0.5, "mm")
  )

  left_ha <- rowAnnotation(
    df = ann_df,
    col = ann_colors,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    annotation_name_gp = gpar(fontfamily = "sans", fontsize = 6.3, fontface = "bold"),
    simple_anno_size = unit(2.6, "mm"),
    gap = unit(0.5, "mm")
  )

  heat_colors <- circlize::colorRamp2(
    c(0, 0.25, 0.50, 0.75, 1),
    c("#F7FCFD", "#D9F0E8", "#A6BDDB", "#756BB1", "#2D004B")
  )

  ht <- Heatmap(
    mat_ord,
    name = "Bray-Curtis\ndissimilarity",
    col = heat_colors,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_order = ordered_ids,
    column_order = ordered_ids,
    row_split = factor(meta_ord$Block, levels = block_levels),
    column_split = factor(meta_ord$Block, levels = block_levels),
    row_gap = unit(0.8, "mm"),
    column_gap = unit(0.8, "mm"),
    row_title = NULL,
    column_title_gp = gpar(fontfamily = "sans", fontsize = 6.2, fontface = "bold"),
    column_title_rot = 0,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    rect_gp = gpar(col = NA),
    border = TRUE,
    top_annotation = top_ha,
    left_annotation = left_ha,
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.50, 0.75, 1),
      labels = c("0", "0.25", "0.50", "0.75", "1"),
      title_gp = gpar(fontfamily = "sans", fontsize = 7, fontface = "bold"),
      labels_gp = gpar(fontfamily = "sans", fontsize = 6.5),
      legend_height = unit(27, "mm"),
      grid_width = unit(3.2, "mm")
    )
  )

  grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = FALSE,
      padding = unit(c(1, 1, 1, 1), "mm")
    ),
    wrap = TRUE
  )

  list(grob = grob, matrix = mat_ord, metadata = meta_ord)
}

save_plot_bundle <- function(plot_obj, stem) {
  width_mm <- 180
  height_mm <- 230
  pdf_file <- file.path(paths$delivery_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(paths$figure_dir, paste0(stem, ".svg"))
  tiff_file <- file.path(paths$delivery_dir, paste0(stem, ".tiff"))
  png_file <- file.path(paths$qc_dir, paste0(stem, "_preview.png"))

  ggsave(pdf_file, plot_obj, device = cairo_pdf, width = width_mm, height = height_mm, units = "mm", bg = "white")
  ggsave(svg_file, plot_obj, device = svglite::svglite, width = width_mm, height = height_mm, units = "mm", bg = "white")

  ragg::agg_tiff(
    tiff_file, width = width_mm, height = height_mm, units = "mm",
    res = 600, compression = "lzw", background = "white", scaling = 1
  )
  print(plot_obj)
  dev.off()

  ragg::agg_png(
    png_file, width = width_mm, height = height_mm, units = "mm",
    res = 240, background = "white", scaling = 1
  )
  print(plot_obj)
  dev.off()

  c(pdf = pdf_file, svg = svg_file, tiff = tiff_file, preview = png_file)
}

qa_rows <- list()
output_paths <- list()

for (marker in c("16S", "18S")) {
  score_df <- scores %>% filter(Marker == marker) %>% arrange(match(SampleID, a60$SampleID))
  stopifnot(nrow(score_df) == 60L, n_distinct(score_df$SampleID) == 60L)

  model_row <- models %>% filter(Marker == marker, SetID == "A60_COMPLETE_FACTORIAL")
  stat_row <- permanova %>% filter(
    Marker == marker,
    Surface == "rarefied_primary",
    ModelID == "BETA_TYPE_CAMPAIGN",
    SetID == "A60_COMPLETE_FACTORIAL",
    Term == "SampleType:Season",
    InferenceReportable %in% c(TRUE, "TRUE")
  )
  disp_row <- permdisp %>% filter(
    Marker == marker,
    Surface == "rarefied_primary",
    ModelID == "BETA_TYPE_CAMPAIGN",
    SetID == "A60_COMPLETE_FACTORIAL",
    TargetFactor == "SampleType:Season",
    InferenceReportable %in% c(TRUE, "TRUE")
  )
  stopifnot(nrow(model_row) == 1L, nrow(stat_row) == 1L, nrow(disp_row) == 1L)

  bray <- reconstruct_bray(marker)
  ordered_ids <- order_within_blocks(bray$matrix, a60)

  p_a <- make_cpcoa(marker, score_df, model_row, stat_row, disp_row)
  h <- make_heatmap_grob(marker, bray$matrix, a60, ordered_ids)
  p_b <- patchwork::wrap_elements(full = h$grob)

  final_plot <- (p_a / p_b) +
    plot_layout(heights = c(0.82, 1.38)) +
    plot_annotation(
      tag_levels = list(c("(a)", "(b)")),
      theme = theme(
        plot.tag = element_text(family = "sans", face = "bold", size = 11, colour = "black"),
        plot.tag.position = c(0.005, 0.995),
        plot.margin = margin(2, 2, 2, 2, unit = "mm")
      )
    )

  supp_no <- if (marker == "16S") "S3" else "S4"
  stem <- paste0("Supplementary_Figure_", supp_no, "_", marker, "_beta_diversity_v1")
  output_paths[[marker]] <- save_plot_bundle(final_plot, stem)

  score_out <- score_df %>%
    select(
      SampleID, Marker, Site, GeologicalHabitat, SampleType, Season,
      FieldReplicate, LocationID, cPCoA1, cPCoA2,
      cPCoA1PercentConstrained, cPCoA2PercentConstrained,
      SetID, Formula, Correction
    )
  readr::write_tsv(score_out, file.path(paths$source_dir, paste0(stem, "_cpcoa_scores.tsv")))

  matrix_out <- tibble::as_tibble(h$matrix, rownames = "SampleID")
  readr::write_tsv(matrix_out, file.path(paths$source_dir, paste0(stem, "_bray_matrix.tsv")))

  order_out <- h$metadata %>%
    mutate(Order = seq_len(n())) %>%
    select(Order, SampleID, Site, GeologicalHabitat, SampleType, Season, FieldReplicate, LocationID, Block)
  readr::write_tsv(order_out, file.path(paths$source_dir, paste0(stem, "_sample_order.tsv")))

  stat_out <- tibble(
    Marker = marker,
    SetID = "A60_COMPLETE_FACTORIAL",
    SampleN = 60L,
    Distance = "Bray-Curtis",
    Normalization = if (marker == "16S") "fixed-depth rarefied counts; 20000 reads per sample" else "fixed-depth rarefied counts; 50000 reads per sample; PR2 v5.1.1-primary 18S",
    Formula = "Bray ~ SampleType * Season + Condition(Site)",
    Correction = "Lingoes",
    ModelR2 = as.numeric(model_row$R2),
    ModelAdjustedR2 = as.numeric(model_row$AdjustedR2),
    Axis1PercentConstrained = unique(score_df$cPCoA1PercentConstrained),
    Axis2PercentConstrained = unique(score_df$cPCoA2PercentConstrained),
    InteractionR2 = as.numeric(stat_row$R2),
    InteractionF = as.numeric(stat_row$F),
    InteractionP = as.numeric(stat_row$PReportable),
    InteractionBHAdjustedP = as.numeric(stat_row$PAdjustedBHTargetWithinMarkerSurface),
    PERMDISPF = as.numeric(disp_row$F),
    PERMDISPP = as.numeric(disp_row$PReportable),
    PERMDISPBHAdjustedP = as.numeric(disp_row$PAdjustedBHTargetWithinMarkerSurface),
    Interpretation = "Ordination and heatmap are descriptive; inference comes from LocationID-restricted PERMANOVA interpreted with matched PERMDISP."
  )
  readr::write_tsv(stat_out, file.path(paths$source_dir, paste0(stem, "_statistics.tsv")))

  qa_rows[[marker]] <- tibble(
    Marker = marker,
    SampleN = nrow(score_df),
    UniqueSampleN = n_distinct(score_df$SampleID),
    PairwiseN = nrow(bray$pairs),
    MatrixRows = nrow(h$matrix),
    MatrixColumns = ncol(h$matrix),
    SymmetryMaxAbsDifference = max(abs(h$matrix - t(h$matrix))),
    DiagonalMaxAbs = max(abs(diag(h$matrix))),
    DistanceMinimum = min(h$matrix),
    DistanceMaximum = max(h$matrix),
    SeaSandN = sum(score_df$SampleType == "Sea-sand"),
    WaterN = sum(score_df$SampleType == "Water"),
    SpringN = sum(score_df$Season == "Spring"),
    AutumnN = sum(score_df$Season == "Autumn"),
    LocationIDN = n_distinct(score_df$LocationID),
    OutputPDF = output_paths[[marker]][["pdf"]],
    OutputSVG = output_paths[[marker]][["svg"]],
    OutputTIFF = output_paths[[marker]][["tiff"]],
    PreviewPNG = output_paths[[marker]][["preview"]]
  )
}

qa <- bind_rows(qa_rows)
readr::write_tsv(qa, file.path(paths$qc_dir, "Supplementary_Figures_S3_S4_numeric_validation_v1.tsv"))

qa_json <- list(
  status = "PASS",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  backend = paste0("R ", getRversion()),
  checks = list(
    sample_population = "60 unique A60 samples per marker",
    pairwise_count = "1770 unique off-diagonal pairs per marker",
    matrix = "60 x 60, symmetric, zero diagonal, Bray-Curtis range 0-1",
    factor_balance = "30 Sea-sand, 30 Water, 30 Spring, 30 Autumn; 15 per Sample type x Season block",
    sites = site_levels,
    excluded_incomplete_sites = c("Dokdo", "Ulleungdo"),
    imputation = "none"
  ),
  output_paths = output_paths
)
jsonlite::write_json(
  qa_json,
  file.path(paths$qc_dir, "Supplementary_Figures_S3_S4_validation_v1.json"),
  pretty = TRUE, auto_unbox = TRUE
)

cat("Supplementary Figures S3-S4 generated successfully.\n")
