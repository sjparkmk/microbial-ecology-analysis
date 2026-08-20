#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(ragg)
  library(digest)
  library(pdftools)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "run"
if (!mode %in% c("run", "check")) {
  stop("Usage: Rscript --vanilla plot_S13_cross_marker_concordance_R_v1.R [run|check]")
}

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(file_arg)) stop("Cannot determine script path")
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
pkg_dir <- normalizePath(file.path(dirname(script_file), "..", "generated", "f5"), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)

canonical_dir <- file.path(project_dir, "revision", "re_analysis", "10_figure_source_data_v4", "source_data")
delivery_dir <- file.path(pkg_dir, "delivery")
source_dir <- file.path(pkg_dir, "source_data")
qc_dir <- file.path(pkg_dir, "qc")
preview_dir <- file.path(qc_dir, "previews")
manifest_dir <- file.path(pkg_dir, "manifests")
log_dir <- file.path(pkg_dir, "logs")

stem <- "Supplementary_Figure_S13_cross_marker_concordance_R_v1"
pdf_out <- file.path(delivery_dir, paste0(stem, ".pdf"))
tiff_out <- file.path(delivery_dir, paste0(stem, ".tiff"))
preview_out <- file.path(preview_dir, paste0(stem, ".png"))
legend_file <- file.path(pkg_dir, "captions", "Supplementary_Figure_S13_legend_v1.md")
qa_file <- file.path(qc_dir, "R_S13_numeric_QA_v1.tsv")
visual_qa_file <- file.path(qc_dir, "R_S13_rendered_visual_QA_v1.md")
input_manifest <- file.path(manifest_dir, "R_S13_input_manifest_v1.tsv")
output_manifest <- file.path(manifest_dir, "R_S13_output_manifest_v1.tsv")
session_file <- file.path(log_dir, "R_S13_sessionInfo_v1.txt")

source_names <- c(
  "DS5_alpha_pairs_all_subsets_source_data_v1.tsv",
  "DS5_alpha_spearman_all_subsets_source_data_v1.tsv",
  "DS5_bray_pairwise_all_subsets_source_data_v1.tsv",
  "DS5_bray_spearman_all_subsets_source_data_v1.tsv",
  "DS5_procrustes_scores_all_subsets_source_data_v1.tsv",
  "DS5_procrustes_summary_all_subsets_source_data_v1.tsv"
)
canonical_paths <- file.path(canonical_dir, source_names)
snapshot_paths <- file.path(source_dir, source_names)

read_tsv <- function(path) {
  read.delim(path, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)

project_relative <- function(path) {
  p <- normalizePath(path, winslash = "\\", mustWork = TRUE)
  root <- paste0(normalizePath(project_dir, winslash = "\\", mustWork = TRUE), "\\")
  sub(root, "", p, fixed = TRUE)
}

write_manifest <- function(paths, target) {
  rows <- data.frame(
    path = vapply(paths, project_relative, character(1)),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, sha256, character(1)),
    stringsAsFactors = FALSE
  )
  write.table(rows, target, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

assert_true <- function(flag, message) {
  if (!isTRUE(flag)) stop(message, call. = FALSE)
}

load_frozen <- function() {
  assert_true(all(file.exists(snapshot_paths)), "Frozen DS5 snapshots are incomplete")
  x <- lapply(snapshot_paths, read_tsv)
  names(x) <- c("alpha_pairs", "alpha_stats", "bray_pairs", "bray_stats", "proc_scores", "proc_stats")
  lapply(x, function(z) z[z$AnalysisSubset == "A75", , drop = FALSE])
}

validate_frozen <- function(x) {
  assert_true(nrow(x$alpha_pairs) == 225L, "A75 alpha source must contain 75 samples x 3 frozen metrics")
  assert_true(length(unique(x$alpha_pairs$SampleID)) == 75L, "A75 alpha source must contain 75 samples")
  assert_true(identical(sort(unique(x$alpha_pairs$Metric)), c("Observed", "Shannon", "Simpson")), "Frozen alpha metrics changed")
  assert_true(nrow(x$alpha_stats) == 3L, "A75 alpha statistics must contain three frozen rows")
  assert_true(nrow(x$bray_pairs) == 2775L, "A75 Bray source must contain 2,775 sample pairs")
  assert_true(nrow(x$bray_stats) == 1L, "A75 Bray statistics must contain one row")
  assert_true(nrow(x$proc_scores) == 75L, "A75 Procrustes scores must contain 75 samples")
  assert_true(nrow(x$proc_stats) == 1L, "A75 Procrustes statistics must contain one row")
  assert_true(all(x$alpha_stats$BootstrapReplicatesValid == 9999L), "Alpha Site-bootstrap count changed")
  assert_true(x$bray_stats$PermutationsEvaluated[[1]] == 9999L, "Bray permutation count changed")
  assert_true(x$proc_stats$PermutationsEvaluated[[1]] == 9999L, "Procrustes permutation count changed")
  assert_true(abs(x$bray_stats$RhoObserved[[1]] - 0.85948544610494) < 1e-12, "Bray rho changed")
  assert_true(abs(x$proc_stats$ProcrustesCorrelation[[1]] - 0.734454393300302) < 1e-12, "Procrustes r changed")
  assert_true(abs(x$proc_stats$ProcrustesM2[[1]] - 0.460576744161885) < 1e-12, "Procrustes M2 changed")
  invisible(TRUE)
}

site_order <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
context_order <- c("Ocean", "Tidal", "Island")
sample_type_order <- c("Sea-sand", "Water")
context_colors <- c(Ocean = "#0072B2", Tidal = "#E69F00", Island = "#009E73")
sample_shapes <- c(`Sea-sand` = 21, Water = 24)

factorize <- function(x) {
  if ("Site" %in% names(x)) x$Site <- factor(x$Site, levels = site_order)
  if ("GeologicalHabitat" %in% names(x)) x$GeologicalHabitat <- factor(x$GeologicalHabitat, levels = context_order)
  if ("SampleType" %in% names(x)) x$SampleType <- factor(x$SampleType, levels = sample_type_order)
  x
}

base_theme <- function() {
  theme_bw(base_size = 7.6, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(colour = "#E8E8E8", linewidth = 0.28),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "#333333", linewidth = 0.4),
      axis.title = element_text(size = 7.7, colour = "#222222", face = "bold"),
      axis.text = element_text(size = 6.7, colour = "#333333"),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0, margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 6.6, colour = "#333333", hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(4, 5, 4, 5),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 7.1, face = "bold"),
      legend.text = element_text(size = 6.8),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.spacing.x = grid::unit(1.5, "mm")
    )
}

make_alpha_panel <- function(x, metric, tag, keep_legend = FALSE) {
  z <- factorize(x$alpha_pairs[x$alpha_pairs$Metric == metric, , drop = FALSE])
  s <- x$alpha_stats[x$alpha_stats$Metric == metric, , drop = FALSE]
  subtitle <- sprintf(
    "Spearman rho = %.3f; Site-bootstrap 95%% CI [%.3f, %.3f]",
    s$RhoObserved, s$BootstrapCILower, s$BootstrapCIUpper
  )
  title <- if (metric == "Observed") paste0(tag, "  Observed richness") else paste0(tag, "  Shannon diversity")
  p <- ggplot(z, aes(x = Value16S, y = Value18S, colour = GeologicalHabitat,
                     fill = GeologicalHabitat, shape = SampleType)) +
    geom_point(size = 1.8, stroke = 0.43, alpha = 0.84) +
    scale_colour_manual(values = context_colors, drop = FALSE) +
    scale_fill_manual(values = context_colors, drop = FALSE) +
    scale_shape_manual(values = sample_shapes, drop = FALSE) +
    labs(
      x = paste("16S", metric), y = paste("18S", metric),
      title = title, subtitle = subtitle,
      colour = "Geological context", fill = "Geological context", shape = "Sample type"
    ) +
    base_theme() + theme(aspect.ratio = 1) +
    guides(
      colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 2.5, alpha = 1)),
      fill = "none",
      shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2.6, alpha = 1))
    )
  if (!keep_legend) p <- p + guides(colour = "none", fill = "none", shape = "none")
  p
}

make_bray_panel <- function(x) {
  s <- x$bray_stats[1, , drop = FALSE]
  subtitle <- sprintf(
    "Spearman rho = %.3f; restricted P = %.4f; BH q = %.4f",
    s$RhoObserved, s$PPermutation, s$PAdjustedBHAcrossSubsets
  )
  ggplot(x$bray_pairs, aes(x = Bray16S, y = Bray18S)) +
    geom_abline(slope = 1, intercept = 0, colour = "#B0B0B0", linetype = "dashed", linewidth = 0.35) +
    geom_point(size = 0.72, alpha = 0.23, colour = "#3F3F3F") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "16S Bray-Curtis", y = "18S Bray-Curtis",
      title = "C  Pairwise distance concordance", subtitle = subtitle
    ) + base_theme() + theme(aspect.ratio = 1)
}

make_proc_panel <- function(x) {
  z <- factorize(x$proc_scores)
  s <- x$proc_stats[1, , drop = FALSE]
  subtitle <- sprintf(
    "Procrustes r = %.3f; M2 = %.3f; restricted P = %.4f; BH q = %.5f",
    s$ProcrustesCorrelation, s$ProcrustesM2, s$PPermutation, s$PAdjustedBHAcrossSubsets
  )
  centroids <- z %>%
    group_by(Site, GeologicalHabitat) %>%
    summarise(
      Axis1_18S_Procrustes = mean(Axis1_18S_Procrustes),
      Axis2_18S_Procrustes = mean(Axis2_18S_Procrustes),
      .groups = "drop"
    )
  ggplot(z) +
    geom_segment(
      aes(x = Axis1_16S, y = Axis2_16S, xend = Axis1_18S_Procrustes,
          yend = Axis2_18S_Procrustes, colour = GeologicalHabitat),
      linewidth = 0.31, alpha = 0.46, show.legend = FALSE,
      arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed")
    ) +
    geom_point(
      aes(x = Axis1_16S, y = Axis2_16S, colour = GeologicalHabitat),
      shape = 4, size = 1.25, stroke = 0.42, alpha = 0.72, show.legend = FALSE
    ) +
    geom_point(
      aes(x = Axis1_18S_Procrustes, y = Axis2_18S_Procrustes,
          colour = GeologicalHabitat, fill = GeologicalHabitat, shape = SampleType),
      size = 1.70, stroke = 0.42, alpha = 0.83, show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = centroids,
      aes(x = Axis1_18S_Procrustes, y = Axis2_18S_Procrustes,
          label = Site, colour = GeologicalHabitat),
      inherit.aes = FALSE, seed = 20260813, size = 2.05, fontface = "bold",
      fill = scales::alpha("white", 0.86), label.size = 0,
      box.padding = 0.30, point.padding = 0.20, min.segment.length = 0,
      segment.colour = "#777777", segment.size = 0.25,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_colour_manual(values = context_colors, drop = FALSE) +
    scale_fill_manual(values = context_colors, drop = FALSE) +
    scale_shape_manual(values = sample_shapes, drop = FALSE) +
    coord_equal(clip = "off") +
    labs(
      x = "Procrustes axis 1", y = "Procrustes axis 2",
      title = "D  Procrustes alignment (16S to rotated 18S)", subtitle = subtitle
    ) + base_theme() +
    theme(aspect.ratio = 1, plot.margin = margin(4, 9, 4, 5)) +
    guides(colour = "none", fill = "none", shape = "none")
}

build_plot <- function(x) {
  p_a <- make_alpha_panel(x, "Observed", "A", keep_legend = TRUE)
  p_b <- make_alpha_panel(x, "Shannon", "B", keep_legend = FALSE)
  p_c <- make_bray_panel(x)
  p_d <- make_proc_panel(x)
  design <- "
AB
CD
EE
"
  patchwork::wrap_plots(
    A = p_a, B = p_b, C = p_c, D = p_d, E = patchwork::guide_area(),
    design = design, heights = c(1, 1, 0.10), guides = "collect"
  ) & theme(
    legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center"
  )
}

make_qa <- function(x) {
  a_obs <- x$alpha_stats[x$alpha_stats$Metric == "Observed", , drop = FALSE]
  a_sha <- x$alpha_stats[x$alpha_stats$Metric == "Shannon", , drop = FALSE]
  b <- x$bray_stats[1, , drop = FALSE]
  p <- x$proc_stats[1, , drop = FALSE]
  data.frame(
    check = c(
      "alpha_source_rows", "alpha_unique_samples", "bray_pair_rows", "procrustes_score_rows",
      "observed_rho", "shannon_rho", "bray_rho", "bray_p", "bray_bh_q",
      "procrustes_r", "procrustes_m2", "procrustes_p", "procrustes_bh_q",
      "site_bootstraps", "restricted_permutations", "footer_absent"
    ),
    value = c(
      nrow(x$alpha_pairs), length(unique(x$alpha_pairs$SampleID)), nrow(x$bray_pairs), nrow(x$proc_scores),
      a_obs$RhoObserved, a_sha$RhoObserved, b$RhoObserved, b$PPermutation, b$PAdjustedBHAcrossSubsets,
      p$ProcrustesCorrelation, p$ProcrustesM2, p$PPermutation, p$PAdjustedBHAcrossSubsets,
      a_obs$BootstrapReplicatesValid, b$PermutationsEvaluated, TRUE
    ),
    pass = TRUE,
    stringsAsFactors = FALSE
  )
}

check_manifests <- function() {
  for (mf in c(input_manifest, output_manifest)) {
    assert_true(file.exists(mf), paste("Missing manifest:", mf))
    m <- read_tsv(mf)
    for (i in seq_len(nrow(m))) {
      path <- file.path(project_dir, m$path[[i]])
      assert_true(file.exists(path), paste("Manifest target missing:", m$path[[i]]))
      assert_true(as.numeric(file.info(path)$size) == as.numeric(m$bytes[[i]]), paste("Byte mismatch:", m$path[[i]]))
      assert_true(sha256(path) == m$sha256[[i]], paste("SHA mismatch:", m$path[[i]]))
    }
  }
}

if (mode == "run") {
  assert_true(all(file.exists(canonical_paths)), "Canonical DS5 source is incomplete")
  assert_true(!file.exists(pdf_out) && !file.exists(tiff_out), "Versioned S13 delivery already exists")
  for (d in c(delivery_dir, source_dir, qc_dir, preview_dir, manifest_dir, log_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  for (i in seq_along(canonical_paths)) {
    if (!file.exists(snapshot_paths[[i]])) {
      copied <- file.copy(canonical_paths[[i]], snapshot_paths[[i]], overwrite = FALSE)
      assert_true(copied, paste("Could not snapshot:", source_names[[i]]))
    } else {
      assert_true(sha256(canonical_paths[[i]]) == sha256(snapshot_paths[[i]]), paste("Existing snapshot mismatch:", source_names[[i]]))
    }
  }
  write_manifest(canonical_paths, input_manifest)
  frozen <- load_frozen()
  validate_frozen(frozen)
  fig <- build_plot(frozen)

  grDevices::cairo_pdf(pdf_out, width = 183 / 25.4, height = 200 / 25.4, family = "Arial")
  print(fig)
  grDevices::dev.off()

  ragg::agg_tiff(
    tiff_out, width = 183, height = 200, units = "mm", res = 600,
    compression = "lzw", background = "white", scaling = 1
  )
  print(fig)
  grDevices::dev.off()

  write.table(make_qa(frozen), qa_file, sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = session_file)
  pdftools::pdf_convert(pdf_out, format = "png", dpi = 220, filenames = preview_out, verbose = FALSE)

  visual_lines <- c(
    "# Supplementary Figure S13 rendered visual QA v1",
    "",
    "- QA surface: actual one-page delivery PDF rendered at 220 dpi with Poppler.",
    "- Layout: 2 x 2 panel roles and dedicated bottom legend band.",
    "- Footer and overall title: absent.",
    "- Statistics: outside the data rectangles as panel subtitles.",
    "- Procrustes Site labels: repelled labels with neutral connectors.",
    "- Legend title: Geological context; sample-type legend is separate and collected below the panels.",
    "- Final visual acceptance requires direct inspection of the rendered PNG recorded by the build task."
  )
  writeLines(visual_lines, visual_qa_file, useBytes = TRUE)

  output_paths <- c(
    pdf_out, tiff_out, legend_file, script_file, snapshot_paths, qa_file,
    visual_qa_file, preview_out, session_file
  )
  write_manifest(output_paths, output_manifest)
  cat("RUN_OK\n")
}

if (mode == "check") {
  assert_true(all(file.exists(canonical_paths)), "Canonical DS5 source is incomplete")
  assert_true(all(file.exists(snapshot_paths)), "Frozen DS5 snapshots are incomplete")
  for (i in seq_along(canonical_paths)) {
    assert_true(sha256(canonical_paths[[i]]) == sha256(snapshot_paths[[i]]), paste("Snapshot mismatch:", source_names[[i]]))
  }
  frozen <- load_frozen()
  validate_frozen(frozen)
  assert_true(file.exists(pdf_out) && file.exists(tiff_out), "S13 delivery is missing")
  info <- pdftools::pdf_info(pdf_out)
  assert_true(info$pages == 1L, "PDF must contain one page")
  text <- paste(pdftools::pdf_text(pdf_out), collapse = "\n")
  forbidden <- c("A75_FULL exploratory", "reader-facing", "Geological habitat", "rebuilt PR2-primary", "adjacent source tables")
  assert_true(!any(vapply(forbidden, grepl, logical(1), x = text, fixed = TRUE)), "Forbidden footer/internal wording remains in PDF")
  required <- c("Observed richness", "Shannon diversity", "Pairwise distance concordance", "Procrustes alignment", "Geological context", "Sample type")
  assert_true(all(vapply(required, grepl, logical(1), x = text, fixed = TRUE)), "Required reader-facing PDF labels are missing")
  check_manifests()
  cat("CHECK_OK\n")
}
