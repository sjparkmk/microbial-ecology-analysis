#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(pdftools))

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
if (length(file_arg) != 1L) stop("Could not resolve script path", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script), "..", "..", "..", ".."),
                              winslash = "/", mustWork = TRUE)
package_root <- file.path(project_root, "revision", "re_analysis",
                          "63_figure_presentation_refinement_R_v1")
generated <- file.path(package_root, "generated")
out <- file.path(package_root, "qc", "rendered_pdf_pages")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pdfs <- c(
  file.path(generated, "figure1", "delivery", "Figure_1_16S_18S_taxonomic_composition_v1.pdf"),
  file.path(generated, "hierarchy", "figures", "Figure_3_16S_18S_alpha_effect_sizes_v1.pdf"),
  file.path(generated, "hierarchy", "figures", "Figure_2_16S_18S_beta_diversity_v1.pdf"),
  file.path(generated, "daa", "artifacts_v1", "figures", "Figure_4_paired_sample_type_differential_abundance_v1.pdf"),
  file.path(generated, "f5", "delivery", "Supplementary_Figure_S13_cross_marker_concordance_R_v1.pdf"),
  file.path(generated, "s2", "delivery", "Supplementary_Figure_S3_taxonomic_bubble_R_v2.pdf"),
  file.path(generated, "core_delivery", "Supplementary_Figure_S5_16S_core_ASV_patterns_R_v5.pdf"),
  file.path(generated, "core_delivery", "Supplementary_Figure_S6_18S_core_ASV_patterns_R_v5.pdf"),
  file.path(generated, "core_delivery", "Supplementary_Figure_S7_core_ASV_phylogeny_abundance_R_v5.pdf")
)
if (!all(file.exists(pdfs))) stop("One or more refinement PDFs are missing", call. = FALSE)

for (i in seq_along(pdfs)) {
  stem <- sprintf("%02d_%s", i, tools::file_path_sans_ext(basename(pdfs[[i]])))
  target <- file.path(out, paste0(stem, "_%d.png"))
  pdf_convert(pdfs[[i]], format = "png", dpi = 180, filenames = target, verbose = FALSE)
}

cat("RENDER_OK\n")
