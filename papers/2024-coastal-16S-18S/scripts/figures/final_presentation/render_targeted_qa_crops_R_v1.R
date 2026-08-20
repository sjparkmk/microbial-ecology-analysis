#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pdftools)
  library(png)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
if (length(file_arg) != 1L) stop("Could not resolve script path", call. = FALSE)
script <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script), "..", "..", "..", ".."),
                              winslash = "/", mustWork = TRUE)
package_root <- file.path(project_root, "revision", "re_analysis",
                          "63_figure_presentation_refinement_R_v1")
out <- file.path(package_root, "qc", "targeted_crops")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

render_and_crop <- function(pdf, full_name, crop_name, x1, x2, y1, y2) {
  full <- file.path(out, full_name)
  pdf_convert(pdf, format = "png", dpi = 600, filenames = full, verbose = FALSE)
  image <- readPNG(full)
  stopifnot(y1 >= 1L, x1 >= 1L, y2 <= dim(image)[1L], x2 <= dim(image)[2L])
  writePNG(image[y1:y2, x1:x2, , drop = FALSE], file.path(out, crop_name))
}

render_and_crop(
  file.path(package_root, "generated", "s2", "delivery",
            "Supplementary_Figure_S3_taxonomic_bubble_R_v2.pdf"),
  "S2_full_600dpi.png", "S2_Ulleungdo_top_crop_600dpi.png",
  1750L, 2500L, 1L, 550L
)

render_and_crop(
  file.path(package_root, "generated", "daa", "artifacts_v1", "figures",
            "Figure_4_paired_sample_type_differential_abundance_v1.pdf"),
  "F4_full_600dpi.png", "F4_Zeaxanthinibacter_crop_600dpi.png",
  300L, 1900L, 220L, 1050L
)

cat("TARGETED_CROPS_OK\n")
