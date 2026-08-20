#!/usr/bin/env Rscript

# Build the manuscript-delivery Figure 1 from accepted C05 composition cell means.
# Numerical composition values are not recalculated here. The script validates
# the C05 sources against the parallel DS3 cell-mean sources, then performs only
# presentation-level recomposition into a two-panel PDF/TIFF figure.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

SCRIPT_VERSION <- "1.0.0"
FIGURE_STEM <- "Figure_1_16S_18S_taxonomic_composition_v1"
PACKAGE_REL <- file.path("revision", "re_analysis", "63_figure_presentation_refinement_R_v1", "generated", "figure1", "package")
DELIVERY_REL <- file.path("revision", "re_analysis", "63_figure_presentation_refinement_R_v1", "generated", "figure1", "delivery")
WIDTH_MM <- 174
HEIGHT_MM <- 240
DPI <- 600

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
  invisible(TRUE)
}

parse_cli <- function(args) {
  out <- list(mode = NULL, project_root = getwd())
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--run", "--dry-run", "--check")) {
      assert_true(is.null(out$mode), "Specify exactly one mode")
      out$mode <- substring(arg, 3L)
    } else if (arg == "--project-root") {
      assert_true(i < length(args), "--project-root requires a value")
      i <- i + 1L
      out$project_root <- args[[i]]
    } else if (arg == "--help") {
      cat(
        "Usage: Rscript build_revision_figure1_v1.R --run|--dry-run|--check ",
        "[--project-root PATH]\n",
        sep = ""
      )
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }
  assert_true(!is.null(out$mode), "One of --run, --dry-run, or --check is required")
  out
}

read_tsv <- function(path) {
  assert_true(file.exists(path), paste0("Missing input: ", path))
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  assert_true(!file.exists(path), paste0("Refusing overwrite: ", path))
  write.table(x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "NA", fileEncoding = "UTF-8")
}

sha256_file <- function(path) {
  assert_true(requireNamespace("digest", quietly = TRUE), "R package 'digest' is required")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

file_record <- function(role, path, project_root) {
  info <- file.info(path)
  rel <- substring(normalizePath(path, winslash = "/", mustWork = TRUE),
                   nchar(normalizePath(project_root, winslash = "/", mustWork = TRUE)) + 2L)
  data.frame(
    Role = role,
    PathRelativeToProject = rel,
    Bytes = as.numeric(info$size),
    SHA256 = sha256_file(path),
    stringsAsFactors = FALSE
  )
}

taxon_palette <- function(labels) {
  labels <- unique(as.character(labels))
  fixed <- c(`Other assigned` = "#BDBDBD", Unassigned = "#666666")
  main <- setdiff(labels, names(fixed))
  hues <- grDevices::hcl(
    h = seq(15, 375, length.out = length(main) + 1L)[seq_along(main)],
    c = 55, l = 70
  )
  out <- c(stats::setNames(hues, main), fixed[intersect(names(fixed), labels)])
  if ("Other" %in% labels) out[["Other"]] <- "#BDBDBD"
  if ("Unassigned" %in% labels) out[["Unassigned"]] <- "#666666"
  out[labels]
}

validate_c05 <- function(x, marker, dataset, rank) {
  required <- c(
    "Marker", "Dataset", "Rank", "Site", "SiteCodeOriginal",
    "GeologicalHabitat", "SampleType", "Season", "TaxonDisplay",
    "CategoryOrder", "MeanRelativeAbundancePct", "SDRelativeAbundancePct", "N"
  )
  assert_true(identical(names(x), required), paste0(marker, " C05 schema mismatch"))
  assert_true(nrow(x) == 350L, paste0(marker, " C05 must contain 350 rows"))
  assert_true(identical(unique(x$Marker), marker), paste0(marker, " marker mismatch"))
  assert_true(identical(unique(x$Dataset), dataset), paste0(marker, " dataset mismatch"))
  assert_true(identical(unique(x$Rank), rank), paste0(marker, " rank mismatch"))
  x$CategoryOrder <- as.integer(x$CategoryOrder)
  x$MeanRelativeAbundancePct <- as.numeric(x$MeanRelativeAbundancePct)
  x$SDRelativeAbundancePct <- as.numeric(x$SDRelativeAbundancePct)
  x$N <- as.integer(x$N)
  assert_true(all(is.finite(x$MeanRelativeAbundancePct)) &&
                all(x$MeanRelativeAbundancePct >= 0),
              paste0(marker, " invalid mean relative abundance"))
  assert_true(all(x$N == 3L), paste0(marker, " every available cell must have N=3"))
  assert_true(setequal(unique(x$CategoryOrder), 1:14),
              paste0(marker, " must contain category orders 1:14"))
  taxon_map <- unique(x[c("TaxonDisplay", "CategoryOrder")])
  assert_true(nrow(taxon_map) == 14L && !anyDuplicated(taxon_map$TaxonDisplay) &&
                !anyDuplicated(taxon_map$CategoryOrder),
              paste0(marker, " category mapping is not one-to-one"))
  taxon_map <- taxon_map[order(taxon_map$CategoryOrder), , drop = FALSE]
  assert_true(identical(tail(taxon_map$TaxonDisplay, 2L), c("Other", "Unassigned")),
              paste0(marker, " orders 13/14 must be Other/Unassigned"))
  key <- interaction(x$Site, x$Season, x$SampleType, drop = TRUE, lex.order = TRUE)
  assert_true(length(unique(key)) == 25L, paste0(marker, " must contain 25 available cells"))
  assert_true(all(as.integer(table(key)) == 14L),
              paste0(marker, " every available cell must contain 14 categories"))
  sums <- tapply(x$MeanRelativeAbundancePct, key, sum)
  assert_true(max(abs(sums - 100)) < 1e-8,
              paste0(marker, " cell composition sums are not 100%"))
  x
}

validate_missing_cells <- function(x16, x18) {
  sites <- c("Dokdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan", "Ulleungdo")
  expected <- expand.grid(
    Site = sites,
    Season = c("Spring", "Autumn"),
    SampleType = c("Sea-sand", "Water"),
    stringsAsFactors = FALSE
  )
  expected$key <- paste(expected$Site, expected$Season, expected$SampleType, sep = "|")
  keys <- function(x) unique(paste(x$Site, x$Season, x$SampleType, sep = "|"))
  missing16 <- sort(setdiff(expected$key, keys(x16)))
  missing18 <- sort(setdiff(expected$key, keys(x18)))
  wanted <- sort(c(
    "Dokdo|Autumn|Water",
    "Ulleungdo|Autumn|Sea-sand",
    "Ulleungdo|Autumn|Water"
  ))
  assert_true(identical(missing16, wanted), "16S structural missing-cell set mismatch")
  assert_true(identical(missing18, wanted), "18S structural missing-cell set mismatch")
  wanted
}

compare_ds3 <- function(c05, ds3, marker) {
  required <- c("Marker", "Dataset", "Rank", "CellID", "Site", "SiteCodeOriginal",
                "GeologicalHabitat", "SampleType", "Campaign", "TaxonDisplay",
                "CategoryOrder", "MeanRelativeAbundance", "SDRelativeAbundance", "N",
                "CampaignSource", "Season")
  assert_true(identical(names(ds3), required), paste0(marker, " DS3 schema mismatch"))
  key_cols_c05 <- c("Site", "SiteCodeOriginal", "GeologicalHabitat", "SampleType",
                    "Season", "TaxonDisplay", "CategoryOrder")
  key_cols_ds3 <- c("Site", "SiteCodeOriginal", "GeologicalHabitat", "SampleType",
                    "Season", "TaxonDisplay", "CategoryOrder")
  a <- c05[order(do.call(order, c05[key_cols_c05])), , drop = FALSE]
  b <- ds3[order(do.call(order, ds3[key_cols_ds3])), , drop = FALSE]
  assert_true(nrow(a) == nrow(b), paste0(marker, " C05/DS3 row-count mismatch"))
  ak <- do.call(paste, c(a[key_cols_c05], sep = "|"))
  bk <- do.call(paste, c(b[key_cols_ds3], sep = "|"))
  assert_true(identical(ak, bk), paste0(marker, " C05/DS3 key mismatch"))
  max_mean <- max(abs(a$MeanRelativeAbundancePct / 100 - as.numeric(b$MeanRelativeAbundance)))
  max_sd <- max(abs(a$SDRelativeAbundancePct / 100 - as.numeric(b$SDRelativeAbundance)))
  assert_true(max_mean < 2e-12 && max_sd < 2e-12,
              paste0(marker, " C05/DS3 numeric mismatch"))
  c(max_mean = max_mean, max_sd = max_sd)
}

make_panel <- function(x, marker, rank, panel_tag) {
  site_levels <- c("Dokdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan", "Ulleungdo")
  taxon_map <- unique(x[c("TaxonDisplay", "CategoryOrder")])
  taxon_map <- taxon_map[order(taxon_map$CategoryOrder), , drop = FALSE]
  taxon_levels <- taxon_map$TaxonDisplay
  x$Site <- factor(x$Site, levels = site_levels)
  x$Season <- factor(x$Season, levels = c("Spring", "Autumn"))
  x$SampleType <- factor(x$SampleType, levels = c("Sea-sand", "Water"))
  x$TaxonDisplay <- factor(x$TaxonDisplay, levels = taxon_levels)
  cols <- taxon_palette(taxon_levels)
  ggplot(x, aes(x = Site, y = MeanRelativeAbundancePct, fill = TaxonDisplay)) +
    geom_col(width = 0.80, colour = "white", linewidth = 0.08) +
    facet_grid(Season ~ SampleType, drop = FALSE) +
    scale_fill_manual(values = cols, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, 100.01), expand = c(0, 0),
      breaks = seq(0, 100, 25), labels = function(v) paste0(v, "%")
    ) +
    labs(
      tag = panel_tag,
      x = NULL,
      y = "Mean relative abundance",
      fill = paste0(marker, " ", tolower(rank))
    ) +
    theme_classic(base_size = 7.5, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(colour = "black", face = "bold"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(2.55, "mm"),
      legend.key.width = grid::unit(3.1, "mm"),
      legend.text = element_text(size = 6.7),
      panel.spacing = grid::unit(2.0, "mm"),
      plot.tag = element_text(face = "bold", size = 10, colour = "black"),
      plot.tag.position = c(0.002, 1.015),
      plot.margin = margin(4, 4, 2, 4, unit = "mm")
    ) +
    guides(fill = guide_legend(ncol = 1, byrow = FALSE))
}

write_lines_new <- function(lines, path) {
  assert_true(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

manifest_for_tree <- function(root, project_root, reported_root = root,
                              role_prefix = "package_output") {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[!dir.exists(files)]
  files <- files[basename(files) != "Figure_1_output_sha256_manifest_v1.tsv"]
  if (!length(files)) return(data.frame())
  do.call(rbind, lapply(files, function(f) {
    rec <- file_record(role_prefix, f, project_root)
    reported <- file.path(reported_root, substring(
      normalizePath(f, winslash = "/", mustWork = TRUE),
      nchar(normalizePath(root, winslash = "/", mustWork = TRUE)) + 2L
    ))
    rec$PathRelativeToProject <- substring(
      normalizePath(reported, winslash = "/", mustWork = FALSE),
      nchar(normalizePath(project_root, winslash = "/", mustWork = TRUE)) + 2L
    )
    rec
  }))
}

validate_signatures <- function(pdf, tiff) {
  pdf_head <- readBin(pdf, what = "raw", n = 5L)
  assert_true(identical(rawToChar(pdf_head), "%PDF-"), "PDF signature mismatch")
  tif_head <- readBin(tiff, what = "raw", n = 4L)
  valid_tif <- identical(as.integer(tif_head), c(73L, 73L, 42L, 0L)) ||
    identical(as.integer(tif_head), c(77L, 77L, 0L, 42L))
  assert_true(valid_tif, "TIFF signature mismatch")
  if (requireNamespace("pdftools", quietly = TRUE)) {
    info <- pdftools::pdf_info(pdf)
    assert_true(info$pages == 1L, "Figure PDF must contain one page")
  }
  invisible(TRUE)
}

load_and_validate <- function(project_root) {
  p <- function(...) file.path(project_root, ...)
  inputs <- list(
    c05_16s = p("revision", "re_analysis", "12_original_pipeline_reanalysis_v1", "core_v1",
                "source_data", "core_C05_16s_composition_cell_source_v1.tsv"),
    c05_18s = p("revision", "re_analysis", "12_original_pipeline_reanalysis_v1", "core_v1",
                "source_data", "core_C05_18s_composition_cell_source_v1.tsv"),
    ds3_16s = p("revision", "re_analysis", "10_figure_source_data_v4", "source_data",
                "DS3a_16S_phylum_cell_mean_source_data_v1.tsv"),
    ds3_18s = p("revision", "re_analysis", "10_figure_source_data_v4", "source_data",
                "DS3b_18S_PR2_class_cell_mean_source_data_v1.tsv"),
    comp_readme = p("revision", "re_analysis", "07_taxonomic_composition_v1", "provenance",
                    "README_taxonomic_composition_v1.md"),
    c05_script = p("revision", "re_analysis", "original_style_pipeline_v1",
                   "01_core_original_figures.R"),
    builder_script = p("revision", "re_analysis", "scripts", "build_revision_figure1_v1.R")
  )
  missing <- names(inputs)[!vapply(inputs, file.exists, logical(1))]
  assert_true(!length(missing), paste0("Missing input(s): ", paste(missing, collapse = ", ")))
  x16 <- validate_c05(read_tsv(inputs$c05_16s), "16S", "primary_rarefied", "Phylum")
  x18 <- validate_c05(read_tsv(inputs$c05_18s), "18S", "full_pr2_primary", "Class")
  missing_cells <- validate_missing_cells(x16, x18)
  d16 <- compare_ds3(x16, read_tsv(inputs$ds3_16s), "16S")
  d18 <- compare_ds3(x18, read_tsv(inputs$ds3_18s), "18S")
  list(inputs = inputs, x16 = x16, x18 = x18, missing_cells = missing_cells,
       ds3_delta_16s = d16, ds3_delta_18s = d18)
}

run_build <- function(project_root, validated) {
  package_final <- file.path(project_root, PACKAGE_REL)
  delivery_dir <- file.path(project_root, DELIVERY_REL)
  pdf_final <- file.path(delivery_dir, paste0(FIGURE_STEM, ".pdf"))
  tiff_final <- file.path(delivery_dir, paste0(FIGURE_STEM, ".tiff"))
  assert_true(!file.exists(package_final), paste0("Refusing existing package: ", package_final))
  assert_true(!file.exists(pdf_final) && !file.exists(tiff_final), "Refusing existing delivery figure")

  stage <- paste0(package_final, ".staging-", Sys.getpid())
  assert_true(!file.exists(stage), paste0("Staging path exists: ", stage))
  dir.create(file.path(stage, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "source_data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "provenance"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "qc"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(stage, "previews"), recursive = TRUE, showWarnings = FALSE)
  promoted_package <- FALSE
  promoted_delivery <- character()
  on.exit({
    if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
    if (!promoted_package && length(promoted_delivery)) unlink(promoted_delivery, force = TRUE)
  }, add = TRUE)

  src16_out <- file.path(stage, "source_data", "Figure_1_16S_phylum_cell_means_v1.tsv")
  src18_out <- file.path(stage, "source_data", "Figure_1_18S_PR2_class_cell_means_v1.tsv")
  file.copy(validated$inputs$c05_16s, src16_out, copy.mode = TRUE)
  file.copy(validated$inputs$c05_18s, src18_out, copy.mode = TRUE)
  assert_true(identical(sha256_file(src16_out), sha256_file(validated$inputs$c05_16s)),
              "16S source copy mismatch")
  assert_true(identical(sha256_file(src18_out), sha256_file(validated$inputs$c05_18s)),
              "18S source copy mismatch")

  palette_rows <- bind_rows(lapply(list(`16S` = validated$x16, `18S` = validated$x18), function(x) {
    m <- unique(x[c("Marker", "Rank", "TaxonDisplay", "CategoryOrder")])
    m <- m[order(m$CategoryOrder), , drop = FALSE]
    m$HexColor <- unname(taxon_palette(m$TaxonDisplay))
    m
  }))
  write_tsv(palette_rows, file.path(stage, "source_data", "Figure_1_taxon_palette_v1.tsv"))

  p16 <- make_panel(validated$x16, "16S", "Phylum", "A")
  p18 <- make_panel(validated$x18, "18S", "Class", "B")
  figure <- p16 / p18 + plot_layout(heights = c(1, 1))
  pdf_stage <- file.path(stage, "figures", paste0(FIGURE_STEM, ".pdf"))
  tiff_stage <- file.path(stage, "figures", paste0(FIGURE_STEM, ".tiff"))
  preview_stage <- file.path(stage, "previews", paste0(FIGURE_STEM, "_preview.png"))
  ggsave(pdf_stage, figure, width = WIDTH_MM, height = HEIGHT_MM, units = "mm",
         device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)
  ggsave(tiff_stage, figure, width = WIDTH_MM, height = HEIGHT_MM, units = "mm",
         device = "tiff", dpi = DPI, compression = "lzw", bg = "white", limitsize = FALSE)
  ggsave(preview_stage, figure, width = WIDTH_MM, height = HEIGHT_MM, units = "mm",
         device = "png", dpi = 200, bg = "white", limitsize = FALSE)
  validate_signatures(pdf_stage, tiff_stage)

  contract <- c(
    "# Figure 1 contract v1",
    "",
    "- Status: manuscript-delivery candidate pending source and rendered visual QA.",
    "- Core conclusion: dominant 16S phylum and full PR2-primary 18S class composition varies descriptively among Sites, Sample types, and Seasons.",
    "- Panel map: (a) 16S Phylum; (b) 18S PR2-primary Class.",
    "- Backend: R/ggplot2 + patchwork only.",
    paste0("- Final size: ", WIDTH_MM, " x ", HEIGHT_MM, " mm; one-page editable PDF and ", DPI, "-dpi LZW TIFF."),
    "- Source: accepted C05 Site x Sample type x Season cell means; DS3 cell means are an independent numeric cross-check.",
    "- Aggregation: sample-wise relative abundance first, followed by the arithmetic mean of three spatial biological replicates within each available cell.",
    "- Season and Sample type are displayed separately and are not averaged together.",
    "- Categories: globally selected top 12 assigned taxa per marker/rank, all remaining assigned taxa as Other, and unresolved rank assignments as Unassigned.",
    "- Missing cells remain blank: Dokdo-Autumn-Water, Ulleungdo-Autumn-Sea-sand, and Ulleungdo-Autumn-Water.",
    "- Evidence hierarchy: descriptive composition only; the figure is not differential-abundance evidence.",
    "- Image integrity: stacked bars represent summary values only; no images or microscopy panels are present.",
    "- Embedded text restriction: panel tags, facet labels, axis labels, tick labels, Site labels, and legend labels only; no overall title or caption is embedded.",
    "- Original-pipeline alignment: retains the stacked-bar visual grammar referenced by tax_stakplot_final_ordered_soft.R while using canonical top-12 categories and design-aware cell means."
  )
  write_lines_new(contract, file.path(stage, "provenance", "Figure_1_contract_v1.md"))

  caption_facts <- c(
    "# Figure 1 caption facts v1",
    "",
    "Panel (a) shows 16S phylum composition from the common 75-sample dataset rarefied to 20,000 reads per sample. Panel (b) shows full PR2 v5.1.1-primary 18S class composition from the same 75 samples rarefied to 50,000 reads per sample.",
    "",
    "For each available Site x Season x Sample type cell, each stacked bar is the arithmetic mean of sample-wise relative abundances from three spatial biological field replicates. Reads were not pooled before relative-abundance calculation, and neither Season nor Sample type was averaged over.",
    "",
    "The 12 assigned taxa with the largest global mean sample-wise relative abundance are shown individually. Remaining assigned taxa are grouped as Other, whereas unresolved assignments at the displayed rank are retained as Unassigned.",
    "",
    "The common design contains 25 available cells. Dokdo Autumn Water and both Ulleungdo Autumn cells are absent and remain blank; no values were imputed. The composition patterns are descriptive and should not be interpreted as differential-abundance tests."
  )
  write_lines_new(caption_facts, file.path(stage, "provenance", "Figure_1_caption_facts_v1.md"))

  input_manifest <- do.call(rbind, Map(
    function(role, path) file_record(role, path, project_root),
    names(validated$inputs), unname(validated$inputs)
  ))
  write_tsv(input_manifest, file.path(stage, "provenance", "Figure_1_input_sha256_manifest_v1.tsv"))

  qa <- list(
    status = "PASS_NUMERIC_AND_FILE_QA",
    script_version = SCRIPT_VERSION,
    figure_stem = FIGURE_STEM,
    common_samples_per_marker = 75L,
    available_cells_per_marker = 25L,
    replicates_per_available_cell = 3L,
    displayed_categories_per_cell = 14L,
    top_assigned_taxa = 12L,
    rarefaction_depth_16s = 20000L,
    rarefaction_depth_18s = 50000L,
    marker_rank_16s = "Phylum",
    marker_rank_18s = "Class",
    missing_cells = validated$missing_cells,
    ds3_max_abs_delta_16s_mean = unname(validated$ds3_delta_16s[["max_mean"]]),
    ds3_max_abs_delta_18s_mean = unname(validated$ds3_delta_18s[["max_mean"]]),
    embedded_title = FALSE,
    embedded_caption = FALSE,
    output_width_mm = WIDTH_MM,
    output_height_mm = HEIGHT_MM,
    tiff_dpi_requested = DPI
  )
  assert_true(requireNamespace("jsonlite", quietly = TRUE), "R package 'jsonlite' is required")
  jsonlite::write_json(
    qa, file.path(stage, "qc", "Figure_1_validation_v1.json"),
    pretty = TRUE, auto_unbox = TRUE
  )

  package_manifest <- manifest_for_tree(stage, project_root, package_final)
  write_tsv(package_manifest, file.path(stage, "provenance", "Figure_1_output_sha256_manifest_v1.tsv"))
  assert_true(dir.create(dirname(package_final), recursive = TRUE, showWarnings = FALSE) ||
                dir.exists(dirname(package_final)), "Cannot create package parent")
  assert_true(file.rename(stage, package_final), "Atomic package promotion failed")
  promoted_package <- TRUE

  dir.create(delivery_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_tmp <- paste0(pdf_final, ".tmp")
  tiff_tmp <- paste0(tiff_final, ".tmp")
  on.exit(unlink(c(pdf_tmp, tiff_tmp), force = TRUE), add = TRUE)
  assert_true(file.copy(file.path(package_final, "figures", basename(pdf_final)), pdf_tmp,
                        copy.mode = TRUE), "PDF delivery staging copy failed")
  assert_true(file.copy(file.path(package_final, "figures", basename(tiff_final)), tiff_tmp,
                        copy.mode = TRUE), "TIFF delivery staging copy failed")
  assert_true(file.rename(pdf_tmp, pdf_final), "PDF delivery promotion failed")
  promoted_delivery <- c(promoted_delivery, pdf_final)
  assert_true(file.rename(tiff_tmp, tiff_final), "TIFF delivery promotion failed")
  promoted_delivery <- c(promoted_delivery, tiff_final)
  assert_true(identical(sha256_file(pdf_final), sha256_file(file.path(package_final, "figures", basename(pdf_final)))),
              "Delivered PDF hash mismatch")
  assert_true(identical(sha256_file(tiff_final), sha256_file(file.path(package_final, "figures", basename(tiff_final)))),
              "Delivered TIFF hash mismatch")
  validate_signatures(pdf_final, tiff_final)
  cat("RUN_OK\n")
  cat("PACKAGE=", normalizePath(package_final, winslash = "/"), "\n", sep = "")
  cat("PDF=", normalizePath(pdf_final, winslash = "/"), "\n", sep = "")
  cat("TIFF=", normalizePath(tiff_final, winslash = "/"), "\n", sep = "")
}

check_build <- function(project_root, validated) {
  package_final <- file.path(project_root, PACKAGE_REL)
  delivery_dir <- file.path(project_root, DELIVERY_REL)
  pdf_final <- file.path(delivery_dir, paste0(FIGURE_STEM, ".pdf"))
  tiff_final <- file.path(delivery_dir, paste0(FIGURE_STEM, ".tiff"))
  required <- c(
    package_final, pdf_final, tiff_final,
    file.path(package_final, "source_data", "Figure_1_16S_phylum_cell_means_v1.tsv"),
    file.path(package_final, "source_data", "Figure_1_18S_PR2_class_cell_means_v1.tsv"),
    file.path(package_final, "source_data", "Figure_1_taxon_palette_v1.tsv"),
    file.path(package_final, "provenance", "Figure_1_contract_v1.md"),
    file.path(package_final, "provenance", "Figure_1_caption_facts_v1.md"),
    file.path(package_final, "provenance", "Figure_1_input_sha256_manifest_v1.tsv"),
    file.path(package_final, "provenance", "Figure_1_output_sha256_manifest_v1.tsv"),
    file.path(package_final, "qc", "Figure_1_validation_v1.json"),
    file.path(package_final, "previews", paste0(FIGURE_STEM, "_preview.png"))
  )
  assert_true(all(file.exists(required)), "Required Figure 1 package file is missing")
  assert_true(identical(sha256_file(pdf_final),
                        sha256_file(file.path(package_final, "figures", basename(pdf_final)))),
              "Delivered/package PDF mismatch")
  assert_true(identical(sha256_file(tiff_final),
                        sha256_file(file.path(package_final, "figures", basename(tiff_final)))),
              "Delivered/package TIFF mismatch")
  assert_true(identical(sha256_file(file.path(package_final, "source_data", "Figure_1_16S_phylum_cell_means_v1.tsv")),
                        sha256_file(validated$inputs$c05_16s)), "16S source drift")
  assert_true(identical(sha256_file(file.path(package_final, "source_data", "Figure_1_18S_PR2_class_cell_means_v1.tsv")),
                        sha256_file(validated$inputs$c05_18s)), "18S source drift")

  input_manifest <- read_tsv(file.path(package_final, "provenance", "Figure_1_input_sha256_manifest_v1.tsv"))
  for (i in seq_len(nrow(input_manifest))) {
    path <- file.path(project_root, input_manifest$PathRelativeToProject[[i]])
    assert_true(file.exists(path), paste0("Manifest input missing: ", path))
    assert_true(as.numeric(file.info(path)$size) == as.numeric(input_manifest$Bytes[[i]]),
                paste0("Manifest input size mismatch: ", path))
    assert_true(identical(sha256_file(path), input_manifest$SHA256[[i]]),
                paste0("Manifest input hash mismatch: ", path))
  }
  output_manifest <- read_tsv(file.path(package_final, "provenance", "Figure_1_output_sha256_manifest_v1.tsv"))
  for (i in seq_len(nrow(output_manifest))) {
    path <- file.path(project_root, output_manifest$PathRelativeToProject[[i]])
    assert_true(file.exists(path), paste0("Manifest output missing: ", path))
    assert_true(as.numeric(file.info(path)$size) == as.numeric(output_manifest$Bytes[[i]]),
                paste0("Manifest output size mismatch: ", path))
    assert_true(identical(sha256_file(path), output_manifest$SHA256[[i]]),
                paste0("Manifest output hash mismatch: ", path))
  }
  validate_signatures(pdf_final, tiff_final)
  cat("CHECK_OK\n")
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(cli$project_root, winslash = "/", mustWork = TRUE)
validated <- load_and_validate(project_root)

if (cli$mode == "dry-run") {
  cat("DRY_RUN_OK\n")
  cat("16S_ROWS=", nrow(validated$x16), "\n", sep = "")
  cat("18S_ROWS=", nrow(validated$x18), "\n", sep = "")
  cat("AVAILABLE_CELLS_PER_MARKER=25\n")
  cat("MISSING_CELLS=", paste(validated$missing_cells, collapse = ";"), "\n", sep = "")
} else if (cli$mode == "run") {
  run_build(project_root, validated)
} else {
  check_build(project_root, validated)
}
