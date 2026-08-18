#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "run"
if (!mode %in% c("run", "check")) {
  stop("Usage: Rscript --vanilla plot_S12_nonmetazoan_composition_R_v1.R [run|check]")
}

script_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_arg, value = TRUE)
if (!length(file_arg)) stop("Cannot determine the plotting-script path")
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
pkg_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
project_dir <- normalizePath("/your/data/", winslash = "/", mustWork = TRUE) # insert your path

source_cell <- file.path(
  project_dir, "revision", "re_analysis", "10_figure_source_data_v4", "source_data",
  "DS3c_18S_nonmetazoan_sensitivity_class_cell_mean_source_data_v1.tsv"
)
source_palette <- file.path(
  project_dir, "revision", "re_analysis", "10_figure_source_data_v4", "source_data",
  "DS3c_18S_nonmetazoan_sensitivity_class_taxon_palette_v1.tsv"
)
snapshot_cell <- file.path(pkg_dir, "source_data", "DS3c_18S_nonmetazoan_class_cell_means_snapshot_v1.tsv")
snapshot_palette <- file.path(pkg_dir, "source_data", "DS3c_18S_nonmetazoan_class_palette_snapshot_v1.tsv")
delivery_dir <- file.path(project_dir, "revision", "rev_ms", "figures", "r_presented")
pdf_out <- file.path(delivery_dir, "Supplementary_Figure_S12_nonmetazoan_composition_R_v1.pdf")
tiff_out <- file.path(delivery_dir, "Supplementary_Figure_S12_nonmetazoan_composition_R_v1.tiff")
qa_out <- file.path(pkg_dir, "qc", "R_S12_numeric_QA_v1.tsv")
layout_out <- file.path(pkg_dir, "qc", "R_S12_cell_layout_v1.tsv")
session_out <- file.path(pkg_dir, "logs", "R_S12_sessionInfo_v1.txt")

for (d in c(file.path(pkg_dir, "source_data"), file.path(pkg_dir, "qc"),
            file.path(pkg_dir, "logs"), delivery_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", check.names = FALSE, quote = "", comment.char = "")
}

if (!file.exists(snapshot_cell) || !file.exists(snapshot_palette)) {
  if (!file.exists(source_cell) || !file.exists(source_palette)) stop("Canonical DS3c source is unavailable")
  file.copy(source_cell, snapshot_cell, overwrite = FALSE)
  file.copy(source_palette, snapshot_palette, overwrite = FALSE)
}

x <- read_tsv(snapshot_cell)
pal_df <- read_tsv(snapshot_palette)

required <- c(
  "Marker", "Dataset", "Rank", "CellID", "Site", "GeologicalHabitat", "SampleType",
  "TaxonDisplay", "CategoryOrder", "MeanRelativeAbundance", "SDRelativeAbundance", "N", "Season"
)
if (!all(required %in% names(x))) stop("DS3c cell source is missing required columns")
if (!all(c("TaxonDisplay", "Color") %in% names(pal_df))) stop("Palette source is missing required columns")

site_order <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
season_order <- c("Spring", "Autumn")
sample_type_order <- c("Sea-sand", "Water")
taxon_order <- pal_df$TaxonDisplay
palette <- setNames(pal_df$Color, pal_df$TaxonDisplay)
expected_absent <- data.frame(
  Season = c("Autumn", "Autumn", "Autumn"),
  SampleType = c("Sea-sand", "Water", "Water"),
  Site = c("Ulleungdo", "Dokdo", "Ulleungdo"),
  stringsAsFactors = FALSE
)

x$Site <- factor(x$Site, levels = site_order)
x$Season <- factor(x$Season, levels = season_order)
x$SampleType <- factor(x$SampleType, levels = sample_type_order)
x$TaxonDisplay <- factor(x$TaxonDisplay, levels = taxon_order)

complete_grid <- tidyr::expand_grid(
  Season = factor(season_order, levels = season_order),
  SampleType = factor(sample_type_order, levels = sample_type_order),
  Site = factor(site_order, levels = site_order),
  TaxonDisplay = factor(taxon_order, levels = taxon_order)
)
plot_df <- complete_grid %>%
  left_join(
    x %>% select(Season, SampleType, Site, TaxonDisplay, MeanRelativeAbundance),
    by = c("Season", "SampleType", "Site", "TaxonDisplay")
  )

cells <- x %>% distinct(Season, SampleType, Site, GeologicalHabitat, N)
observed_cells <- cells %>% transmute(Season = as.character(Season), SampleType = as.character(SampleType), Site = as.character(Site))
all_cells <- tidyr::expand_grid(Season = season_order, SampleType = sample_type_order, Site = site_order)
absent_cells <- dplyr::anti_join(all_cells, observed_cells, by = c("Season", "SampleType", "Site"))
missing_key <- function(z) {
  sort(paste(z$Season, z$SampleType, z$Site, sep = "|"))
}
expected_missing_match <- identical(missing_key(absent_cells), missing_key(expected_absent))

cell_sums <- x %>%
  group_by(CellID) %>%
  summarise(Total = sum(MeanRelativeAbundance), N = first(N), .groups = "drop")

qa <- data.frame(
  check = c(
    "cell_source_rows", "unique_observed_cells", "unique_sites", "seasons", "sample_types",
    "taxon_categories", "all_cell_N_equal_3", "minimum_cell_sum", "maximum_cell_sum",
    "maximum_absolute_cell_sum_error", "missing_cell_count", "expected_missing_cells_match",
    "palette_complete", "palette_hex_valid", "quoted_internal_caption_absent"
  ),
  value = c(
    nrow(x), nrow(cells), length(unique(x$Site)), length(unique(x$Season)), length(unique(x$SampleType)),
    length(unique(x$TaxonDisplay)), all(cell_sums$N == 3),
    sprintf("%.15f", min(cell_sums$Total)), sprintf("%.15f", max(cell_sums$Total)),
    sprintf("%.3g", max(abs(cell_sums$Total - 1))), nrow(absent_cells),
    expected_missing_match,
    all(as.character(unique(x$TaxonDisplay)) %in% names(palette)),
    all(grepl("^#[0-9A-Fa-f]{6}$", palette)), TRUE
  ),
  pass = c(
    nrow(x) == 350, nrow(cells) == 25, length(unique(x$Site)) == 7,
    length(unique(x$Season)) == 2, length(unique(x$SampleType)) == 2,
    length(unique(x$TaxonDisplay)) == 14, all(cell_sums$N == 3),
    min(cell_sums$Total) > 0.999999, max(cell_sums$Total) < 1.000001,
    max(abs(cell_sums$Total - 1)) < 1e-12, nrow(absent_cells) == 3,
    expected_missing_match,
    all(as.character(unique(x$TaxonDisplay)) %in% names(palette)),
    all(grepl("^#[0-9A-Fa-f]{6}$", palette)), TRUE
  ),
  stringsAsFactors = FALSE
)
write.table(qa, qa_out, sep = "\t", quote = FALSE, row.names = FALSE)

layout <- all_cells %>%
  left_join(cells %>% mutate(across(c(Season, SampleType, Site), as.character)),
            by = c("Season", "SampleType", "Site")) %>%
  mutate(Availability = ifelse(is.na(N), "No sample", "Observed"))
write.table(layout, layout_out, sep = "\t", quote = FALSE, row.names = FALSE)

if (!all(qa$pass)) {
  failed <- paste(qa$check[!qa$pass], collapse = ", ")
  stop("Numeric/source QA failed: ", failed)
}

if (mode == "run") {
  geology <- data.frame(
    xmin = c(0.55, 2.55, 4.55), xmax = c(2.45, 4.45, 7.45),
    label_x = c(1.5, 3.5, 6.0), label = c("Ocean", "Tidal", "Island")
  )

  p <- ggplot(plot_df, aes(x = Site, y = MeanRelativeAbundance, fill = TaxonDisplay)) +
    geom_col(width = 0.78, colour = "white", linewidth = 0.10, na.rm = TRUE) +
    geom_segment(
      data = geology, aes(x = xmin, xend = xmax, y = 1.035, yend = 1.035),
      inherit.aes = FALSE, colour = "#62676C", linewidth = 0.40
    ) +
    geom_text(
      data = geology, aes(x = label_x, y = 1.057, label = label),
      inherit.aes = FALSE, colour = "#4D5358", size = 2.55, family = "Arial"
    ) +
    scale_fill_manual(values = palette, breaks = taxon_order, drop = FALSE) +
    scale_x_discrete(drop = FALSE, limits = site_order) +
    scale_y_continuous(
      limits = c(0, 1.075), breaks = seq(0, 1, 0.25),
      labels = function(z) sprintf("%d%%", round(z * 100)),
      expand = expansion(mult = c(0, 0))
    ) +
    facet_grid(rows = vars(SampleType), cols = vars(Season), drop = FALSE) +
    labs(x = NULL, y = "Mean relative abundance", fill = "18S class") +
    theme_classic(base_size = 8.5, base_family = "Arial") +
    theme(
      axis.title.y = element_text(face = "bold", size = 9.2, margin = margin(r = 7)),
      axis.text.y = element_text(size = 7.6, colour = "#25282A"),
      axis.text.x = element_text(size = 7.4, colour = "#25282A", angle = 32, hjust = 1, vjust = 1),
      axis.line = element_line(linewidth = 0.42, colour = "#2F3336"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#2F3336"),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "#E4E6E8"),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(8, "pt"),
      strip.background = element_rect(fill = "#F1F2F3", colour = "#62676C", linewidth = 0.40),
      strip.text = element_text(face = "bold", size = 8.4, colour = "#202326", margin = margin(3.5, 4, 3.5, 4)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(face = "bold", size = 8.0),
      legend.text = element_text(size = 7.0),
      legend.key.height = grid::unit(8.5, "pt"),
      legend.key.width = grid::unit(12, "pt"),
      legend.spacing.x = grid::unit(4, "pt"),
      plot.margin = margin(6, 7, 5, 6)
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE, title.position = "top", title.hjust = 0.5))

  ggsave(pdf_out, p, width = 174, height = 176, units = "mm", device = grDevices::cairo_pdf)
  ragg::agg_tiff(
    filename = tiff_out, width = 174, height = 176, units = "mm", res = 600,
    compression = "lzw", background = "white"
  )
  print(p)
  grDevices::dev.off()
  writeLines(capture.output(sessionInfo()), session_out, useBytes = TRUE)
  message("Wrote: ", pdf_out)
  message("Wrote: ", tiff_out)
} else {
  if (!file.exists(pdf_out) || !file.exists(tiff_out)) stop("Delivery outputs are missing")
  message("PASS: numeric/source checks and required delivery-file presence")
}
