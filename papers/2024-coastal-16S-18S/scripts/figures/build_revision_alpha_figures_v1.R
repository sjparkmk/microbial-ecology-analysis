suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

project_root <- "/user/link/path" # insert your path
alpha_root <- file.path(project_root, "revision", "re_analysis", "05_alpha_diversity_v3", "results")
package_root <- file.path(project_root, "revision", "re_analysis", "15_revision_alpha_results_v1")
source_dir <- file.path(package_root, "source_data")
internal_figure_dir <- file.path(package_root, "figures")
qa_dir <- file.path(package_root, "qc")
delivery_dir <- file.path(project_root, "revision", "rev_ms", "figures")

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(internal_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(delivery_dir, recursive = TRUE, showWarnings = FALSE)

sample_path <- file.path(alpha_root, "alpha_sample_metrics_v3.tsv")
summary_path <- file.path(alpha_root, "alpha_cell_summaries_v3.tsv")
test_path <- file.path(alpha_root, "alpha_exact_signflip_results_v3.tsv")

marker_cols <- cols(Marker = col_character(), .default = col_guess())

samples <- read_tsv(sample_path, col_types = marker_cols, show_col_types = FALSE, progress = FALSE) %>%
  filter(CountSurface == "rarefied_primary")
summaries <- read_tsv(summary_path, col_types = marker_cols, show_col_types = FALSE, progress = FALSE)
tests <- read_tsv(test_path, col_types = marker_cols, show_col_types = FALSE, progress = FALSE)

stopifnot(nrow(samples) == 900L)
stopifnot(nrow(summaries) == 300L)
stopifnot(nrow(tests) == 60L)
stopifnot(all(samples %>% count(Marker, Metric) %>% pull(n) == 75L))
stopifnot(all(samples %>% distinct(Marker, SampleID) %>% count(Marker) %>% pull(n) == 75L))
stopifnot(all(summaries$N == 3L))

site_levels <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
metric_levels <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
season_levels <- c("Spring", "Autumn")
sample_type_levels <- c("Sea-sand", "Water")

panel_keys <- unlist(lapply(season_levels, function(season) paste(season, metric_levels, sep = "|")))
panel_metric_labels <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Inverse Simpson")
panel_labels <- character(length(panel_keys))
letters12 <- letters[seq_along(panel_keys)]
for (i in seq_along(panel_keys)) {
  pieces <- strsplit(panel_keys[i], "\\|", fixed = FALSE)[[1]]
  metric_index <- match(pieces[2], metric_levels)
  panel_labels[i] <- paste0(letters12[i], "  ", pieces[1], ": ", panel_metric_labels[metric_index])
}
names(panel_labels) <- panel_keys

prepare_factors <- function(data) {
  data %>%
    mutate(
      Site = factor(Site, levels = site_levels),
      Metric = factor(Metric, levels = metric_levels),
      Season = factor(Season, levels = season_levels),
      SampleType = factor(SampleType, levels = sample_type_levels),
      GeologicalHabitat = factor(GeologicalHabitat, levels = c("Ocean", "Tidal", "Island")),
      PanelKey = paste(as.character(Season), as.character(Metric), sep = "|"),
      Panel = factor(PanelKey, levels = panel_keys, labels = unname(panel_labels))
    )
}

samples <- prepare_factors(samples)
summaries <- prepare_factors(summaries)

write_tsv(samples, file.path(source_dir, "alpha_figure_sample_metrics_v1.tsv"), na = "")
write_tsv(summaries, file.path(source_dir, "alpha_figure_cell_mean_sd_v1.tsv"), na = "")
write_tsv(tests, file.path(source_dir, "alpha_figure_exact_tests_v1.tsv"), na = "")

habitat_palette <- c(Ocean = "#0072B2", Tidal = "#E69F00", Island = "#009E73")
sample_fill_palette <- c(`Sea-sand` = "#D8C3A5", Water = "#C8DCEB")
sample_shapes <- c(`Sea-sand` = 21, Water = 24)

theme_revision_alpha <- function() {
  theme_classic(base_size = 7.2, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#2B2B2B"),
      axis.ticks = element_line(linewidth = 0.3, colour = "#2B2B2B"),
      axis.text.x = element_text(size = 5.7, angle = 48, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 6.0, colour = "#2B2B2B"),
      axis.title.y = element_text(size = 7.2, margin = margin(r = 4)),
      strip.background = element_rect(fill = "#F4F4F4", colour = "#4D4D4D", linewidth = 0.35),
      strip.text = element_text(size = 6.4, face = "bold", colour = "#222222", margin = margin(2.1, 1.5, 2.1, 1.5)),
      panel.grid.major.y = element_line(linewidth = 0.22, colour = "#E1E1E1"),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(2.6, "mm"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 6.5, face = "bold"),
      legend.text = element_text(size = 6.2),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      plot.margin = margin(3.5, 3.5, 2.5, 3.5, unit = "mm")
    )
}

make_alpha_plot <- function(marker) {
  d <- samples %>% filter(Marker == marker)
  s <- summaries %>% filter(Marker == marker)
  stopifnot(nrow(d) == 450L)
  stopifnot(nrow(s) == 150L)

  ggplot(d, aes(x = Site, y = Value)) +
    geom_errorbar(
      data = s,
      aes(
        x = Site, ymin = pmax(Mean - SD, 0), ymax = Mean + SD,
        group = SampleType, colour = GeologicalHabitat
      ),
      position = position_dodge(width = 0.62),
      width = 0.16, linewidth = 0.38, alpha = 0.95,
      inherit.aes = FALSE
    ) +
    geom_point(
      aes(
        colour = GeologicalHabitat, fill = SampleType,
        shape = SampleType, group = SampleType
      ),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.62),
      size = 1.0, stroke = 0.28, alpha = 0.52
    ) +
    geom_point(
      data = s,
      aes(
        x = Site, y = Mean, colour = GeologicalHabitat, fill = SampleType,
        shape = SampleType, group = SampleType
      ),
      position = position_dodge(width = 0.62),
      size = 1.85, stroke = 0.62, alpha = 1,
      inherit.aes = FALSE
    ) +
    facet_wrap(vars(Panel), scales = "free_y", ncol = 3, drop = FALSE) +
    scale_colour_manual(name = "Geological context", values = habitat_palette, drop = FALSE) +
    scale_fill_manual(name = "Sample type", values = sample_fill_palette, drop = FALSE) +
    scale_shape_manual(name = "Sample type", values = sample_shapes, drop = FALSE) +
    labs(x = NULL, y = "Alpha-diversity value") +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(shape = 16, fill = NA, size = 2.0, alpha = 1, linewidth = 0)
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(colour = "#333333", size = 2.0, alpha = 1)
      ),
      fill = "none"
    ) +
    theme_revision_alpha()
}

save_marker <- function(marker, figure_number) {
  plot <- make_alpha_plot(marker)
  stem <- paste0("Figure_", figure_number, "_", marker, "_alpha_diversity_v1")
  pdf_path <- file.path(delivery_dir, paste0(stem, ".pdf"))
  tiff_path <- file.path(delivery_dir, paste0(stem, ".tiff"))
  svg_path <- file.path(internal_figure_dir, paste0(stem, ".svg"))
  png_path <- file.path(qa_dir, paste0(stem, "_preview.png"))

  ggsave(
    filename = pdf_path, plot = plot,
    device = grDevices::cairo_pdf,
    width = 174, height = 220, units = "mm", bg = "white"
  )

  ragg::agg_tiff(
    filename = tiff_path,
    width = 174, height = 220, units = "mm", res = 600,
    background = "white", compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    file = svg_path,
    width = 174 / 25.4, height = 220 / 25.4,
    bg = "white", system_fonts = list(sans = "Arial")
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    filename = png_path,
    width = 174, height = 220, units = "mm", res = 300,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()

  invisible(c(pdf = pdf_path, tiff = tiff_path, svg = svg_path, preview = png_path))
}

outputs_16s <- save_marker("16S", 2)
outputs_18s <- save_marker("18S", 3)

qa_summary <- tibble(
  Marker = c("16S", "18S"),
  Figure = c("Figure 2", "Figure 3"),
  SampleRows = c(sum(samples$Marker == "16S"), sum(samples$Marker == "18S")),
  UniqueSamples = c(n_distinct(samples$SampleID[samples$Marker == "16S"]), n_distinct(samples$SampleID[samples$Marker == "18S"])),
  SummaryRows = c(sum(summaries$Marker == "16S"), sum(summaries$Marker == "18S")),
  Metrics = c(n_distinct(samples$Metric[samples$Marker == "16S"]), n_distinct(samples$Metric[samples$Marker == "18S"])),
  RarefactionDepth = c(unique(samples$RarefactionDepth[samples$Marker == "16S"]), unique(samples$RarefactionDepth[samples$Marker == "18S"])),
  Panels = c(12L, 12L),
  WidthMM = c(174L, 174L),
  HeightMM = c(220L, 220L),
  TiffDPI = c(600L, 600L)
)
write_tsv(qa_summary, file.path(qa_dir, "alpha_figure_generation_summary_v1.tsv"))

message("Generated: ", paste(c(outputs_16s, outputs_18s), collapse = "\n"))
