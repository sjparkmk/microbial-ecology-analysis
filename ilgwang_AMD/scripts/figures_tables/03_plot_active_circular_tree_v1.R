suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(ggtreeExtra)
  library(ggnewscale)
  library(cowplot)
  library(ragg)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args) >= 1L) args[[1]] else "/input/youre/data/"
output_root <- if (length(args) >= 2L) args[[2]] else "/output/your/data"
tree_path <- file.path(
  source_root, "output_data", "asv_tree_all663.newick"
)
annotation_path <- file.path(
  source_root, "source_data", "Fig4_circular_ASV_tree_annotation_A4_v1.tsv"
)
figure_dir <- file.path(output_root, "build", "tree_stage", "figures")
log_dir <- file.path(output_root, "logs")
v2_log_dir <- file.path(output_root, "build", "tree_stage", "logs")
source_dir <- file.path(output_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(v2_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "ape", "dplyr", "ggplot2", "ggtree", "ggtreeExtra", "ggnewscale",
  "cowplot", "ragg", "svglite"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages (no automatic installation attempted): ", paste(missing_packages, collapse = ", "))
}
if (!file.exists(tree_path) || !file.exists(annotation_path)) {
  stop("Missing tree or annotation input")
}

tree <- ape::read.tree(tree_path)
annotation <- read.delim(
  annotation_path, header = TRUE, sep = "\t", check.names = FALSE,
  stringsAsFactors = FALSE, na.strings = c("NA", "")
)

if (length(tree$tip.label) != 663L) {
  stop("Expected 663 tips in the full tree; found ", length(tree$tip.label))
}
if (nrow(annotation) != 56L || anyDuplicated(annotation$OTUID)) {
  stop("Expected 56 unique annotation rows")
}
if (!all(annotation$OTUID %in% tree$tip.label)) {
  stop("At least one displayed ASV is absent from the full tree")
}
if (any(annotation$pct_pooled < 0.2)) {
  stop("At least one displayed ASV falls below the pooled 0.2% threshold")
}

pruned_tree <- ape::keep.tip(tree, annotation$OTUID)
if (length(pruned_tree$tip.label) != 56L) {
  stop("Tree pruning did not produce exactly 56 tips")
}
pruned_path <- file.path(source_dir, "Fig4_circular_ASV_tree_pruned56_A4_v2.newick")
ape::write.tree(pruned_tree, file = pruned_path)

phylum_colors <- c(
  "Acidobacteriota" = "#C0392B",
  "Actinomycetota" = "#1E8449",
  "Bacillota--other" = "#E6B35A",
  "Chloroflexota--other" = "#A98467",
  "Elusimicrobiota" = "#E8A0BF",
  "Nanobdellota" = "#6B4A3A",
  "Nitrospirota" = "#E07B1A",
  "Planctomycetota" = "#8FCE8F",
  "Pseudomonadota" = "#2166AC",
  "Thermoplasmatota" = "#A9C6E8",
  "Thermoproteota" = "#F08A7E",
  "Unclassified Bacteria" = "#8A4A1A",
  "Vulcanimicrobiota" = "#B9A0D4"
)
observed_phyla <- unique(annotation$Phylum_display)
unknown_phyla <- setdiff(observed_phyla, names(phylum_colors))
if (length(unknown_phyla) > 0L) {
  stop("No fixed color defined for: ", paste(unknown_phyla, collapse = ", "))
}

type_order <- c("P_biofilm", "W_biofilm", "water")
type_colors <- c(
  "P_biofilm" = "#1BAF7A",
  "W_biofilm" = "#EB6834",
  "water" = "#2A78D6"
)
ring_position <- c(
  "P_biofilm" = "inner",
  "W_biofilm" = "middle",
  "water" = "outer"
)

normalize_01 <- function(values) {
  value_range <- range(values, finite = TRUE)
  if (!all(is.finite(value_range))) {
    stop("Non-finite abundance in a Type ring")
  }
  if (value_range[[2]] == value_range[[1]]) {
    return(rep(0, length(values)))
  }
  (values - value_range[[1]]) / (value_range[[2]] - value_range[[1]])
}

tax_data <- annotation[, c("OTUID", "Phylum_display")]
tax_data$label <- tax_data$OTUID
tax_data$label_display <- ifelse(tax_data$OTUID == "ASV_1", "ASV_1 *", tax_data$OTUID)
tree_joined <- dplyr::full_join(pruned_tree, tax_data, by = "label")

tree_plot <- ggtree(
  tree_joined,
  layout = "circular",
  open.angle = 8,
  branch.length = "none",
  linewidth = 0.36,
  color = "black",
  show.legend = FALSE
) +
  geom_tiplab(
    aes(label = label_display, color = Phylum_display),
    size = 2.65,
    fontface = 2,
    offset = 0.46,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = phylum_colors,
    breaks = names(phylum_colors),
    na.value = "#555555",
    guide = "none"
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

offsets <- c(0.42, 0.11, 0.11)
for (index in seq_along(type_order)) {
  type_name <- type_order[[index]]
  ring_data <- annotation[, c("OTUID", type_name)]
  names(ring_data) <- c("OTUID", "raw_abundance")
  ring_data$scaled_abundance <- normalize_01(ring_data$raw_abundance)

  tree_plot <- tree_plot +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      data = ring_data,
      geom = geom_tile,
      mapping = aes(y = OTUID, fill = scaled_abundance),
      offset = offsets[[index]],
      pwidth = 0.50,
      width = 1.5,
      height = 1.0,
      color = "#383838",
      linewidth = 0.25
    ) +
    scale_fill_gradient(
      low = "#F7F6F2",
      high = unname(type_colors[[type_name]]),
      limits = c(0, 1),
      guide = "none"
    )
}


legend_phyla <- names(phylum_colors)[names(phylum_colors) %in% observed_phyla]
legend_index <- seq_along(legend_phyla) - 1L
legend_tax <- data.frame(
  Phylum = legend_phyla,
  column = legend_index %% 3L,
  row = 4L - (legend_index %/% 3L),
  stringsAsFactors = FALSE
)
tax_legend <- ggplot(legend_tax) +
  geom_point(
    aes(x = column * 1.15, y = row, fill = Phylum),
    shape = 22, size = 3.45, color = "#2D2D2D", stroke = 0.40
  ) +
  geom_text(
    aes(x = column * 1.15 + 0.10, y = row, label = Phylum),
    hjust = 0, size = 3.05, family = "Arial"
  ) +
  scale_fill_manual(values = phylum_colors, guide = "none") +
  annotate("text", x = 0, y = 5.0, label = "Phylum", hjust = 0, size = 3.45, fontface = "bold", family = "Arial") +
  coord_cartesian(xlim = c(-0.05, 3.45), ylim = c(-0.35, 5.25), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0, 0.5, 0, 0.5, unit = "mm"))


gradient_rows <- list()
for (row_index in seq_along(type_order)) {
  type_name <- type_order[[row_index]]
  palette <- grDevices::colorRampPalette(c("#F7F6F2", type_colors[[type_name]]))(101)
  gradient_rows[[row_index]] <- data.frame(
    x = seq(0, 1, length.out = 101),
    y = 4 - row_index,
    fill_value = palette,
    Type = type_name,
    stringsAsFactors = FALSE
  )
}
gradient_data <- do.call(rbind, gradient_rows)
ring_labels <- data.frame(
  Type = type_order,
  y = 3:1,
  text = paste0(type_order, " (", unname(ring_position[type_order]), ")"),
  stringsAsFactors = FALSE
)
ring_ticks <- do.call(rbind, lapply(seq_along(type_order), function(index) {
  data.frame(x = c(0, 0.5, 1), y = 4 - index - 0.35, label = c("0", "0.5", "1"))
}))
ring_legend <- ggplot() +
  geom_tile(
    data = gradient_data,
    aes(x = x, y = y, fill = fill_value),
    width = 0.011, height = 0.34
  ) +
  scale_fill_identity() +
  geom_text(
    data = ring_labels,
    aes(x = -0.06, y = y, label = text),
    hjust = 1, size = 3.05, family = "Arial", fontface = "bold"
  ) +
  geom_text(
    data = ring_ticks,
    aes(x = x, y = y, label = label),
    size = 2.75, family = "Arial"
  ) +
  annotate(
    "text", x = -0.55, y = 4.05,
    label = "Ring intensity (within-Type min-max scale)",
    hjust = 0, size = 3.45, fontface = "bold", family = "Arial"
  ) +
  coord_cartesian(xlim = c(-0.55, 1.04), ylim = c(0.25, 4.25), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0, 0.5, 0, 0.5, unit = "mm"))

legend_panel <- cowplot::plot_grid(
  tax_legend, ring_legend,
  ncol = 2,
  rel_widths = c(0.60, 0.40),
  align = "h"
)
final_plot <- cowplot::plot_grid(
  tree_plot,
  legend_panel,
  ncol = 1,
  rel_heights = c(0.70, 0.30),
  align = "v"
)

width_mm <- 210
height_mm <- 297
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
base_path <- file.path(figure_dir, "Fig4_circular_ASV_tree_A4_v2")

pdf_path <- paste0(base_path, ".pdf")
grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, bg = "white", family = "Arial")
print(final_plot)
grDevices::dev.off()

svg_path <- paste0(base_path, ".svg")
svglite::svglite(svg_path, width = width_in, height = height_in, bg = "white", system_fonts = list(sans = "Arial"))
print(final_plot)
grDevices::dev.off()

png_path <- paste0(base_path, ".png")
ragg::agg_png(
  png_path, width = width_mm, height = height_mm, units = "mm", res = 300,
  background = "white", scaling = 1
)
print(final_plot)
grDevices::dev.off()

tiff_path <- paste0(base_path, ".tiff")
ragg::agg_tiff(
  tiff_path, width = width_mm, height = height_mm, units = "mm", res = 400,
  background = "white", compression = "lzw", scaling = 1
)
print(final_plot)
grDevices::dev.off()

session_path <- file.path(v2_log_dir, "04_tree_R_sessionInfo_A4_v2.txt")
session_lines <- capture.output(sessionInfo())
writeLines(session_lines, session_path, useBytes = TRUE)

register_path <- file.path(v2_log_dir, "04_tree_outputs_A4_v2.tsv")
output_paths <- c(pdf_path, svg_path, png_path, tiff_path, pruned_path)
register <- data.frame(
  output = sub(paste0("^", gsub("([\\.\\+\\*\\?\\[\\^\\]\\$\\(\\)\\{\\}=!<>|:\\-])", "\\\\\\1", output_root), "/?"), "", output_paths),
  bytes = file.info(output_paths)$size,
  stringsAsFactors = FALSE
)
write.table(register, register_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("FULL TREE TIPS:", length(tree$tip.label), "\n")
cat("DISPLAYED TIPS:", length(pruned_tree$tip.label), "\n")
cat("POOLED ABUNDANCE THRESHOLD: >=0.2%\n")
cat("TREE DISPLAY: circular topology only; branch lengths suppressed; inherited ASV_1 root not interpreted\n")
cat("RINGS:", paste(paste0(type_order, "=", unname(ring_position[type_order])), collapse = "; "), "\n")
cat("OUTPUTS:\n", paste(output_paths, collapse = "\n"), "\n")
cat("R CIRCULAR TREE GENERATION: PASS\n")
