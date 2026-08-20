#!/usr/bin/env Rscript

# Presentation-only update of Supplementary Figure S3.
# No ecological model is fitted. Panel c decomposes the assigned PR2 classes
# grouped as "Other" in Figure 1b using the accepted common-75 class table.

options(stringsAsFactors = FALSE, warn = 1)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
if (!length(file_arg)) stop("Run with Rscript", call. = FALSE)
SCRIPT <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(SCRIPT), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)

PKG <- file.path(ROOT, "revision", "re_analysis", "63_figure_presentation_refinement_R_v1", "generated", "s2")
SRC29 <- file.path(ROOT, "revision", "re_analysis", "29_dual_backend_figure_sources_v1", "source_data")
COMP <- file.path(ROOT, "revision", "re_analysis", "07_taxonomic_composition_v1", "results")
TAX18 <- file.path(ROOT, "revision", "re_analysis", "03_18s_rarefied_50000_v1", "results",
                   "18s_taxonomy_pr2_primary_rarefied_50000_analysis_v1.tsv")
DELIVERY <- file.path(PKG, "delivery")
STEM <- "Supplementary_Figure_S3_taxonomic_bubble_R_v2"

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ragg)
})

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
read_tsv <- function(path) {
  assert(file.exists(path), paste0("Missing source: ", path))
  read.delim(path, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
}

SITE_ORDER <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
SAMPLE_TYPE <- c(`Sea-sand` = "#D98C2B", Water = "#2C7FB8")

theme_paper <- function(base = 7.2) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = rel(1.25)),
      plot.margin = margin(3, 5, 3, 5, unit = "pt")
    )
}

bubble16_path <- file.path(SRC29, "Supplementary_Figure_S2_16S_sample_taxonomy_bubble_source_v1.tsv")
bubble18_path <- file.path(SRC29, "Supplementary_Figure_S2_18S_sample_taxonomy_bubble_source_v1.tsv")
class_long_path <- file.path(COMP, "18s_full_class_sample_relative_abundance_long_v1.tsv")
top_ledger_path <- file.path(COMP, "top_taxa_selection_ledger_v1.tsv")

bubble <- bind_rows(read_tsv(bubble16_path), read_tsv(bubble18_path))
class_long <- read_tsv(class_long_path)
top_ledger <- read_tsv(top_ledger_path)
tax18 <- read_tsv(TAX18)

needed_bubble <- c("Marker", "SampleID", "Site", "SampleType", "Season", "Taxon", "TaxonOrder", "RelativeAbundancePct")
assert(all(needed_bubble %in% names(bubble)), "Bubble source columns changed")
assert(all(c("Marker", "Dataset", "Rank", "SampleID", "Taxon", "RelativeAbundance") %in% names(class_long)),
       "18S class-long source columns changed")
assert(all(c("Marker", "Dataset", "Rank", "TaxonDisplay", "CategoryOrder", "CategoryRole", "GlobalMeanRelativeAbundance") %in% names(top_ledger)),
       "Top-taxa ledger columns changed")
assert(all(c("ASV_ID", "Class") %in% names(tax18)), "18S taxonomy source columns changed")

assert(length(unique(class_long$SampleID)) == 75L, "Expected 75 samples in the accepted 18S class table")

top12 <- top_ledger %>%
  filter(Marker == "18S", Dataset == "full_pr2_primary", Rank == "Class", CategoryRole == "top_assigned_taxon") %>%
  arrange(as.integer(CategoryOrder)) %>%
  pull(TaxonDisplay)
assert(length(top12) == 12L, "Expected exactly 12 individually displayed 18S classes in Figure 1b")

other_taxa <- class_long %>%
  filter(Marker == "18S", Dataset == "full_pr2_primary", Rank == "Class", Taxon != "Unassigned", !Taxon %in% top12) %>%
  group_by(Taxon) %>%
  summarise(MeanRelativeAbundance = mean(RelativeAbundance), .groups = "drop") %>%
  arrange(desc(MeanRelativeAbundance), Taxon)

assert(nrow(other_taxa) == 172L, "Figure 1b Other must contain 172 assigned PR2 classes")
other_total <- sum(other_taxa$MeanRelativeAbundance)
assert(abs(other_total - 0.4452552) < 1e-12, "Figure 1b Other mean relative abundance drifted")

leading <- other_taxa %>% slice_head(n = 20L)
residual <- other_taxa %>% slice_tail(n = nrow(other_taxa) - 20L)
asv_counts <- tax18 %>% count(Class, name = "ASVCount")
leading <- leading %>% left_join(asv_counts, by = c("Taxon" = "Class"))
residual_asv_n <- sum(asv_counts$ASVCount[asv_counts$Class %in% residual$Taxon])
panel_c <- bind_rows(
  leading %>% transmute(
    DisplayClass = Taxon,
    Category = "Leading constituent class",
    ConstituentClassCount = 1L,
    ASVCount = ASVCount,
    MeanRelativeAbundancePct = 100 * MeanRelativeAbundance
  ),
  data.frame(
    DisplayClass = "Remaining 152 assigned classes",
    Category = "Residual Other",
    ConstituentClassCount = nrow(residual),
    ASVCount = residual_asv_n,
    MeanRelativeAbundancePct = 100 * sum(residual$MeanRelativeAbundance),
    stringsAsFactors = FALSE
  )
) %>%
  mutate(
    ShareOfFigure1OtherPct = 100 * MeanRelativeAbundancePct / (100 * other_total),
    DisplayOrder = rank(-MeanRelativeAbundancePct, ties.method = "first")
  ) %>%
  arrange(DisplayOrder)

assert(sum(panel_c$ConstituentClassCount) == 172L, "Panel c class accounting failed")
assert(sum(panel_c$ASVCount) == 9302L, "Panel c ASV accounting failed")
assert(sum(leading$ASVCount) == 3943L && residual_asv_n == 5359L, "Panel c leading/residual ASV accounting failed")
assert(abs(sum(panel_c$MeanRelativeAbundancePct) - 44.52552) < 1e-10, "Panel c abundance accounting failed")
assert(abs(sum(leading$MeanRelativeAbundance) * 100 - 28.8548) < 1e-10, "Top-20 abundance drifted")
assert(abs(sum(residual$MeanRelativeAbundance) * 100 - 15.67072) < 1e-10, "Residual abundance drifted")

source_out <- file.path(PKG, "source_data", "Supplementary_Figure_S3c_18S_Other_decomposition_v1.tsv")
assert(!file.exists(source_out), paste0("Refusing overwrite: ", source_out))
write.table(panel_c[, c("DisplayOrder", "DisplayClass", "Category", "ConstituentClassCount", "ASVCount",
                         "MeanRelativeAbundancePct", "ShareOfFigure1OtherPct")],
            source_out, sep = "\t", quote = FALSE, row.names = FALSE)

plot_bubble_marker <- function(x, marker) {
  d <- x %>%
    filter(Marker == marker) %>%
    mutate(
      Site = factor(Site, SITE_ORDER),
      SampleID = factor(SampleID, unique(SampleID[order(Site, Season, SampleType, SampleID)])),
      TaxonKey = paste(sprintf("%02d", as.integer(TaxonOrder)), Taxon, sep = "__")
    )
  d$TaxonKey <- factor(d$TaxonKey, levels = rev(unique(d$TaxonKey[order(as.integer(d$TaxonOrder))])))
  rank_label <- if (marker == "16S") "16S genus" else "18S PR2 class"

  ggplot(d, aes(SampleID, TaxonKey)) +
    geom_point(aes(size = RelativeAbundancePct, fill = SampleType), shape = 21,
               color = "black", stroke = 0.12, alpha = 0.86) +
    facet_grid(Season ~ Site, scales = "free", space = "free_x") +
    scale_y_discrete(labels = function(z) sub("^[0-9]+__", "", z)) +
    scale_fill_manual(values = SAMPLE_TYPE) +
    scale_size_area(max_size = 5.1, limits = c(0, max(bubble$RelativeAbundancePct))) +
    labs(x = NULL, y = rank_label, size = "Relative abundance (%)", fill = "Sample type") +
    theme_paper(6.4) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4.5),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom",
      strip.text.x = element_text(size = 5.7),
      strip.text.y = element_text(size = 6.2)
    )
}

plot_other <- function(x) {
  d <- x %>%
    arrange(MeanRelativeAbundancePct) %>%
    mutate(DisplayClass = factor(DisplayClass, levels = DisplayClass))
  cols <- c(`Leading constituent class` = "#3C8DAD", `Residual Other` = "#BDBDBD")

  ggplot(d, aes(MeanRelativeAbundancePct, DisplayClass, fill = Category)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.18) +
    geom_text(aes(label = sprintf("%.2f", MeanRelativeAbundancePct)),
              hjust = -0.12, size = 2.1, color = "black") +
    scale_fill_manual(values = cols, breaks = names(cols), guide = "none") +
    scale_x_continuous(limits = c(0, 17.7), expand = expansion(mult = c(0, 0))) +
    labs(
      x = "Mean relative abundance across 75 samples (%)",
      y = NULL,
      subtitle = "Figure 1b Other = 44.53%: 20 leading classes = 28.85%; remaining 152 classes = 15.67%"
    ) +
    theme_paper(7.0) +
    theme(
      plot.subtitle = element_text(face = "bold", size = 7.0),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(3, 5, 3, 95, unit = "pt")
    )
}

figure <- (
  plot_bubble_marker(bubble, "16S") /
    plot_bubble_marker(bubble, "18S") /
    plot_other(panel_c)
) +
  plot_layout(heights = c(1, 1, 0.82), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

dir.create(file.path(PKG, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(DELIVERY, recursive = TRUE, showWarnings = FALSE)

save_outputs <- function(root) {
  pdf <- file.path(root, paste0(STEM, ".pdf"))
  tif <- file.path(root, paste0(STEM, ".tiff"))
  assert(!file.exists(pdf) && !file.exists(tif), paste0("Refusing overwrite in ", root))

  cairo_pdf(pdf, width = 210 / 25.4, height = 285 / 25.4, family = "Arial", onefile = TRUE)
  print(figure)
  dev.off()

  ragg::agg_tiff(tif, width = 210, height = 285, units = "mm", res = 600,
                 compression = "lzw", background = "white", bitsize = 8)
  print(figure)
  dev.off()
  assert(all(file.info(c(pdf, tif))$size > 0), "Empty figure output")
  c(pdf, tif)
}

pkg_outputs <- save_outputs(file.path(PKG, "figures"))
delivery_outputs <- file.path(DELIVERY, basename(pkg_outputs))
assert(!any(file.exists(delivery_outputs)), "Refusing overwrite in delivery directory")
copied <- mapply(file.copy, pkg_outputs, delivery_outputs,
                 MoreArgs = list(overwrite = FALSE, copy.date = TRUE))
assert(all(copied), "Delivery copy failed")
assert(all(file.info(pkg_outputs)$size == file.info(delivery_outputs)$size), "Delivery size mismatch")

qa <- data.frame(
  Check = c(
    "Accepted sample count", "Figure 1b Other classes", "Leading classes displayed",
    "Residual classes", "Figure 1b Other ASVs", "Leading-class ASVs", "Residual-class ASVs",
    "Figure 1b Other mean abundance (%)",
    "Leading 20 mean abundance (%)", "Residual mean abundance (%)",
    "PDF bytes", "TIFF bytes"
  ),
  Value = c(
    length(unique(class_long$SampleID)), nrow(other_taxa), nrow(leading), nrow(residual),
    sum(panel_c$ASVCount), sum(leading$ASVCount), residual_asv_n,
    sprintf("%.6f", 100 * other_total), sprintf("%.6f", 100 * sum(leading$MeanRelativeAbundance)),
    sprintf("%.6f", 100 * sum(residual$MeanRelativeAbundance)),
    file.info(pkg_outputs[1])$size, file.info(pkg_outputs[2])$size
  ),
  Status = "PASS",
  stringsAsFactors = FALSE
)
write.table(qa, file.path(PKG, "qc", "Supplementary_Figure_S3_other_decomposition_numeric_QA_v1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Created: ", paste(c(pkg_outputs, delivery_outputs, source_out), collapse = "\n"))
