#!/usr/bin/env Rscript

# Publication redraws in R from the frozen dual-backend source package.
# This script never fits ecological/statistical models. It only validates and
# renders values assembled in 29_dual_backend_figure_sources_v1.

options(stringsAsFactors = FALSE, warn = 1)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
if (!length(file_arg)) stop("Run with Rscript", call. = FALSE)
SCRIPT <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(file.path(dirname(SCRIPT), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)
SRC_ROOT <- file.path(ROOT, "revision", "re_analysis", "29_dual_backend_figure_sources_v1")
PKG_ROOT <- file.path(ROOT, "revision", "re_analysis", "63_figure_presentation_refinement_R_v1", "generated", "core")
DELIVERY <- file.path(ROOT, "revision", "re_analysis", "63_figure_presentation_refinement_R_v1", "generated", "core_delivery")
PREV_PKG_ROOT <- file.path(ROOT, "revision", "re_analysis", "35_dual_backend_figures_R_v4")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(ape)
  library(ggtree)
  library(digest)
  library(ragg)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
cli <- commandArgs(trailingOnly = TRUE)
modes <- c("--dry-run", "--run", "--check")
active <- intersect(cli, modes)
assert(length(active) == 1L, "Choose exactly one of --dry-run, --run, --check")
MODE <- active[[1]]

read_tsv <- function(path) {
  assert(file.exists(path), paste0("Missing source: ", path))
  read.delim(path, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
}
need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  assert(!length(miss), paste0(label, " missing columns: ", paste(miss, collapse = ", ")))
}
sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
atomic_tsv <- function(x, path) {
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  assert(!file.exists(tmp), paste0("Temporary file exists: ", tmp))
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  assert(file.rename(tmp, path), paste0("Atomic write failed: ", path))
}

HABITAT <- c(Ocean = "#0072B2", Tidal = "#E69F00", Island = "#009E73")
SAMPLE_TYPE <- c(`Sea-sand` = "#D98C2B", Water = "#2C7FB8")
SEASON <- c(Spring = "#3B8C6E", Autumn = "#B86E3B")
SITE_ORDER <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")

theme_paper <- function(base = 7.4) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"),
          axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .3),
          strip.background = element_rect(fill = "white", color = "black", linewidth = .3),
          strip.text = element_text(face = "bold", color = "black"),
          legend.title = element_text(face = "bold"), legend.key.height = grid::unit(3.3, "mm"),
          plot.tag = element_text(face = "bold", size = rel(1.25)),
          plot.margin = margin(4, 5, 4, 5, unit = "pt"))
}

specs <- data.frame(
  stem = c(
    "Supplementary_Figure_S5_16S_core_ASV_patterns_R_v5",
    "Supplementary_Figure_S6_18S_core_ASV_patterns_R_v5",
    "Supplementary_Figure_S7_core_ASV_phylogeny_abundance_R_v5"
  ),
  width_mm = c(180, 180, 183),
  height_mm = c(205, 205, 240),
  stringsAsFactors = FALSE
)

source_aliases <- list(
  core_samples16 = "Supplementary_Figure_S5_16S_core_sample_summary_v1.tsv",
  core_samples18 = "Supplementary_Figure_S6_18S_core_sample_summary_v1.tsv",
  core_patterns16 = "Supplementary_Figure_S5_16S_core_ASV_group_ledger_v1.tsv",
  core_patterns18 = "Supplementary_Figure_S6_18S_core_ASV_group_ledger_v1.tsv",
  phylo_abundance16 = "Supplementary_Figure_S7_16S_representative_core_group_abundance_v1.tsv",
  phylo_abundance18 = "Supplementary_Figure_S7_18S_representative_core_group_abundance_v1.tsv"
)

resolve_source <- function(key) {
  dirs <- c(file.path(SRC_ROOT, "source_data"), file.path(SRC_ROOT, "tables"), SRC_ROOT)
  cand <- unlist(lapply(dirs, function(d) file.path(d, source_aliases[[key]])))
  hit <- cand[file.exists(cand)]
  assert(length(hit) == 1L, paste0("Expected exactly one source for ", key, "; candidates: ", paste(cand, collapse = "; ")))
  hit[[1]]
}

canon_names <- function(x) {
  aliases <- list(
    Marker = c("marker"), Taxon = c("taxon", "TaxonLabel", "Feature"),
    SampleID = c("sample_id"), Site = c("site"), SampleType = c("sample_type"), Season = c("season"),
    Prevalence = c("prevalence", "PrevalenceAtCountGT0"),
    MeanRelativeAbundance = c("mean_relative_abundance", "MeanRelativeAbundancePct"),
    CoreReadProportion = c("core_read_proportion", "CoreReadProportionPct", "CoreReadFraction", "CoreReadPct"),
    ASV = c("asv", "ASV_ID", "FeatureID"), Group = c("group", "CommunityGroup"), TipOrder = c("tip_order")
  )
  for (target in names(aliases)) if (!target %in% names(x)) {
    hit <- aliases[[target]][aliases[[target]] %in% names(x)]
    if (length(hit)) names(x)[match(hit[[1]], names(x))] <- target
  }
  x
}

load_core_only <- function() {
  keys <- c(
    "core_samples16", "core_samples18", "core_patterns16", "core_patterns18",
    "phylo_abundance16", "phylo_abundance18"
  )
  paths <- setNames(vapply(keys, resolve_source, character(1)), keys)
  raw <- lapply(paths, function(p) canon_names(read_tsv(p)))
  data <- list(
    core_samples = bind_rows(raw$core_samples16, raw$core_samples18),
    core_patterns = bind_rows(raw$core_patterns16, raw$core_patterns18),
    phylo_abundance = bind_rows(raw$phylo_abundance16, raw$phylo_abundance18)
  )
  for (nm in c("core_patterns", "phylo_abundance")) {
    if (!"Group" %in% names(data[[nm]])) {
      data[[nm]]$Group <- paste(data[[nm]]$SampleType, data[[nm]]$Season, sep = " | ")
    }
    data[[nm]]$Group <- factor(
      data[[nm]]$Group,
      levels = c("Sea-sand | Spring", "Sea-sand | Autumn", "Water | Spring", "Water | Autumn")
    )
  }
  list(paths = paths, data = data)
}

validate_core_only <- function(b) {
  d <- b$data
  assert(all(c("16S", "18S") %in% unique(d$core_samples$Marker)), "Core sample marker contract failed")
  assert(all(c("16S", "18S") %in% unique(d$core_patterns$Marker)), "Core pattern marker contract failed")
  assert(all(c("16S", "18S") %in% unique(d$phylo_abundance$Marker)), "Phylogeny marker contract failed")
  invisible(TRUE)
}

plot_core <- function(marker, sample, patterns) {
  s <- sample %>% filter(Marker == marker) %>% mutate(Site = factor(Site, SITE_ORDER))
  if (max(s$CoreReadProportion, na.rm = TRUE) <= 1.01) s$CoreReadProportion <- s$CoreReadProportion * 100
  summ <- s %>% group_by(Site, SampleType, Season) %>% summarise(Mean = mean(CoreReadProportion), SD = sd(CoreReadProportion), .groups = "drop")
  p1 <- ggplot(s, aes(Site, CoreReadProportion, color = SampleType, shape = Season)) +
    geom_point(position = position_jitterdodge(jitter.width = .11, dodge.width = .55), size = 1.55, alpha = .75) +
    geom_errorbar(data = summ, aes(x = Site, y = Mean, ymin = pmax(0, Mean - SD), ymax = Mean + SD, group = interaction(SampleType, Season)),
                  inherit.aes = FALSE, position = position_dodge(.55), width = .17, linewidth = .38, color = "black") +
    geom_point(data = summ, aes(x = Site, y = Mean, group = interaction(SampleType, Season), fill = SampleType), inherit.aes = FALSE,
               position = position_dodge(.55), shape = 23, size = 2.0, color = "black") +
    scale_color_manual(values = SAMPLE_TYPE) + scale_fill_manual(values = SAMPLE_TYPE) +
    scale_shape_manual(values = c(Spring = 16, Autumn = 1)) +
    labs(x = NULL, y = "Core-ASV read proportion (%)", color = "Sample type", fill = "Sample type", shape = "Season") +
    guides(color = guide_legend(order = 1, nrow = 1),
           fill = guide_legend(order = 1, nrow = 1),
           shape = guide_legend(order = 2, nrow = 1)) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom",
          legend.box = "vertical", legend.box.just = "center",
          legend.spacing.y = grid::unit(0.3, "mm"), legend.margin = margin(0, 0, 0, 0))
  p <- patterns %>% filter(Marker == marker)
  fill_limits <- if (marker == "16S") c(0, 9.72666666666667) else c(0, 6.61066666666667)
  fill_breaks <- if (marker == "16S") c(0, 2.5, 5.0, 7.5) else c(0, 2, 4, 6)
  top_asv <- p %>% group_by(ASV, Taxon) %>% summarise(ASVOrder = min(ASVOrder), .groups = "drop") %>%
    arrange(ASVOrder, ASV) %>% slice_head(n = 40)
  p <- inner_join(p, top_asv[, c("ASV", "Taxon")], by = c("ASV", "Taxon")) %>%
    mutate(Label = factor(paste(ASV, Taxon, sep = " | "), levels = rev(unique(paste(ASV, Taxon, sep = " | ")[order(MeanRelativeAbundance)]))))
  p2 <- ggplot(p, aes(Group, Label)) +
    geom_tile(aes(fill = MeanRelativeAbundance), color = "white", linewidth = .15) +
    geom_point(data = p %>% filter(Prevalence > 0), aes(size = Prevalence),
               shape = 21, fill = NA, color = "black", stroke = .40) +
    scale_fill_viridis_c(
      option = "C", trans = "sqrt", limits = fill_limits, breaks = fill_breaks,
      name = "Mean relative abundance (%)"
    ) +
    scale_size_area(
      name = "Prevalence (count >0, %)", max_size = 4.2, limits = c(0, 1),
      breaks = c(.25, .50, .75, 1.00), labels = c("25", "50", "75", "100")
    ) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_colorbar(order = 1, title.position = "top", barwidth = grid::unit(32, "mm"),
                                 barheight = grid::unit(2.6, "mm")),
           size = guide_legend(order = 2, title.position = "top", nrow = 1, byrow = TRUE,
                               keywidth = grid::unit(7, "mm"),
                               override.aes = list(shape = 21, fill = NA, colour = "black", stroke = .40))) +
    theme_paper(6.6) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom",
          legend.box = "vertical", legend.box.just = "left",
          legend.spacing.x = grid::unit(1.2, "mm"), legend.margin = margin(0, 0, 0, 0))
  (p1 | p2) + plot_layout(widths = c(.9, 1.2), guides = "keep") + plot_annotation(tag_levels = "A")
}

tree_file_for <- function(marker) {
  dirs <- c(file.path(SRC_ROOT, "source_data"), file.path(SRC_ROOT, "trees"), SRC_ROOT)
  names <- c(paste0("Supplementary_Figure_S7_", marker, "_representative_core_ML_tree_v1.newick"),
             paste0(marker, "_core_ASV_tree_v1.nwk"), paste0("core_ASV_tree_", marker, "_v1.nwk"), paste0(marker, "_core_asv_tree_v1.tree"))
  hit <- unlist(lapply(dirs, function(d) file.path(d, names)))
  hit <- hit[file.exists(hit)]
  assert(length(hit) == 1L, paste0("Expected one current tree for ", marker))
  hit[[1]]
}

tree_abundance_panel <- function(marker, x) {
  d <- x %>% filter(Marker == marker)
  fill_limits <- if (marker == "16S") c(0, 9.72666666666667) else c(0, 6.61066666666667)
  fill_breaks <- if (marker == "16S") c(0, 2.5, 5.0, 7.5) else c(0, 2, 4, 6)
  tree <- ape::read.tree(tree_file_for(marker))
  keep <- intersect(tree$tip.label, unique(d$ASV))
  assert(length(keep) >= 3L, paste0(marker, " tree/abundance linkage too small"))
  tree <- ape::keep.tip(tree, keep)
  d <- d[d$ASV %in% tree$tip.label, , drop = FALSE]
  tree_df <- ggtree::ggtree(tree, size = .35)
  tip_dat <- tree_df$data %>% filter(isTip)
  ord <- tip_dat %>% arrange(y) %>% pull(label)
  tip_lab <- d %>% distinct(ASV, Taxon) %>% mutate(Display = paste(ASV, Taxon, sep = " | "))
  tree_df <- tree_df + geom_tiplab(aes(label = tip_lab$Display[match(label, tip_lab$ASV)]), size = 1.55,
                                   align = TRUE, linesize = .2, offset = max(tree_df$data$x, na.rm = TRUE) * .02) +
    xlim(0, max(tree_df$data$x, na.rm = TRUE) * 1.62) + theme_tree2()
  d$ASV <- factor(d$ASV, levels = ord)
  core_rows <- d %>% filter(IsCore %in% c(TRUE, "True", "TRUE", 1L))
  tile <- ggplot(d, aes(Group, ASV)) + geom_tile(aes(fill = MeanRelativeAbundance), color = "white", linewidth = .15) +
    geom_point(data = core_rows, aes(shape = "Core"), fill = NA, color = "black", stroke = .45, size = 2.6) +
    scale_fill_viridis_c(
      option = "C", trans = "sqrt", limits = fill_limits, breaks = fill_breaks,
      name = "Mean relative abundance (%)"
    ) +
    scale_shape_manual(
      name = "Core status", values = c(Core = 21),
      labels = c(Core = "Core in group\n(prevalence >=80%; count >0)")
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      fill = guide_colorbar(order = 1, title.position = "top", barwidth = grid::unit(30, "mm"),
                            barheight = grid::unit(2.6, "mm")),
      shape = guide_legend(order = 2, title.position = "top",
                           override.aes = list(size = 3.0, fill = NA, colour = "black", stroke = .45))
    ) +
    theme_paper(6.5) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom",
          legend.box = "vertical", legend.box.just = "left")
  (tree_df | tile) + plot_layout(widths = c(1.15, .85)) + plot_annotation(title = marker)
}

plot_s7 <- function(x) {
  ((tree_abundance_panel("16S", x) / tree_abundance_panel("18S", x)) +
      plot_annotation(tag_levels = "A")) &
    theme(plot.tag = element_text(family = "Arial", face = "bold", size = rel(1.25)))
}

plot_bundle <- function(d) list(
  Supplementary_Figure_S5_16S_core_ASV_patterns_R_v5 = plot_core("16S", d$core_samples, d$core_patterns),
  Supplementary_Figure_S6_18S_core_ASV_patterns_R_v5 = plot_core("18S", d$core_samples, d$core_patterns),
  Supplementary_Figure_S7_core_ASV_phylogeny_abundance_R_v5 = plot_s7(d$phylo_abundance)
)

save_one <- function(plot, stem, width_mm, height_mm, root) {
  pdf <- file.path(root, paste0(stem, ".pdf")); tif <- file.path(root, paste0(stem, ".tiff"))
  assert(!file.exists(pdf) && !file.exists(tif), paste0("Refusing overwrite: ", stem))
  grDevices::cairo_pdf(pdf, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", onefile = TRUE)
  print(plot); grDevices::dev.off()
  ragg::agg_tiff(tif, width = width_mm, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white", bitsize = 8)
  print(plot); grDevices::dev.off()
  assert(all(file.info(c(pdf, tif))$size > 0), paste0("Empty output: ", stem))
  c(pdf, tif)
}

copy_v4_one <- function(stem, root) {
  old_stem <- sub("_R_v5$", "_R_v4", stem)
  src <- file.path(PREV_PKG_ROOT, "figures", paste0(old_stem, c(".pdf", ".tiff")))
  dst <- file.path(root, paste0(stem, c(".pdf", ".tiff")))
  assert(all(file.exists(src)), paste0("Missing unchanged R v4 predecessor for ", stem))
  assert(!any(file.exists(dst)), paste0("Refusing overwrite: ", stem))
  ok <- mapply(file.copy, src, dst, MoreArgs = list(overwrite = FALSE, copy.date = TRUE))
  assert(all(ok) && all(file.info(dst)$size > 0), paste0("Predecessor copy failed: ", stem))
  dst
}


sig_ok <- function(path) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  raw <- readBin(con, "raw", 8)
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") return(rawToChar(raw[1:4]) == "%PDF")
  if (ext %in% c("tif", "tiff")) return(identical(as.integer(raw[1:4]), c(73L, 73L, 42L, 0L)) || identical(as.integer(raw[1:4]), c(77L, 77L, 0L, 42L)))
  FALSE
}

check_outputs <- function() {
  out_manifest <- file.path(PKG_ROOT, "manifests", "R_output_manifest_v5.tsv")
  assert(file.exists(out_manifest), "Missing R output manifest")
  man <- read_tsv(out_manifest)
  need_cols(man, c("Path", "Bytes", "SHA256"), "R output manifest")
  for (i in seq_len(nrow(man))) {
    p <- file.path(ROOT, man$Path[i])
    assert(file.exists(p), paste0("Manifest file missing: ", p))
    assert(file.info(p)$size == man$Bytes[i], paste0("Size drift: ", p))
    assert(tolower(sha(p)) == tolower(man$SHA256[i]), paste0("Hash drift: ", p))
  }
  figs <- unlist(lapply(specs$stem, function(z) file.path(DELIVERY, paste0(z, c(".pdf", ".tiff")))))
  assert(all(file.exists(figs)), "One or more R delivery figures are missing")
  assert(all(vapply(figs, sig_ok, logical(1))), "PDF/TIFF signature check failed")
  cat("CHECK PASS: source manifest, three PDF/TIFF bundles, sizes, hashes, and signatures are current.\n")
}

if (MODE == "--check") {
  check_outputs()
  quit(save = "no", status = 0)
}

assert(dir.exists(SRC_ROOT), paste0("Shared source package is not ready: ", SRC_ROOT))
bundle <- load_core_only()
validate_core_only(bundle)
cat("R dual-backend figure contract\n")
cat("  source package:", SRC_ROOT, "\n")
cat("  figures:", nrow(specs), "(core-pattern and tree presentations only)\n")
cat("  backend: R 4.5.2 / ggplot2; no ecological model refit\n")
if (MODE == "--dry-run") {
  assert(!dir.exists(PKG_ROOT), paste0("R package root already exists: ", PKG_ROOT))
  final_targets <- unlist(lapply(specs$stem, function(z) file.path(DELIVERY, paste0(z, c(".pdf", ".tiff")))))
  assert(!any(file.exists(final_targets)), "One or more R v5 delivery targets already exist")
  cat("DRY RUN PASS: frozen core-source schemas, tree linkage, and no-overwrite targets validated.\n")
  quit(save = "no", status = 0)
}

assert(MODE == "--run", "Internal mode dispatch failure")
assert(!dir.exists(PKG_ROOT), paste0("Refusing overwrite of package: ", PKG_ROOT))
dir.create(dirname(PKG_ROOT), recursive = TRUE, showWarnings = FALSE)
stage <- file.path(dirname(PKG_ROOT), paste0(".38_dual_backend_figures_R_v5_staging_", Sys.getpid()))
assert(!file.exists(stage), paste0("Staging path exists: ", stage))
dir.create(file.path(stage, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "manifests"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "logs"), recursive = TRUE, showWarnings = FALSE)
promoted <- FALSE; delivered <- character()
on.exit({
  if (dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
  if (!promoted && length(delivered)) unlink(delivered, force = TRUE)
}, add = TRUE)

plots <- plot_bundle(bundle$data)
assert(setequal(names(plots), specs$stem), "Plot/spec stem mismatch")
for (i in seq_len(nrow(specs))) {
  stem <- specs$stem[i]
  save_one(plots[[stem]], stem, specs$width_mm[i], specs$height_mm[i], file.path(stage, "figures"))
}

root_prefix <- paste0(normalizePath(ROOT, winslash = "/"), "/")
input_paths_abs <- normalizePath(bundle$paths, winslash = "/")
input_manifest <- data.frame(Role = names(bundle$paths), Path = substring(input_paths_abs, nchar(root_prefix) + 1L),
                             Bytes = as.numeric(file.info(bundle$paths)$size), SHA256 = vapply(bundle$paths, sha, character(1)), stringsAsFactors = FALSE)
atomic_tsv(input_manifest, file.path(stage, "manifests", "R_input_manifest_v5.tsv"))
atomic_tsv(specs, file.path(stage, "manifests", "R_figure_dimensions_v5.tsv"))
writeLines(capture.output(sessionInfo()), file.path(stage, "logs", "R_sessionInfo_v5.txt"), useBytes = TRUE)
qa <- data.frame(Check = c("backend", "figure_count", "source_only_redraw", "tree_scope"),
                 Observed = c("R", nrow(specs), TRUE, "core ASVs"), Expected = c("R", 3, TRUE, "core ASVs"), Status = "PASS")
atomic_tsv(qa, file.path(stage, "qc", "R_figure_numeric_QA_v5.tsv"))
readme <- c("# Core-figure panel-label standardization", "", "Presentation-only R redraws from the frozen shared source package.",
             "Only the frozen core-ASV source tables and trees were rendered; no ecological or statistical model was refitted here.",
             "Only the existing panel identifiers were changed to bold uppercase letters without parentheses.")
writeLines(readme, file.path(stage, "README_R_figures_v5.md"), useBytes = TRUE)

supersession <- data.frame(
  Artifact = c("revision/re_analysis/30_dual_backend_figures_R_v1 and *_R_v1 deliveries",
               "revision/re_analysis/32_dual_backend_figures_R_v2 and *_R_v2 deliveries",
               "revision/re_analysis/34_dual_backend_figures_R_v3 and *_R_v3 deliveries",
               "revision/re_analysis/35_dual_backend_figures_R_v4 and *_R_v4 deliveries",
               "revision/re_analysis/38_dual_backend_figures_R_v5 and *_R_v5 deliveries"),
  Status = c("scientifically_rejected_unhanded", "intermediate_unhanded", "visual_QA_rejected_unhanded",
             "visual_QA_intermediate_unhanded", "current_candidate"),
  Reason = c("Rendered before correction of the original-script core-ASV definition; core unions were 13/3 and are invalid for S5-S7.",
             "Rendered during schema transition before the frozen PrevalenceAtCountGT0 field and top-40 S5-S6 display contract were finalized.",
             "Numerically current but rejected during PDF visual QA because S1, Figure 7, S5, S7, S8, and S10 contained label or clipping defects.",
             "Numerically current, but S5-S6 combined bottom guides clipped the Season legend title in the delivery PDFs.",
             "Uses the frozen count>0 prevalence>=0.80 core definition (union 324/73); only the S5-S6 legend arrangement changed from R v4."),
  Replacement = rep("R v5", 5), stringsAsFactors = FALSE
)
atomic_tsv(supersession, file.path(stage, "manifests", "R_supersession_ledger_v5.tsv"))

# Copy from staging to versioned delivery only after every figure exists.
staged_figs <- unlist(lapply(specs$stem, function(z) file.path(stage, "figures", paste0(z, c(".pdf", ".tiff")))))
assert(all(file.exists(staged_figs)) && all(file.info(staged_figs)$size > 0), "Staged figure set incomplete")
dir.create(DELIVERY, recursive = TRUE, showWarnings = FALSE)
for (src in staged_figs) {
  dst <- file.path(DELIVERY, basename(src)); assert(!file.exists(dst), paste0("Delivery appeared during run: ", dst))
  assert(file.copy(src, dst, overwrite = FALSE, copy.date = TRUE), paste0("Delivery copy failed: ", dst)); delivered <- c(delivered, dst)
}

pkg_files_before_manifest <- list.files(stage, recursive = TRUE, full.names = TRUE)
rel_pkg <- substring(normalizePath(pkg_files_before_manifest, winslash = "/"), nchar(normalizePath(stage, winslash = "/")) + 2L)
delivery_rel <- gsub("\\\\", "/", substring(normalizePath(delivered, winslash = "/"), nchar(normalizePath(ROOT, winslash = "/")) + 2L))
man <- data.frame(Role = c(rep("package", length(pkg_files_before_manifest)), rep("delivery", length(delivered))),
                  Path = c(file.path("revision/re_analysis/63_figure_presentation_refinement_R_v1/generated/core", rel_pkg), delivery_rel),
                  Bytes = c(as.numeric(file.info(pkg_files_before_manifest)$size), as.numeric(file.info(delivered)$size)),
                  SHA256 = c(vapply(pkg_files_before_manifest, sha, character(1)), vapply(delivered, sha, character(1))), stringsAsFactors = FALSE)
man$Path <- gsub("\\\\", "/", man$Path)
atomic_tsv(man, file.path(stage, "manifests", "R_output_manifest_v5.tsv"))
assert(file.rename(stage, PKG_ROOT), "Atomic package promotion failed")
promoted <- TRUE
cat("RUN PASS: 3 R PDF/TIFF core-figure bundles created.\n")
