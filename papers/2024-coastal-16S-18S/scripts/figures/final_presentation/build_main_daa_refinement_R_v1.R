#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, digits = 17)

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
if (!length(file_arg)) stop("Run with Rscript", call. = FALSE)
SCRIPT <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(cowplot)
  library(digest)
  library(ragg)
  library(jsonlite)
})

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
slash <- function(x) gsub("\\\\", "/", x)

cli <- commandArgs(trailingOnly = TRUE)
modes <- c("--dry-run", "--run", "--check")
active <- intersect(cli, modes)
assert(length(active) == 1L, "Choose exactly one of --dry-run, --run, --check")
MODE <- active[[1]]

opt_value <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), cli, value = TRUE)
  if (!length(hit)) return(default)
  assert(length(hit) == 1L, paste0("Duplicate option: ", flag))
  sub(paste0("^", flag, "="), "", hit[[1]])
}

DEFAULT_ROOT <- normalizePath(file.path(dirname(SCRIPT), "..", "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(opt_value("--project-root", DEFAULT_ROOT), winslash = "/", mustWork = TRUE)
PACKAGE_ROOT <- normalizePath(file.path(dirname(SCRIPT), "..", "generated", "daa"), winslash = "/", mustWork = TRUE)
TARGET <- slash(file.path(PACKAGE_ROOT, "artifacts_v1"))
SOURCE <- slash(file.path(ROOT, "revision", "re_analysis", "29_dual_backend_figure_sources_v1", "source_data",
                          "Supplementary_Figure_S8_paired_ANCOMBC2_source_v1.tsv"))
CONTRACT <- slash(file.path(PACKAGE_ROOT, "FIGURE_CONTRACT_v1.md"))

read_char_tsv <- function(path) {
  assert(file.exists(path), paste0("Missing TSV: ", path))
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, colClasses = "character", na.strings = character())
}

need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  assert(!length(miss), paste0(label, " missing columns: ", paste(miss, collapse = ", ")))
}

relative_to_root <- function(path) {
  p <- slash(normalizePath(path, winslash = "/", mustWork = FALSE))
  prefix <- paste0(slash(ROOT), "/")
  assert(startsWith(tolower(p), tolower(prefix)) || identical(tolower(p), tolower(slash(ROOT))),
         paste0("Path outside project root: ", p))
  if (identical(tolower(p), tolower(slash(ROOT)))) "." else substring(p, nchar(prefix) + 1L)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  assert(!file.exists(path), paste0("Refusing overwrite: ", path))
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, digits = 16, na = "null")
}

prepare_source <- function() {
  x <- read_char_tsv(SOURCE)
  need_cols(x, c("Taxon", "Effect", "QValue", "Marker", "Rank", "SurfaceID", "ManuscriptStatus",
                 "ContrastLabel", "Baseline", "Target", "NSamples", "NPairs"), "canonical DA source")
  expected <- c("16S::A72_Substrate", "18S::A72_Substrate",
                "16S::A36_SeaSand_Season", "16S::A30_Water_Season",
                "18S::A36_SeaSand_Season", "18S::A30_Water_Season")
  observed <- unique(paste(x$Marker, x$SurfaceID, sep = "::"))
  assert(setequal(observed, expected), paste0("Unexpected six-surface inventory: ", paste(observed, collapse = ", ")))
  assert(all(x$Marker %in% c("16S", "18S")), "Unexpected marker")
  assert(all(x$Rank[x$Marker == "16S"] == "Genus"), "16S rank is not uniformly Genus")
  assert(all(x$Rank[x$Marker == "18S"] == "Class"), "18S rank is not uniformly Class")
  assert(all(as.numeric(x$NSamples[x$SurfaceID == "A72_Substrate"]) == 72), "A72 sample count mismatch")
  assert(all(as.numeric(x$NPairs[x$SurfaceID == "A72_Substrate"]) == 36), "A72 pair count mismatch")
  assert(all(as.numeric(x$NSamples[x$SurfaceID == "A36_SeaSand_Season"]) == 36), "A36 sample count mismatch")
  assert(all(as.numeric(x$NPairs[x$SurfaceID == "A36_SeaSand_Season"]) == 18), "A36 pair count mismatch")
  assert(all(as.numeric(x$NSamples[x$SurfaceID == "A30_Water_Season"]) == 30), "A30 sample count mismatch")
  assert(all(as.numeric(x$NPairs[x$SurfaceID == "A30_Water_Season"]) == 15), "A30 pair count mismatch")
  main <- x[x$SurfaceID == "A72_Substrate", , drop = FALSE]
  seasonal <- x[x$SurfaceID %in% c("A36_SeaSand_Season", "A30_Water_Season"), , drop = FALSE]
  assert(nrow(main) + nrow(seasonal) == nrow(x), "Six-surface partition is incomplete or duplicated")
  list(all = x, main = main, seasonal = seasonal)
}

STATUS_COLORS <- c(Main = "#B2182B", Supplementary = "#EF8A62", `Not informative` = "#BDBDBD")

plot_frame <- function(x) {
  y <- x
  y$EffectPlot <- suppressWarnings(as.numeric(y$Effect))
  y$QPlot <- as.numeric(y$QValue)
  assert(all(is.finite(y$QPlot) & y$QPlot >= 0 & y$QPlot <= 1), "Invalid QValue")
  y$PlotEligible <- is.finite(y$EffectPlot)
  assert(all(y$QPlot[!y$PlotEligible] == 1), "Undefined effects must have q = 1")
  assert(all(grepl("not_informative", y$ManuscriptStatus[!y$PlotEligible], ignore.case = TRUE)),
         "Undefined effects must be classified not informative")
  y$MinusLog10Q <- -log10(pmax(y$QPlot, .Machine$double.xmin))
  y$EvidenceStatus <- ifelse(grepl("main|three", y$ManuscriptStatus, ignore.case = TRUE), "Main",
                             ifelse(grepl("supp|two", y$ManuscriptStatus, ignore.case = TRUE),
                                    "Supplementary", "Not informative"))
  y$EvidenceStatus <- factor(y$EvidenceStatus, levels = names(STATUS_COLORS))
  y
}

theme_paper <- function(base = 8.0) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.30),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(5.0, "mm"),
      plot.tag = element_text(face = "bold", size = rel(1.25)),
      plot.title = element_text(face = "bold", hjust = 0, size = rel(0.98), margin = margin(b = 4)),
      plot.margin = margin(5, 6, 5, 6, unit = "pt")
    )
}

surface_title <- function(marker, surface) {
  rank <- ifelse(marker == "16S", "16S genus", "18S PR2 class")
  detail <- switch(surface,
                   A72_Substrate = "Water vs Sea-sand (n = 72; 36 pairs)",
                   A36_SeaSand_Season = "Sea-sand: Autumn vs Spring (n = 36; 18 pairs)",
                   A30_Water_Season = "Water: Autumn vs Spring (n = 30; 15 pairs)",
                   stop("Unsupported surface", call. = FALSE))
  paste(rank, detail, sep = "\n")
}

volcano_panel <- function(x, marker, surface, label_n = 5L) {
  d <- plot_frame(x)
  dplot <- d %>% filter(PlotEligible)
  labels <- dplot %>%
    filter(as.character(EvidenceStatus) != "Not informative") %>%
    arrange(QPlot, desc(abs(EffectPlot)), Taxon) %>%
    slice_head(n = label_n) %>%
    mutate(LabelHjust = ifelse(
      marker == "16S" & surface == "A72_Substrate" & Taxon == "Zeaxanthinibacter",
      0,
      0.5
    ))
  label_nudge_x <- ifelse(
    marker == "16S" & surface == "A72_Substrate" & labels$Taxon == "Zeaxanthinibacter",
    0.35,
    0
  )
  ggplot(dplot, aes(EffectPlot, MinusLog10Q)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey35", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey72", linewidth = 0.30) +
    geom_point(aes(fill = EvidenceStatus), shape = 21, color = "black", stroke = 0.16,
               size = 1.55, alpha = 0.90) +
    ggrepel::geom_text_repel(data = labels, aes(label = Taxon, hjust = LabelHjust), size = 2.05,
                             max.overlaps = Inf, min.segment.length = 0,
                             seed = 20260814, box.padding = 0.24, point.padding = 0.12,
                             nudge_x = label_nudge_x,
                             segment.color = "grey35", segment.size = 0.22) +
    scale_fill_manual(values = STATUS_COLORS, breaks = names(STATUS_COLORS), drop = FALSE) +
    labs(x = "ANCOM-BC2 log fold change", y = expression(bold(-log[10](BH~q))),
         title = surface_title(marker, surface), fill = "Evidence") +
    coord_cartesian(clip = "off") +
    theme_paper() + theme(legend.position = "bottom")
}

make_main <- function(x) {
  p16 <- volcano_panel(x[x$Marker == "16S", , drop = FALSE], "16S", "A72_Substrate")
  p18 <- volcano_panel(x[x$Marker == "18S", , drop = FALSE], "18S", "A72_Substrate")
  (p16 | p18) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
}

make_seasonal <- function(x) {
  order <- data.frame(
    Marker = c("16S", "16S", "18S", "18S"),
    SurfaceID = c("A36_SeaSand_Season", "A30_Water_Season", "A36_SeaSand_Season", "A30_Water_Season"),
    stringsAsFactors = FALSE
  )
  plots <- lapply(seq_len(nrow(order)), function(i) {
    z <- x[x$Marker == order$Marker[[i]] & x$SurfaceID == order$SurfaceID[[i]], , drop = FALSE]
    assert(nrow(z) > 0, paste0("Missing seasonal panel: ", order$Marker[[i]], " ", order$SurfaceID[[i]]))
    volcano_panel(z, order$Marker[[i]], order$SurfaceID[[i]]) + labs(tag = LETTERS[[i]])
  })
  shared_legend <- cowplot::get_legend(
    plots[[1]] + theme(legend.position = "bottom", legend.box = "horizontal")
  )
  plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
  body <- wrap_plots(plots, ncol = 2)
  body / patchwork::wrap_elements(full = shared_legend) + plot_layout(heights = c(1, 0.075))
}

save_plot <- function(plot, stem, width_mm, height_mm, root) {
  fig_dir <- file.path(root, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf <- file.path(fig_dir, paste0(stem, ".pdf"))
  tif <- file.path(fig_dir, paste0(stem, ".tiff"))
  assert(!file.exists(pdf) && !file.exists(tif), paste0("Refusing overwrite: ", stem))
  grDevices::cairo_pdf(pdf, width = width_mm / 25.4, height = height_mm / 25.4,
                       family = "Arial", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(tif, width = width_mm, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white", bitsize = 8)
  print(plot)
  grDevices::dev.off()
  assert(all(file.info(c(pdf, tif))$size > 0), paste0("Empty figure output: ", stem))
  c(pdf, tif)
}

signature_ok <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw <- readBin(con, "raw", n = 8L)
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) return(rawToChar(raw[seq_len(5)]) == "%PDF-")
  if (grepl("\\.tiff?$", path, ignore.case = TRUE)) {
    return(identical(as.integer(raw[1:4]), c(73L, 73L, 42L, 0L)) ||
             identical(as.integer(raw[1:4]), c(77L, 77L, 0L, 42L)))
  }
  FALSE
}

surface_ledger <- function(x) {
  y <- plot_frame(x)
  y %>%
    group_by(Marker, Rank, SurfaceID, ContrastLabel, Baseline, Target, NSamples, NPairs) %>%
    summarise(
      TaxaTested = n(),
      TaxaWithFiniteEffect = sum(PlotEligible),
      TaxaWithUndefinedEffect = sum(!PlotEligible),
      MainCandidates = sum(EvidenceStatus == "Main"),
      SupplementaryCandidates = sum(EvidenceStatus == "Supplementary"),
      NotInformative = sum(EvidenceStatus == "Not informative"),
      MinQ = min(QPlot),
      .groups = "drop"
    )
}

manifest_for <- function(paths, roles, output_paths = NULL) {
  assert(length(paths) == length(roles), "Manifest length mismatch")
  data.frame(
    Role = roles,
    Path = if (is.null(output_paths)) vapply(paths, relative_to_root, character(1)) else output_paths,
    Bytes = as.numeric(file.info(paths)$size),
    SHA256 = vapply(paths, sha256, character(1)),
    stringsAsFactors = FALSE
  )
}

expected_tiff_pixels <- function(width_mm, height_mm) round(c(width_mm, height_mm) / 25.4 * 600)

check_artifacts <- function(target, src) {
  assert(dir.exists(target), paste0("Missing artifact root: ", target))
  input_manifest_path <- file.path(target, "provenance", "input_sha256_manifest_v1.tsv")
  output_manifest_path <- file.path(target, "provenance", "output_sha256_manifest_v1.tsv")
  assert(file.exists(input_manifest_path) && file.exists(output_manifest_path), "Missing manifests")
  im <- read.delim(input_manifest_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  om <- read.delim(output_manifest_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(im))) {
    p <- file.path(ROOT, im$Path[[i]])
    assert(file.exists(p), paste0("Missing recorded input: ", im$Path[[i]]))
    assert(as.numeric(file.info(p)$size) == as.numeric(im$Bytes[[i]]), paste0("Input byte drift: ", im$Path[[i]]))
    assert(sha256(p) == im$SHA256[[i]], paste0("Input hash drift: ", im$Path[[i]]))
  }
  for (i in seq_len(nrow(om))) {
    p <- file.path(ROOT, om$Path[[i]])
    assert(file.exists(p), paste0("Missing recorded output: ", om$Path[[i]]))
    assert(as.numeric(file.info(p)$size) == as.numeric(om$Bytes[[i]]), paste0("Output byte drift: ", om$Path[[i]]))
    assert(sha256(p) == om$SHA256[[i]], paste0("Output hash drift: ", om$Path[[i]]))
  }
  main_path <- file.path(target, "source_data", "Figure_4_paired_sample_type_DA_source_v1.tsv")
  season_path <- file.path(target, "source_data", "Supplementary_Figure_seasonal_DA_source_v1.tsv")
  main_out <- read_char_tsv(main_path)
  season_out <- read_char_tsv(season_path)
  assert(identical(names(main_out), names(src$main)) &&
           identical(unname(as.matrix(main_out)), unname(as.matrix(src$main))),
         "Main source-data cell mismatch")
  assert(identical(names(season_out), names(src$seasonal)) &&
           identical(unname(as.matrix(season_out)), unname(as.matrix(src$seasonal))),
         "Seasonal source-data cell mismatch")
  assert(nrow(main_out) + nrow(season_out) == nrow(src$all), "Source partition row-count mismatch")
  assert(length(unique(paste(rbind(main_out, season_out)$Marker, rbind(main_out, season_out)$SurfaceID, sep = "::"))) == 6L,
         "Six-surface completeness failed")
  figs <- c(
    file.path(target, "figures", "Figure_4_paired_sample_type_differential_abundance_v1.pdf"),
    file.path(target, "figures", "Figure_4_paired_sample_type_differential_abundance_v1.tiff"),
    file.path(target, "figures", "Supplementary_Figure_seasonal_differential_abundance_v1.pdf"),
    file.path(target, "figures", "Supplementary_Figure_seasonal_differential_abundance_v1.tiff")
  )
  assert(all(file.exists(figs)) && all(file.info(figs)$size > 0), "Missing/empty figure")
  assert(all(vapply(figs, signature_ok, logical(1))), "Figure signature failure")
  if (requireNamespace("pdftools", quietly = TRUE)) {
    assert(pdftools::pdf_info(figs[[1]])$pages == 1L && pdftools::pdf_info(figs[[3]])$pages == 1L,
           "PDF page-count failure")
  }
  if (requireNamespace("magick", quietly = TRUE)) {
    i1 <- magick::image_info(magick::image_read(figs[[2]]))[1, ]
    i2 <- magick::image_info(magick::image_read(figs[[4]]))[1, ]
    e1 <- expected_tiff_pixels(183, 105)
    e2 <- expected_tiff_pixels(183, 180)
    assert(abs(i1$width - e1[[1]]) <= 1 && abs(i1$height - e1[[2]]) <= 1, "Main TIFF pixel dimensions mismatch")
    assert(abs(i2$width - e2[[1]]) <= 1 && abs(i2$height - e2[[2]]) <= 1, "Seasonal TIFF pixel dimensions mismatch")
  }
  TRUE
}

src <- prepare_source()

if (MODE == "--dry-run") {
  assert(!dir.exists(TARGET), paste0("Target already exists; use --check: ", TARGET))
  led <- surface_ledger(src$all)
  cat("DRY_RUN_OK\n")
  cat("canonical_rows=", nrow(src$all), " main_rows=", nrow(src$main), " seasonal_rows=", nrow(src$seasonal), "\n", sep = "")
  cat("surfaces=", nrow(led), " main_panels=2 seasonal_panels=4\n", sep = "")
  quit(save = "no", status = 0)
}

if (MODE == "--check") {
  check_artifacts(TARGET, src)
  cat("CHECK_OK\n")
  quit(save = "no", status = 0)
}

assert(MODE == "--run", "Internal mode error")
assert(!dir.exists(TARGET), paste0("Refusing overwrite; target exists: ", TARGET))
stage <- paste0(TARGET, ".staging-", Sys.getpid())
assert(!dir.exists(stage), paste0("Staging path exists: ", stage))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
promoted <- FALSE
on.exit(if (!promoted && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

for (d in c("figures", "source_data", "tables", "qc", "provenance", "logs", "captions"))
  dir.create(file.path(stage, d), recursive = TRUE, showWarnings = FALSE)

write_tsv(src$main, file.path(stage, "source_data", "Figure_4_paired_sample_type_DA_source_v1.tsv"))
write_tsv(src$seasonal, file.path(stage, "source_data", "Supplementary_Figure_seasonal_DA_source_v1.tsv"))
ledger <- surface_ledger(src$all)
write_tsv(ledger, file.path(stage, "tables", "differential_abundance_surface_ledger_v1.tsv"))

save_plot(make_main(src$main), "Figure_4_paired_sample_type_differential_abundance_v1", 183, 105, stage)
save_plot(make_seasonal(src$seasonal), "Supplementary_Figure_seasonal_differential_abundance_v1", 183, 180, stage)

caption <- c(
  "# Differential-abundance figure caption facts v1",
  "",
  "## Main Figure 4",
  "",
  "Paired differential-abundance results for (a) 16S genera and (b) PR2-primary 18S classes comparing Water with Sea-sand. Each marker surface contains 72 samples arranged as 36 contemporaneous, co-located Sea-sand/Water pairs. Points show ANCOM-BC2 log fold changes and within-surface BH-adjusted q values; colors retain the accepted multi-method and leave-one-Site-out evidence classification.",
  "",
  "## Supplementary seasonal figure",
  "",
  "Within-substrate seasonal differential-abundance results for 16S genera and PR2-primary 18S classes. Sea-sand panels contain 36 samples arranged as 18 fixed-position Spring/Autumn pairs; Water panels contain 30 samples arranged as 15 fixed-position Spring/Autumn pairs. The four seasonal surfaces are shown together to avoid selective emphasis on one substrate. Seasonal surfaces were analyzed separately and were not a formal test that the seasonal effect differed between substrates.",
  "",
  "## Shared qualifier",
  "",
  "Differential-abundance results are associations within the prespecified paired/blocked analysis sets. They do not establish causality or universal biomarkers. Taxa with an undefined ANCOM-BC2 coefficient (all classified as not informative with q = 1) remain in the source-data tables but cannot be positioned on the log-fold-change axis and are therefore omitted from the plotted points."
)
writeLines(caption, file.path(stage, "captions", "differential_abundance_caption_facts_v1.md"), useBytes = TRUE)

status <- data.frame(
  Artifact = c("Figure_4_paired_sample_type_differential_abundance_v1",
               "Supplementary_Figure_seasonal_differential_abundance_v1"),
  AnalysisChanged = c(FALSE, FALSE),
  SourceRows = c(nrow(src$main), nrow(src$seasonal)),
  Panels = c(2L, 4L),
  ManuscriptRole = c("main", "supplementary"),
  GateStatus = c("accepted_candidate_pending_independent_visual_QA", "accepted_candidate_pending_independent_visual_QA"),
  Reason = c("Only the two central paired SampleType surfaces are retained in the main figure.",
             "All four within-substrate seasonal surfaces are presented together with equal grammar."),
  stringsAsFactors = FALSE
)
write_tsv(status, file.path(stage, "qc", "figure_status_gate_v1.tsv"))

run_config <- list(
  version = "1.0.0",
  backend = "R",
  analysis_refit = FALSE,
  canonical_source = relative_to_root(SOURCE),
  main_surfaces = c("16S::A72_Substrate", "18S::A72_Substrate"),
  supplementary_surfaces = c("16S::A36_SeaSand_Season", "16S::A30_Water_Season",
                             "18S::A36_SeaSand_Season", "18S::A30_Water_Season"),
  evidence_colors = as.list(STATUS_COLORS),
  seed_for_label_placement = 20260814,
  output = list(main_mm = c(183, 105), supplementary_mm = c(183, 180), tiff_dpi = 600, tiff_compression = "LZW")
)
write_json(run_config, file.path(stage, "provenance", "run_config_v1.json"))
capture.output(sessionInfo(), file = file.path(stage, "provenance", "R_session_info_v1.txt"))

readme <- c(
  "# Balanced differential-abundance presentation v1",
  "",
  "This immutable package reassembles accepted ANCOM-BC2 source values; it does not refit any differential-abundance model.",
  "",
  "- Main Figure 4 contains the two paired SampleType surfaces.",
  "- The supplementary seasonal figure contains all four Sea-sand and Water seasonal surfaces.",
  "- Source tables preserve the canonical rows and columns exactly.",
  "- PDF and 600-dpi LZW TIFF files were produced in R.",
  "- The full visual acceptance decision is recorded separately after rendered inspection."
)
writeLines(readme, file.path(stage, "README_v1.md"), useBytes = TRUE)
writeLines(c("2026-08-14", "R-only presentation build", "No statistical model refitted"),
           file.path(stage, "logs", "build_v1.log"), useBytes = TRUE)

input_paths <- c(SOURCE, SCRIPT, CONTRACT)
input_roles <- c("canonical_DA_source", "figure_builder", "figure_contract")
write_tsv(manifest_for(input_paths, input_roles), file.path(stage, "provenance", "input_sha256_manifest_v1.tsv"))

qa <- list(
  status = "PASS_PENDING_VISUAL_INSPECTION",
  canonical_rows = nrow(src$all),
  main_rows = nrow(src$main),
  seasonal_rows = nrow(src$seasonal),
  surfaces = nrow(ledger),
  main_panels = 2,
  supplementary_panels = 4,
  source_partition_complete = nrow(src$main) + nrow(src$seasonal) == nrow(src$all),
  analysis_refit = FALSE,
  figure_signatures_valid = TRUE
)
write_json(qa, file.path(stage, "qc", "figure_numeric_file_QA_v1.json"))

all_stage_files <- list.files(stage, recursive = TRUE, full.names = TRUE, all.files = FALSE)
manifest_path <- file.path(stage, "provenance", "output_sha256_manifest_v1.tsv")
all_stage_files <- setdiff(all_stage_files, manifest_path)
virtual_paths <- vapply(all_stage_files, function(p) {
  rel <- substring(slash(p), nchar(slash(stage)) + 2L)
  relative_to_root(file.path(TARGET, rel))
}, character(1))
ord <- order(virtual_paths)
write_tsv(manifest_for(all_stage_files[ord], rep("output", length(ord)), virtual_paths[ord]), manifest_path)

assert(file.rename(stage, TARGET), "Atomic promotion failed")
promoted <- TRUE
check_artifacts(TARGET, src)
cat("RUN_OK\n")
