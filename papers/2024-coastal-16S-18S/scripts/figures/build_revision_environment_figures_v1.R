#!/usr/bin/env Rscript

SCRIPT_VERSION <- "1.0.0"
PACKAGE_VERSION <- "v1"
FIG5_WIDTH_MM <- 183
FIG5_HEIGHT_MM <- 174
FIG6_WIDTH_MM <- 183
FIG6_HEIGHT_MM <- 100
TIFF_DPI <- 600
PREVIEW_DPI <- 300

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
assert_true <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) stopf("Required R package is missing: %s", pkg)
}
invisible(lapply(c("ggplot2", "ggrepel", "patchwork", "cowplot", "scales", "digest", "jsonlite", "pdftools"), require_package))

parse_args <- function(args) {
  out <- list(mode = NULL, project_root = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--run", "--dry-run", "--check")) {
      if (!is.null(out$mode)) stopf("Specify exactly one mode: --run, --dry-run, or --check")
      out$mode <- sub("^--", "", arg)
      i <- i + 1L
    } else if (arg == "--project-root") {
      if (i == length(args)) stopf("--project-root requires a value")
      out$project_root <- normalizePath(args[[i + 1L]], winslash = "/", mustWork = TRUE)
      i <- i + 2L
    } else if (arg %in% c("-h", "--help")) {
      cat(paste(
        "Usage:",
        "  Rscript build_revision_environment_figures_v1.R --dry-run [--project-root PATH]",
        "  Rscript build_revision_environment_figures_v1.R --run [--project-root PATH]",
        "  Rscript build_revision_environment_figures_v1.R --check [--project-root PATH]",
        sep = "\n"
      ), "\n")
      quit(status = 0L)
    } else {
      stopf("Unknown argument: %s", arg)
    }
  }
  if (is.null(out$mode)) stopf("Specify exactly one mode: --run, --dry-run, or --check")
  out
}

path_join <- function(...) normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
read_tsv <- function(path) {
  assert_true(file.exists(path), sprintf("Missing input: %s", path))
  utils::read.delim(
    path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = "", na.strings = c("NA")
  )
}
sha256_file <- function(path) {
  assert_true(file.exists(path), sprintf("Cannot hash missing file: %s", path))
  tolower(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}
file_bytes <- function(path) as.numeric(file.info(path)$size)

write_tsv_new <- function(x, path) {
  if (file.exists(path)) stopf("Refusing to overwrite: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
    na = "NA", fileEncoding = "UTF-8", eol = "\n"
  )
}
write_text_new <- function(lines, path) {
  if (file.exists(path)) stopf("Refusing to overwrite: %s", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}
copy_file_new <- function(from, to) {
  assert_true(file.exists(from), sprintf("Cannot copy missing file: %s", from))
  if (file.exists(to)) stopf("Refusing to overwrite: %s", to)
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  assert_true(file.copy(from, to, overwrite = FALSE, copy.date = TRUE), sprintf("Copy failed: %s", from))
}

format_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(sub("^0", "", sprintf("%.5f", x)))
  sub("^0", "", sprintf("%.4f", x))
}
format_num <- function(x, digits = 3L) sprintf(paste0("%.", digits, "f"), x)

resolve_paths <- function(root) {
  src4 <- path_join(root, "revision/re_analysis/10_figure_source_data_v4/source_data")
  env8 <- path_join(root, "revision/re_analysis/08_environment_exploratory_v1/results")
  package_root <- path_join(root, "revision/re_analysis/16_revision_environment_figures_v1")
  delivery_dir <- path_join(root, "revision/rev_ms/figures")
  script_path <- path_join(root, "revision/re_analysis/scripts/build_revision_environment_figures_v1.R")
  contract_source <- path_join(root, "revision/re_analysis/figure_contracts/environment_figure_contract_v1.md")
  list(
    root = root,
    package_root = package_root,
    delivery_dir = delivery_dir,
    script = script_path,
    contract_source = contract_source,
    palette = path_join(root, "revision/re_analysis/00_design_v4/palette_contract_v4.tsv"),
    current_ds4_pdf = path_join(root, "revision/re_analysis/10_figure_source_data_v4/figures/Fig_DS4_spring_seasand_environment_v4.pdf"),
    current_ds4_png = path_join(root, "revision/re_analysis/10_figure_source_data_v4/figures/Fig_DS4_spring_seasand_environment_v4.png"),
    ds4 = c(
      pca_sites = path_join(src4, "DS4_environment_pca_sites_source_data_v1.tsv"),
      pca_loadings = path_join(src4, "DS4_environment_pca_loadings_source_data_v1.tsv"),
      rda_sites = path_join(src4, "DS4_environment_rda_sites_source_data_v1.tsv"),
      rda_biplot = path_join(src4, "DS4_environment_rda_biplot_source_data_v1.tsv"),
      rda_tests = path_join(src4, "DS4_environment_rda_tests_source_data_v1.tsv"),
      mantel = path_join(src4, "DS4_environment_mantel_summary_source_data_v1.tsv"),
      mantel_loso = path_join(src4, "DS4_environment_mantel_loso_source_data_v1.tsv")
    ),
    canonical = c(
      pca_sites = path_join(env8, "environment_pca_site_scores_v1.tsv"),
      pca_loadings = path_join(env8, "environment_pca_variable_loadings_v1.tsv"),
      rda_sites = path_join(env8, "environment_rda_site_scores_v1.tsv"),
      rda_biplot = path_join(env8, "environment_rda_biplot_scores_v1.tsv"),
      rda_tests = path_join(env8, "environment_rda_model_tests_v1.tsv"),
      mantel = path_join(env8, "environment_mantel_summary_v1.tsv"),
      mantel_loso = path_join(env8, "environment_mantel_loso_v1.tsv")
    ),
    deliverables = c(
      figure5_pdf = path_join(delivery_dir, "Figure_5_16S_environment_associations_v1.pdf"),
      figure5_tiff = path_join(delivery_dir, "Figure_5_16S_environment_associations_v1.tiff"),
      figure6_pdf = path_join(delivery_dir, "Figure_6_18S_environment_associations_v1.pdf"),
      figure6_tiff = path_join(delivery_dir, "Figure_6_18S_environment_associations_v1.tiff")
    )
  )
}

require_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing)) stopf("%s missing columns: %s", label, paste(missing, collapse = ", "))
}

align_and_compare <- function(canonical, ds4, keys, label, tolerance = 1e-12) {
  require_columns(canonical, keys, paste(label, "canonical"))
  require_columns(ds4, keys, paste(label, "DS4"))
  make_key <- function(x) do.call(paste, c(x[keys], sep = "\r"))
  kc <- make_key(canonical)
  kd <- make_key(ds4)
  assert_true(!anyDuplicated(kc), sprintf("%s canonical keys are duplicated", label))
  assert_true(!anyDuplicated(kd), sprintf("%s DS4 keys are duplicated", label))
  assert_true(setequal(kc, kd), sprintf("%s canonical and DS4 key sets differ", label))
  d <- ds4[match(kc, kd), , drop = FALSE]
  shared <- intersect(names(canonical), names(d))
  for (nm in shared) {
    a <- canonical[[nm]]
    b <- d[[nm]]
    if (is.numeric(a) || is.numeric(b)) {
      aa <- suppressWarnings(as.numeric(a)); bb <- suppressWarnings(as.numeric(b))
      same_na <- is.na(aa) == is.na(bb)
      delta <- abs(aa - bb)
      delta[is.na(delta)] <- 0
      assert_true(all(same_na) && all(delta <= tolerance), sprintf("%s numeric mismatch in %s", label, nm))
    } else {
      aa <- ifelse(is.na(a), "<NA>", as.character(a))
      bb <- ifelse(is.na(b), "<NA>", as.character(b))
      assert_true(identical(aa, bb), sprintf("%s text mismatch in %s", label, nm))
    }
  }
  invisible(TRUE)
}

load_and_validate <- function(paths) {
  for (p in c(paths$ds4, paths$canonical, paths$palette, paths$current_ds4_pdf, paths$current_ds4_png, paths$script, paths$contract_source)) {
    assert_true(file.exists(p), sprintf("Missing required input: %s", p))
  }
  ds4 <- lapply(paths$ds4, read_tsv)
  canonical <- lapply(paths$canonical, read_tsv)
  key_map <- list(
    pca_sites = "SiteCodeOriginal",
    pca_loadings = "Variable",
    rda_sites = c("Marker", "SiteCodeOriginal"),
    rda_biplot = c("Marker", "Label"),
    rda_tests = c("Marker", "TestType", "Term"),
    mantel = "Marker",
    mantel_loso = c("Marker", "LeftOutSiteCode")
  )
  for (nm in names(ds4)) align_and_compare(canonical[[nm]], ds4[[nm]], key_map[[nm]], nm)

  expected_sites <- c("DD", "ULD", "GS", "PH", "PS", "SH", "SS")
  expected_variables <- c("pH", "NaCl", "T-N", "T-S", "Ca", "Mg", "K", "Na", "T-P", "water")
  assert_true(nrow(ds4$pca_sites) == 7L && setequal(ds4$pca_sites$SiteCodeOriginal, expected_sites), "PCA must contain exactly seven Sites")
  assert_true(nrow(ds4$pca_loadings) == 10L && setequal(ds4$pca_loadings$Variable, expected_variables), "PCA must contain exactly the ten predeclared variables")
  assert_true(!any(toupper(ds4$pca_loadings$Variable) == "TOC"), "TOC must be excluded from the primary PCA")
  assert_true(all(table(ds4$rda_sites$Marker) == 7L), "RDA must contain seven Site scores per marker")
  assert_true(all(table(ds4$rda_biplot$Marker) == 2L), "RDA must contain EnvPC1 and EnvPC2 arrows per marker")
  assert_true(nrow(ds4$mantel) == 2L && setequal(ds4$mantel$Marker, c("16S", "18S")), "Mantel summary must contain two marker rows")
  assert_true(all(table(ds4$mantel_loso$Marker) == 7L), "Mantel LOSO must contain seven rows per marker")

  scope_tables <- ds4
  for (nm in names(scope_tables)) {
    x <- scope_tables[[nm]]
    for (col in c("SetID", "SampleType", "InferenceScope")) {
      require_columns(x, col, nm)
    }
    season_col <- if ("Season" %in% names(x)) "Season" else "SeasonDisplay"
    assert_true(all(x$SetID == "A7_SITE_CENTROIDS"), sprintf("%s contains non-A7 rows", nm))
    assert_true(all(x$SampleType == "Sea-sand"), sprintf("%s contains non-Sea-sand rows", nm))
    assert_true(all(x[[season_col]] == "Spring"), sprintf("%s contains non-Spring rows", nm))
    assert_true(all(x$InferenceScope == "exploratory_site_level_n7"), sprintf("%s scope differs", nm))
  }

  overall <- ds4$rda_tests[ds4$rda_tests$TestType == "overall" & ds4$rda_tests$Term == "Model", , drop = FALSE]
  overall <- overall[match(c("16S", "18S"), overall$Marker), , drop = FALSE]
  mantel <- ds4$mantel[match(c("16S", "18S"), ds4$mantel$Marker), , drop = FALSE]
  expected <- data.frame(
    Marker = c("16S", "18S"),
    F = c(3.30507586323131, 1.64263334407708),
    RdaP = c(0.00694444444444444, 0.0527777777777778),
    R2 = c(0.623002563665168, 0.450946661087369),
    AdjustedR2 = c(0.434503845497752, 0.176419991631054),
    MantelRho = c(0.644155844155844, 0.666233766233766),
    MantelP = c(0.00615079365079365, 0.000992063492063492),
    MantelQ = c(0.00615079365079365, 0.00198412698412698),
    stringsAsFactors = FALSE
  )
  observed <- data.frame(
    Marker = overall$Marker, F = overall$F, RdaP = overall$P,
    R2 = overall$R2, AdjustedR2 = overall$AdjustedR2,
    MantelRho = mantel$Statistic, MantelP = mantel$P,
    MantelQ = mantel$PAdjustedBHPrimaryAcrossMarkers,
    stringsAsFactors = FALSE
  )
  numeric_cols <- setdiff(names(expected), "Marker")
  for (nm in numeric_cols) assert_true(max(abs(observed[[nm]] - expected[[nm]])) <= 1e-12, sprintf("Frozen statistic mismatch: %s", nm))

  loso_ranges <- do.call(rbind, lapply(c("16S", "18S"), function(marker) {
    z <- ds4$mantel_loso[ds4$mantel_loso$Marker == marker, , drop = FALSE]
    data.frame(Marker = marker, LOSOMin = min(z$Statistic), LOSOMax = max(z$Statistic), stringsAsFactors = FALSE)
  }))
  expected_min <- c(0.446428571428571, 0.521428571428571)
  expected_max <- c(0.810714285714286, 0.732142857142857)
  assert_true(max(abs(loso_ranges$LOSOMin - expected_min)) <= 1e-12, "LOSO minima differ from frozen values")
  assert_true(max(abs(loso_ranges$LOSOMax - expected_max)) <= 1e-12, "LOSO maxima differ from frozen values")
  assert_true(all(ds4$mantel_loso$SignMatchesFull %in% c(TRUE, "TRUE", 1, "1")), "A LOSO estimate changes sign")

  palette <- read_tsv(paths$palette)
  p <- palette[palette$Aesthetic == "color" & palette$Variable == "GeologicalHabitat", , drop = FALSE]
  p <- p[match(c("Ocean", "Tidal", "Island"), p$Level), , drop = FALSE]
  assert_true(!anyNA(p$Value), "Frozen geological palette is incomplete")
  colors <- stats::setNames(p$Value, p$Level)
  expected_colors <- c(Ocean = "#0072B2", Tidal = "#E69F00", Island = "#009E73")
  assert_true(identical(colors, expected_colors), "Frozen geological palette differs")

  list(ds4 = ds4, overall = overall, mantel = mantel, loso_ranges = loso_ranges, colors = colors)
}

theme_environment <- function(base_size = 8.6) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.6, hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.6, hjust = 0, colour = "#222222", lineheight = 1.05, margin = ggplot2::margin(b = 3)),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 1.2),
      plot.tag.position = c(0.01, 0.99),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.8, colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "#E6E6E6", linewidth = 0.28),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", linewidth = 0.45),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 0.8),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.margin = ggplot2::margin(t = 1, r = 1, b = 0, l = 1),
      plot.margin = ggplot2::margin(5, 6, 5, 6)
    )
}

habitat_scales <- function(colors) list(
  ggplot2::scale_fill_manual(values = colors, breaks = names(colors), drop = FALSE, name = "Geological context"),
  ggplot2::scale_colour_manual(values = colors, breaks = names(colors), drop = FALSE, name = "Geological context")
)

build_pca <- function(d, colors) {
  sites <- d$pca_sites
  loadings <- d$pca_loadings
  sites$GeologicalHabitat <- factor(sites$GeologicalHabitat, levels = names(colors))
  p1 <- 100 * unique(sites$PC1ProportionVariance)
  p2 <- 100 * unique(sites$PC2ProportionVariance)
  assert_true(length(p1) == 1L && length(p2) == 1L, "PCA variance labels are not unique")
  set.seed(20260811)
  ggplot2::ggplot(sites, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.30) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8B8B8", linewidth = 0.30) +
    ggplot2::geom_segment(
      data = loadings, ggplot2::aes(x = 0, y = 0, xend = PC1Plot, yend = PC2Plot),
      inherit.aes = FALSE, colour = "black", linewidth = 0.42,
      arrow = grid::arrow(length = grid::unit(1.55, "mm"), type = "closed")
    ) +
    ggrepel::geom_text_repel(
      data = loadings, ggplot2::aes(x = PC1Plot, y = PC2Plot, label = Variable),
      inherit.aes = FALSE, colour = "black", size = 2.45, family = "sans",
      seed = 101, box.padding = 0.18, point.padding = 0.08,
      min.segment.length = 0, segment.colour = "#555555", segment.size = 0.22,
      max.overlaps = Inf, max.time = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = GeologicalHabitat), shape = 21, colour = "black",
      size = 2.55, stroke = 0.48
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Site), colour = "black", size = 2.55,
      family = "sans", fontface = "bold", seed = 102,
      box.padding = 0.28, point.padding = 0.28, min.segment.length = 0,
      segment.colour = "#777777", segment.size = 0.25, max.overlaps = Inf,
      max.time = 2, show.legend = FALSE
    ) +
    habitat_scales(colors) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.20)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.20)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = "Standardized environmental PCA",
      x = sprintf("Environmental PC1 (%.1f%%)", p1),
      y = sprintf("Environmental PC2 (%.1f%%)", p2), tag = "a"
    ) +
    theme_environment() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.7)), colour = "none") +
    ggplot2::theme(aspect.ratio = 1)
}

build_rda <- function(d, colors, marker, tag) {
  sites <- d$rda_sites[d$rda_sites$Marker == marker, , drop = FALSE]
  arrows <- d$rda_biplot[d$rda_biplot$Marker == marker & as.character(d$rda_biplot$ArrowSelected) %in% c("TRUE", "1"), , drop = FALSE]
  tests <- d$rda_tests[d$rda_tests$Marker == marker & d$rda_tests$TestType == "overall" & d$rda_tests$Term == "Model", , drop = FALSE]
  assert_true(nrow(sites) == 7L && nrow(arrows) == 2L && nrow(tests) == 1L, sprintf("Invalid %s RDA plot grain", marker))
  sites$GeologicalHabitat <- factor(sites$GeologicalHabitat, levels = names(colors))
  feature <- unique(sites$FeatureLevel)
  axis1 <- unique(sites$RDA1PercentConstrained)
  axis2 <- unique(sites$RDA2PercentConstrained)
  assert_true(length(feature) == 1L && length(axis1) == 1L && length(axis2) == 1L, sprintf("%s RDA labels are not unique", marker))
  annotation <- sprintf(
    "F = %.3f; permutation P = %s\nR2 = %.3f; adjusted R2 = %.3f; n = 7",
    tests$F, format_p(tests$P), tests$R2, tests$AdjustedR2
  )
  set.seed(if (marker == "16S") 202616L else 202618L)
  ggplot2::ggplot(sites, ggplot2::aes(RDA1, RDA2)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B8B8B8", linewidth = 0.30) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8B8B8", linewidth = 0.30) +
    ggplot2::geom_segment(
      data = arrows, ggplot2::aes(x = 0, y = 0, xend = RDA1Plot, yend = RDA2Plot),
      inherit.aes = FALSE, colour = "black", linewidth = 0.45,
      arrow = grid::arrow(length = grid::unit(1.55, "mm"), type = "closed")
    ) +
    ggrepel::geom_text_repel(
      data = arrows, ggplot2::aes(x = RDA1Plot, y = RDA2Plot, label = Label),
      inherit.aes = FALSE, colour = "black", size = 2.55, family = "sans",
      seed = 201, box.padding = 0.15, point.padding = 0.05,
      min.segment.length = 0, segment.colour = "#555555", segment.size = 0.22,
      max.overlaps = Inf, max.time = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = GeologicalHabitat), shape = 21, colour = "black",
      size = 2.65, stroke = 0.48
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Site), colour = "black", size = 2.55,
      family = "sans", fontface = "bold", seed = 202,
      box.padding = 0.28, point.padding = 0.28, min.segment.length = 0,
      segment.colour = "#777777", segment.size = 0.25, max.overlaps = Inf,
      max.time = 2, show.legend = FALSE
    ) +
    ggplot2::annotate(
      "label", x = -Inf, y = Inf, label = annotation,
      hjust = -0.03, vjust = 1.08, family = "sans", size = 2.35,
      linewidth = 0, fill = scales::alpha("white", 0.88), colour = "#222222"
    ) +
    habitat_scales(colors) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.23)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.23)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = sprintf("%s Hellinger RDA (%s)", marker, tolower(feature)),
      x = sprintf("RDA1 (%.1f%% constrained)", axis1),
      y = sprintf("RDA2 (%.1f%% constrained)", axis2), tag = tag
    ) +
    theme_environment() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.7)), colour = "none") +
    ggplot2::theme(aspect.ratio = 1)
}

build_mantel <- function(d, colors, marker, tag, aspect_ratio = NULL) {
  full <- d$mantel[d$mantel$Marker == marker, , drop = FALSE]
  loso <- d$mantel_loso[d$mantel_loso$Marker == marker, , drop = FALSE]
  assert_true(nrow(full) == 1L && nrow(loso) == 7L, sprintf("Invalid %s Mantel grain", marker))
  site_order <- c("Dokdo", "Ulleungdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
  loso <- loso[match(site_order, loso$LeftOutSite), , drop = FALSE]
  assert_true(!anyNA(loso$LeftOutSite), sprintf("%s Mantel Site ordering failed", marker))
  loso$RowLabel <- paste("Omit", loso$LeftOutSite, "(n = 6)")
  y_levels <- rev(c("Full data (n = 7)", loso$RowLabel))
  loso$RowLabel <- factor(loso$RowLabel, levels = y_levels)
  loso$GeologicalHabitat <- factor(loso$GeologicalHabitat, levels = names(colors))
  full_row <- data.frame(RowLabel = factor("Full data (n = 7)", levels = y_levels), Statistic = full$Statistic)
  range_min <- min(loso$Statistic)
  range_max <- max(loso$Statistic)
  annotation <- sprintf(
    "Full rho = %.3f; permutation P = %s; BH q = %s\nLOSO rho range = %.3f-%.3f",
    full$Statistic, format_p(full$P), format_p(full$PAdjustedBHPrimaryAcrossMarkers), range_min, range_max
  )
  x_limits <- c(0.00, 0.90)
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = full$Statistic, linetype = "dashed", colour = "#777777", linewidth = 0.38) +
    ggplot2::geom_segment(
      data = loso,
      ggplot2::aes(x = full$Statistic, xend = Statistic, y = RowLabel, yend = RowLabel),
      inherit.aes = FALSE, colour = "#B5B5B5", linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = loso, ggplot2::aes(x = Statistic, y = RowLabel, fill = GeologicalHabitat),
      shape = 21, colour = "black", size = 2.5, stroke = 0.45
    ) +
    ggplot2::geom_point(
      data = full_row, ggplot2::aes(x = Statistic, y = RowLabel),
      shape = 23, fill = "black", colour = "black", size = 2.9
    ) +
    habitat_scales(colors) +
    ggplot2::scale_x_continuous(limits = x_limits, breaks = seq(0.0, 0.9, by = 0.1), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s Mantel association and LOSO influence", marker),
      subtitle = annotation, x = "Mantel Spearman rho", y = NULL, tag = tag
    ) +
    theme_environment() +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 2.7)), colour = "none") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.25),
      axis.text.y = ggplot2::element_text(size = 7.3),
      aspect.ratio = aspect_ratio
    )
}

build_figures <- function(validated) {
  d <- validated$ds4
  colors <- validated$colors
  pca <- build_pca(d, colors)
  rda16 <- build_rda(d, colors, "16S", "b")
  rda18 <- build_rda(d, colors, "18S", "a")
  mantel16 <- build_mantel(d, colors, "16S", "c", aspect_ratio = NULL)
  mantel18 <- build_mantel(d, colors, "18S", "b", aspect_ratio = 0.70)
  pca <- pca + ggplot2::labs(tag = "a")
  legend_grob <- cowplot::get_legend(
    pca + ggplot2::theme(legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center")
  )
  legend_patch <- patchwork::wrap_elements(full = legend_grob)
  hide_legend <- function(p) p + ggplot2::theme(legend.position = "none")
  figure5 <- patchwork::wrap_plots(
    A = hide_legend(pca), B = hide_legend(rda16),
    C = hide_legend(mantel16), D = legend_patch,
    design = "AB\nCC\nDD"
  ) + patchwork::plot_layout(heights = c(1.0, 0.72, 0.075))
  figure6 <- patchwork::wrap_plots(
    A = hide_legend(rda18), B = hide_legend(mantel18), C = legend_patch,
    design = "AB\nCC"
  ) + patchwork::plot_layout(widths = c(1.02, 0.98), heights = c(1.0, 0.105))
  list(figure5 = figure5, figure6 = figure6)
}

save_pdf <- function(plot, path, width_mm, height_mm) {
  if (file.exists(path)) stopf("Refusing to overwrite: %s", path)
  grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4, family = "sans", onefile = TRUE)
  print(plot)
  grDevices::dev.off()
}
save_tiff <- function(plot, path, width_mm, height_mm) {
  if (file.exists(path)) stopf("Refusing to overwrite: %s", path)
  grDevices::tiff(
    path, width = width_mm, height = height_mm, units = "mm", res = TIFF_DPI,
    compression = "lzw", type = "cairo", family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}
save_preview <- function(plot, path, width_mm, height_mm) {
  if (file.exists(path)) stopf("Refusing to overwrite: %s", path)
  grDevices::png(
    path, width = width_mm, height = height_mm, units = "mm", res = PREVIEW_DPI,
    type = "cairo", family = "sans", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
}

signature_ok <- function(path, kind) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  raw <- readBin(con, what = "raw", n = 8L)
  if (kind == "pdf") return(rawToChar(raw[1:5]) == "%PDF-")
  if (kind == "tiff") {
    sig <- as.integer(raw[1:4])
    return(identical(sig, c(73L, 73L, 42L, 0L)) || identical(sig, c(77L, 77L, 0L, 42L)))
  }
  FALSE
}

make_input_manifest <- function(paths) {
  roles <- c(
    palette_contract = paths$palette,
    frozen_figure_contract = paths$contract_source,
    accepted_ds4_pdf_visual_reference = paths$current_ds4_pdf,
    accepted_ds4_png_visual_reference = paths$current_ds4_png,
    builder_script = paths$script,
    setNames(paths$ds4, paste0("accepted_ds4_", names(paths$ds4))),
    setNames(paths$canonical, paste0("canonical_environment_", names(paths$canonical)))
  )
  data.frame(
    InputRole = names(roles),
    Path = unname(roles),
    Bytes = vapply(unname(roles), file_bytes, numeric(1L)),
    SHA256 = vapply(unname(roles), sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
}

caption_lines <- function() c(
  "# Technical caption facts and proposed legends",
  "",
  "## Figure 5. Environmental associations with the 16S community",
  "",
  "Environmental associations were evaluated exploratorily using Spring Sea-sand samples only. Three spatial field replicates were converted to sample-wise relative abundances and arithmetically averaged to one microbial centroid per Site (n = 7 Sites); each Site contributed one measured environmental record. (a) Standardized PCA of 10 environmental variables (pH, NaCl, T-N, T-S, Ca, Mg, K, Na, T-P, and water content). TOC was excluded because the Sehwa value was recorded as ND/legacy zero without a documented LOD/LOQ. (b) Genus-level 16S Hellinger RDA constrained by EnvPC1 and EnvPC2 (F = 3.305, permutation P = 0.00694, R2 = 0.623, adjusted R2 = 0.435); overall and sequential EnvPC term tests used all 5,039 nonidentity Site-label permutations. (c) Mantel association between ASV site-mean Bray-Curtis dissimilarity and scaled 10-variable environmental Euclidean distance (Spearman rho = 0.644, P = BH q = 0.00615). The full n = 7 Mantel test used all 5,039 nonidentity Site-label permutations; each n = 6 leave-one-Site-out estimate used all 719 nonidentity permutations of the remaining Sites while retaining the full-data environmental standardization. The two primary marker-specific full Mantel P values were BH-adjusted together. The black diamond is the full seven-Site estimate; colored circles are the seven six-Site leave-one-Site-out estimates (rho = 0.446-0.811), not a confidence interval, and their colors denote the omitted Site's descriptive geological context. Environmental values were not imputed or copied to Water or Autumn samples. The analysis is exploratory and does not establish causal environmental drivers or independently replicated geological-context effects.",
  "",
  "## Figure 6. Environmental associations with the 18S community",
  "",
  "Environmental associations were evaluated exploratorily using the same seven Spring Sea-sand Site centroids and environmental records as in Figure 5. The shared environmental PCA used 10 variables and excluded TOC because the Sehwa value was recorded as ND/legacy zero without a documented LOD/LOQ. (a) PR2 class-level 18S Hellinger RDA constrained by EnvPC1 and EnvPC2, as defined by the standardized environmental PCA in Figure 5a (F = 1.643, permutation P = 0.0528, R2 = 0.451, adjusted R2 = 0.176); overall and sequential EnvPC term tests used all 5,039 nonidentity Site-label permutations. (b) Mantel association between ASV site-mean Bray-Curtis dissimilarity and scaled 10-variable environmental Euclidean distance (Spearman rho = 0.666, P = 0.000992, BH q = 0.00198). The full n = 7 Mantel test used all 5,039 nonidentity Site-label permutations; each n = 6 leave-one-Site-out estimate used all 719 nonidentity permutations of the remaining Sites while retaining the full-data environmental standardization. The two primary marker-specific full Mantel P values were BH-adjusted together. The black diamond is the full seven-Site estimate; colored circles are the seven six-Site leave-one-Site-out estimates (rho = 0.521-0.732), not a confidence interval, and their colors denote the omitted Site's descriptive geological context. Environmental values were not imputed or copied to Water or Autumn samples. The analysis is exploratory and does not establish causal environmental drivers or independently replicated geological-context effects.",
  "",
  "## Shared technical notes",
  "",
  "- Full n = 7 Mantel tests used all 5,039 nonidentity Site-label permutations. Each n = 6 LOSO estimate used all 719 nonidentity permutations of the remaining six Sites while retaining the full-data environmental standardization.",
  "- The two primary marker-specific full Mantel P values were BH-adjusted together.",
  "- RDA inference includes the overall model test and sequential term-level tests for EnvPC1 and EnvPC2, each using all 5,039 nonidentity Site-label permutations.",
  "- Environmental arrows are black and Site points use the frozen geological-context palette.",
  "- The exact correspondence of the Dokdo microbial sampling position to the environmental station labeled Dokdo-1 was not documented; no exact co-location claim is made.",
  "- Water and Autumn samples have no environmental-association values in these figures."
)

readme_lines <- function() c(
  "# Revision environmental association figures v1",
  "",
  sprintf("Builder version: %s", SCRIPT_VERSION),
  "",
  "Accepted manuscript-delivery replacements for Reviewer 1, Comment 3.",
  "",
  "- Figure 5: shared 10-variable environmental PCA, 16S genus-level Hellinger RDA, and 16S Mantel/LOSO influence.",
  "- Figure 6: 18S PR2 class-level Hellinger RDA and 18S Mantel/LOSO influence; EnvPC1/2 refer to Figure 5a.",
  "- Scope: exactly seven Spring Sea-sand Site centroids derived from 21 samples.",
  "- Numerical analyses were not rerun. Plots use accepted DS4 source values only.",
  "- Water and Autumn environmental values were neither imputed nor copied.",
  "- TOC is excluded from the primary PCA because Sehwa was ND/legacy zero without documented LOD/LOQ.",
  "- Full n=7 Mantel tests use 5,039 nonidentity permutations; each n=6 LOSO estimate uses 719 permutations while retaining full-data environmental standardization; the two full Mantel P values form one BH family.",
  "- RDA reporting includes overall and sequential EnvPC1/EnvPC2 term-level tests.",
  "- Geological context is descriptive; causal environmental-driver and habitat-replicated claims are prohibited.",
  "- PDF and TIFF delivery files are under revision/rev_ms/figures; PNG files in qc/previews are QA-only.",
  "- The output manifest intentionally excludes itself to avoid recursive hashing."
)

make_statistics_table <- function(validated) {
  o <- validated$overall[match(c("16S", "18S"), validated$overall$Marker), , drop = FALSE]
  m <- validated$mantel[match(c("16S", "18S"), validated$mantel$Marker), , drop = FALSE]
  r <- validated$loso_ranges[match(c("16S", "18S"), validated$loso_ranges$Marker), , drop = FALSE]
  data.frame(
    Marker = o$Marker,
    Population = "Spring Sea-sand Site centroids",
    NSites = 7L,
    RDAFeatureLevel = o$FeatureLevel,
    RDAFormula = "Hellinger community ~ EnvPC1 + EnvPC2",
    RDAF = o$F,
    RDAP = o$P,
    RDAPermutations = o$Permutations,
    RDAR2 = o$R2,
    RDAAdjustedR2 = o$AdjustedR2,
    RDATestContract = "overall_model_plus_sequential_EnvPC1_EnvPC2_tests",
    MantelMethod = m$Method,
    MantelStatistic = m$Statistic,
    MantelP = m$P,
    MantelBHq = m$PAdjustedBHPrimaryAcrossMarkers,
    MantelPermutations = m$Permutations,
    LOSONSites = 6L,
    LOSOPermutationsEach = 719L,
    LOSOEnvironmentalStandardization = "retained_from_full_seven_site_data",
    LOSOMinStatistic = r$LOSOMin,
    LOSOMaxStatistic = r$LOSOMax,
    Interpretation = c("exploratory_association", "exploratory_association_rda_not_conventionally_significant"),
    stringsAsFactors = FALSE
  )
}

make_panel_registry <- function() data.frame(
  Figure = c("Figure 5", "Figure 5", "Figure 5", "Figure 6", "Figure 6"),
  Panel = c("a", "b", "c", "a", "b"),
  Marker = c("shared", "16S", "16S", "18S", "18S"),
  Content = c("standardized environmental PCA and loadings", "genus-level Hellinger RDA", "full Mantel and site-labeled LOSO influence", "PR2 class-level Hellinger RDA", "full Mantel and site-labeled LOSO influence"),
  Population = "Spring Sea-sand; 7 Site centroids from 21 samples",
  NumericalAnalysisRerun = FALSE,
  GeologicalContextRole = "descriptive_only",
  WaterOrAutumnEnvironmentImputed = FALSE,
  stringsAsFactors = FALSE
)

make_qa <- function(validated) {
  list(
    status = "PASS",
    script_version = SCRIPT_VERSION,
    package_version = PACKAGE_VERSION,
    numerical_analysis_rerun = FALSE,
    scope = list(samples = 21L, sites = 7L, season = "Spring", sample_type = "Sea-sand"),
    pca = list(variable_count = 10L, toc_included = FALSE, pc1_percent = 44.8589651245106, pc2_percent = 36.4749453361212),
    rda_rows = list(`16S` = 7L, `18S` = 7L),
    mantel_rows = list(full = 2L, loso_16S = 7L, loso_18S = 7L),
    mantel_permutation_contract = list(full_n7 = 5039L, loso_n6_each = 719L, loso_standardization = "retained_full_data", bh_family_size = 2L),
    rda_permutation_contract = list(test_types = c("overall", "sequential_term_level"), permutations = 5039L),
    water_or_autumn_rows = 0L,
    geological_context_inference = FALSE,
    causal_driver_claim = FALSE,
    dokdo_station_note = "Exact correspondence of the microbial sampling position to environmental station Dokdo-1 was not documented.",
    figure5 = list(width_mm = FIG5_WIDTH_MM, height_mm = FIG5_HEIGHT_MM, panels = 3L),
    figure6 = list(width_mm = FIG6_WIDTH_MM, height_mm = FIG6_HEIGHT_MM, panels = 2L),
    tiff_dpi = TIFF_DPI
  )
}

run_build <- function(paths, validated) {
  assert_true(!dir.exists(paths$package_root), sprintf("Refusing to overwrite package: %s", paths$package_root))
  for (p in paths$deliverables) assert_true(!file.exists(p), sprintf("Refusing to overwrite deliverable: %s", p))
  dir.create(paths$delivery_dir, recursive = TRUE, showWarnings = FALSE)
  stage_package <- paste0(paths$package_root, ".staging-", Sys.getpid())
  stage_delivery <- path_join(paths$delivery_dir, paste0(".environment_figures_v1.staging-", Sys.getpid()))
  assert_true(!file.exists(stage_package) && !file.exists(stage_delivery), "A staging path already exists")
  dir.create(stage_package, recursive = TRUE, showWarnings = FALSE)
  dir.create(stage_delivery, recursive = TRUE, showWarnings = FALSE)
  promoted <- character()
  success <- FALSE
  on.exit({
    if (dir.exists(stage_package)) unlink(stage_package, recursive = TRUE, force = TRUE)
    if (dir.exists(stage_delivery)) unlink(stage_delivery, recursive = TRUE, force = TRUE)
    if (!success && length(promoted)) unlink(promoted[file.exists(promoted)], force = TRUE)
  }, add = TRUE)

  figures <- build_figures(validated)
  staged_delivery <- file.path(stage_delivery, basename(paths$deliverables))
  names(staged_delivery) <- names(paths$deliverables)
  save_pdf(figures$figure5, staged_delivery[["figure5_pdf"]], FIG5_WIDTH_MM, FIG5_HEIGHT_MM)
  save_tiff(figures$figure5, staged_delivery[["figure5_tiff"]], FIG5_WIDTH_MM, FIG5_HEIGHT_MM)
  save_pdf(figures$figure6, staged_delivery[["figure6_pdf"]], FIG6_WIDTH_MM, FIG6_HEIGHT_MM)
  save_tiff(figures$figure6, staged_delivery[["figure6_tiff"]], FIG6_WIDTH_MM, FIG6_HEIGHT_MM)

  preview_dir <- file.path(stage_package, "qc", "previews")
  dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)
  save_preview(figures$figure5, file.path(preview_dir, "Figure_5_16S_environment_associations_v1_preview.png"), FIG5_WIDTH_MM, FIG5_HEIGHT_MM)
  save_preview(figures$figure6, file.path(preview_dir, "Figure_6_18S_environment_associations_v1_preview.png"), FIG6_WIDTH_MM, FIG6_HEIGHT_MM)

  for (nm in names(paths$ds4)) {
    copy_file_new(paths$ds4[[nm]], file.path(stage_package, "source_data", basename(paths$ds4[[nm]])))
  }
  copy_file_new(paths$contract_source, file.path(stage_package, "provenance", "environment_figure_contract_v1.md"))
  write_tsv_new(make_statistics_table(validated), file.path(stage_package, "tables", "environment_figure_statistics_v1.tsv"))
  write_tsv_new(make_panel_registry(), file.path(stage_package, "tables", "environment_figure_panel_registry_v1.tsv"))
  write_text_new(caption_lines(), file.path(stage_package, "captions", "environment_figure_caption_facts_v1.md"))
  write_text_new(readme_lines(), file.path(stage_package, "provenance", "README_environment_figures_v1.md"))
  write_tsv_new(make_input_manifest(paths), file.path(stage_package, "provenance", "input_sha256_manifest_v1.tsv"))
  jsonlite::write_json(make_qa(validated), file.path(stage_package, "qc", "environment_figure_validation_v1.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)

  for (p in staged_delivery[c("figure5_pdf", "figure6_pdf")]) {
    assert_true(signature_ok(p, "pdf"), sprintf("Invalid PDF signature: %s", p))
    info <- pdftools::pdf_info(p)
    assert_true(info$pages == 1L, sprintf("PDF must have one page: %s", p))
  }
  for (p in staged_delivery[c("figure5_tiff", "figure6_tiff")]) assert_true(signature_ok(p, "tiff"), sprintf("Invalid TIFF signature: %s", p))

  package_files <- list.files(stage_package, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  package_files <- package_files[file.info(package_files)$isdir %in% FALSE]
  output_manifest <- rbind(
    data.frame(
      Scope = "package",
      RelativePath = substring(gsub("\\\\", "/", package_files), nchar(gsub("\\\\", "/", stage_package)) + 2L),
      Bytes = vapply(package_files, file_bytes, numeric(1L)),
      SHA256 = vapply(package_files, sha256_file, character(1L)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      Scope = "delivery",
      RelativePath = substring(gsub("\\\\", "/", paths$deliverables), nchar(gsub("\\\\", "/", paths$root)) + 2L),
      Bytes = vapply(staged_delivery, file_bytes, numeric(1L)),
      SHA256 = vapply(staged_delivery, sha256_file, character(1L)),
      stringsAsFactors = FALSE
    )
  )
  output_manifest <- output_manifest[order(output_manifest$Scope, output_manifest$RelativePath), , drop = FALSE]
  write_tsv_new(output_manifest, file.path(stage_package, "provenance", "output_sha256_manifest_v1.tsv"))

  for (nm in names(paths$deliverables)) {
    assert_true(file.rename(staged_delivery[[nm]], paths$deliverables[[nm]]), sprintf("Could not promote deliverable: %s", paths$deliverables[[nm]]))
    promoted <- c(promoted, paths$deliverables[[nm]])
  }
  unlink(stage_delivery, recursive = TRUE, force = TRUE)
  assert_true(file.rename(stage_package, paths$package_root), "Could not promote the versioned package")
  success <- TRUE
  invisible(TRUE)
}

check_build <- function(paths) {
  assert_true(dir.exists(paths$package_root), sprintf("Missing package: %s", paths$package_root))
  input_manifest_path <- file.path(paths$package_root, "provenance", "input_sha256_manifest_v1.tsv")
  input_manifest <- read_tsv(input_manifest_path)
  require_columns(input_manifest, c("InputRole", "Path", "Bytes", "SHA256"), "input manifest")
  assert_true(all(file.exists(input_manifest$Path)), "Input manifest contains missing source files")
  current_input_bytes <- vapply(input_manifest$Path, file_bytes, numeric(1L))
  current_input_hash <- vapply(input_manifest$Path, sha256_file, character(1L))
  assert_true(identical(as.numeric(input_manifest$Bytes), unname(current_input_bytes)), "Input manifest size mismatch")
  assert_true(identical(tolower(input_manifest$SHA256), unname(current_input_hash)), "Input manifest SHA256 mismatch")
  manifest_path <- file.path(paths$package_root, "provenance", "output_sha256_manifest_v1.tsv")
  manifest <- read_tsv(manifest_path)
  require_columns(manifest, c("Scope", "RelativePath", "Bytes", "SHA256"), "output manifest")
  resolved <- ifelse(
    manifest$Scope == "package",
    file.path(paths$package_root, manifest$RelativePath),
    file.path(paths$root, manifest$RelativePath)
  )
  assert_true(all(file.exists(resolved)), "Output manifest contains missing files")
  observed_bytes <- vapply(resolved, file_bytes, numeric(1L))
  observed_hash <- vapply(resolved, sha256_file, character(1L))
  assert_true(identical(as.numeric(manifest$Bytes), unname(observed_bytes)), "Output manifest size mismatch")
  assert_true(identical(tolower(manifest$SHA256), unname(observed_hash)), "Output manifest SHA256 mismatch")
  expected_package <- sort(c(
    manifest$RelativePath[manifest$Scope == "package"],
    "provenance/output_sha256_manifest_v1.tsv"
  ))
  actual_package <- list.files(paths$package_root, recursive = TRUE, full.names = FALSE, all.files = FALSE)
  actual_package <- sort(gsub("\\\\", "/", actual_package[file.info(file.path(paths$package_root, actual_package))$isdir %in% FALSE]))
  assert_true(identical(expected_package, actual_package), "Package contains unexpected or missing files")
  for (p in paths$deliverables[c("figure5_pdf", "figure6_pdf")]) {
    assert_true(signature_ok(p, "pdf"), sprintf("Invalid PDF signature: %s", p))
    assert_true(pdftools::pdf_info(p)$pages == 1L, sprintf("PDF must contain one page: %s", p))
  }
  for (p in paths$deliverables[c("figure5_tiff", "figure6_tiff")]) assert_true(signature_ok(p, "tiff"), sprintf("Invalid TIFF signature: %s", p))
  qa <- jsonlite::read_json(file.path(paths$package_root, "qc", "environment_figure_validation_v1.json"), simplifyVector = TRUE)
  assert_true(identical(qa$status, "PASS"), "QA status is not PASS")
  cat(sprintf("CHECK_OK package=%s manifest_rows=%d\n", paths$package_root, nrow(manifest)))
  invisible(TRUE)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  paths <- resolve_paths(args$project_root)
  validated <- load_and_validate(paths)
  if (args$mode == "dry-run") {
    assert_true(!dir.exists(paths$package_root), sprintf("Target package already exists: %s", paths$package_root))
    for (p in paths$deliverables) assert_true(!file.exists(p), sprintf("Target deliverable already exists: %s", p))
    cat(sprintf(
      "DRY_RUN_OK sites=7 samples=21 pca_variables=10 rda_rows=%d mantel_loso_rows=%d output_package=%s\n",
      nrow(validated$ds4$rda_sites), nrow(validated$ds4$mantel_loso), paths$package_root
    ))
  } else if (args$mode == "run") {
    run_build(paths, validated)
    tryCatch(
      check_build(paths),
      error = function(e) {
        for (p in paths$deliverables) if (file.exists(p)) unlink(p, force = TRUE)
        if (dir.exists(paths$package_root)) unlink(paths$package_root, recursive = TRUE, force = TRUE)
        stopf("Post-build check failed; newly created outputs were rolled back: %s", conditionMessage(e))
      }
    )
    cat("RUN_OK\n")
  } else {
    check_build(paths)
  }
}

main()
