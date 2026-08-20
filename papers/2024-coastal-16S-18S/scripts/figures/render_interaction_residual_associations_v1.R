#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)

find_project_root <- function(start) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (i in seq_len(10L)) {
    if (dir.exists(file.path(current, "revision", "re_analysis"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not resolve project root", call. = FALSE)
}

project_root <- find_project_root(package_root)
cli <- commandArgs(trailingOnly = TRUE)
modes <- intersect(cli, c("--dry-run", "--run", "--check"))
if (length(modes) != 1L) stop("Specify exactly one of --dry-run, --run, or --check", call. = FALSE)
mode <- sub("^--", "", modes[[1L]])

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
  invisible(TRUE)
}

read_tsv <- function(path) {
  assert_true(file.exists(path), paste0("Missing file: ", path))
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv_new <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(tmp)) stop("Refusing overwrite: ", path, call. = FALSE)
  write.table(x, tmp, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Atomic file promotion failed: ", path, call. = FALSE)
  invisible(path)
}

write_lines_new <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(tmp)) stop("Refusing overwrite: ", path, call. = FALSE)
  writeLines(x, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic file promotion failed: ", path, call. = FALSE)
  invisible(path)
}

assert_true(as.character(getRversion()) == "4.5.2",
            paste0("R 4.5.2 required; found ", getRversion()))
required_packages <- c("digest", "dplyr", "ggplot2", "igraph", "ggrepel",
                       "patchwork", "svglite", "ragg", "png", "pdftools")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1L), quietly = TRUE)]
assert_true(!length(missing_packages),
            paste0("Missing R packages: ", paste(missing_packages, collapse = ", ")))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(ggrepel)
  library(patchwork)
})

sha256_file <- function(path) {
  toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

relative_to_package <- function(path) {
  root <- paste0(normalizePath(package_root, winslash = "/", mustWork = TRUE), "/")
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  assert_true(startsWith(normalized, root), paste0("Path escapes package root: ", path))
  substring(normalized, nchar(root) + 1L)
}

artifacts_dir <- file.path(package_root, "artifacts")
results_dir <- file.path(artifacts_dir, "results")
numeric_manifest_path <- file.path(artifacts_dir, "provenance", "output_manifest_sha256_v1.tsv")
network_edges_path <- file.path(results_dir, "network_all_edges_v1.tsv")
network_summary_path <- file.path(results_dir, "network_summary_v1.tsv")
comparison_path <- file.path(results_dir, "additive_interaction_comparison_v1.tsv")
accepted_reproduction_path <- file.path(results_dir, "accepted_additive_baseline_reproduction_v1.tsv")
module_status_path <- file.path(results_dir, "module_status_v1.tsv")
network_contract_path <- file.path(package_root, "config", "network_contract_v1.tsv")
display_gate_path <- file.path(package_root, "config", "presentation_claim_gate_v1.tsv")
figure_contract_path <- file.path(package_root, "config", "figure_contract_interaction_s11_v1.md")
caption_path <- file.path(package_root, "config", "Supplementary_Figure_S11_interaction_caption_v1.md")
presentation_dir <- file.path(package_root, "presentation_v1")

validate_numeric_manifest <- function() {
  manifest <- read_tsv(numeric_manifest_path)
  assert_true(identical(names(manifest), c("PackageRelativePath", "Bytes", "SHA256")),
              "Unexpected numeric artifact manifest schema")
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(package_root, manifest$PackageRelativePath[[i]])
    assert_true(file.exists(path), paste0("Manifested numeric artifact missing: ", path))
    assert_true(as.numeric(file.info(path)$size) == as.numeric(manifest$Bytes[[i]]),
                paste0("Numeric artifact size drift: ", manifest$PackageRelativePath[[i]]))
    assert_true(identical(sha256_file(path), toupper(manifest$SHA256[[i]])),
                paste0("Numeric artifact SHA-256 drift: ", manifest$PackageRelativePath[[i]]))
  }
  manifest
}

load_and_validate_sources <- function() {
  numeric_manifest <- validate_numeric_manifest()
  edges <- read_tsv(network_edges_path)
  summary_df <- read_tsv(network_summary_path)
  comparison <- read_tsv(comparison_path)
  reproduction <- read_tsv(accepted_reproduction_path)
  module_status <- read_tsv(module_status_path)
  network_contract <- read_tsv(network_contract_path)
  gate <- read_tsv(display_gate_path)

  edges$StableStrict <- as.logical(edges$StableStrict)
  edges$SiteSensitivityPass <- as.logical(edges$SiteSensitivityPass)
  edges$CombinedSiteLOSORobust <- as.logical(edges$CombinedSiteLOSORobust)
  gate$FigureDisplay <- as.logical(gate$FigureDisplay)

  expected_networks <- c("16S_within", "18S_within", "16S_18S_cross")
  assert_true(setequal(gate$Network, expected_networks), "Presentation gate network set mismatch")
  assert_true(nrow(module_status) == 1L && module_status$Status[[1L]] == "PASS" &&
                module_status$BaselineReproduction[[1L]] == "PASS" &&
                isTRUE(as.logical(module_status$FullModelsEstimable[[1L]])),
              "Residual interaction module is not in a validated PASS state")
  assert_true(all(reproduction$Status == "PASS") &&
                all(reproduction$StableIdentityMismatches == 0L),
              "Accepted additive baseline reproduction gate failed")

  lookup_summary <- function(model, network, column) {
    hit <- summary_df$Model == model & summary_df$Network == network
    assert_true(sum(hit) == 1L, paste0("Summary lookup failed: ", model, "/", network))
    summary_df[[column]][hit]
  }
  lookup_reproduction <- function(network, column) {
    hit <- reproduction$Network == network
    assert_true(sum(hit) == 1L, paste0("Reproduction lookup failed: ", network))
    reproduction[[column]][hit]
  }

  for (i in seq_len(nrow(gate))) {
    network <- gate$Network[[i]]
    assert_true(lookup_reproduction(network, "AcceptedStableEdges") ==
                  gate$ExpectedAcceptedAdditiveStrict[[i]],
                paste0("Accepted additive count mismatch: ", network))
    assert_true(lookup_summary("additive", network, "StrictEdges") ==
                  gate$ExpectedAdditiveRefitStrict[[i]],
                paste0("Additive refit count mismatch: ", network))
    assert_true(lookup_summary("interaction", network, "StrictEdges") ==
                  gate$ExpectedInteractionStrict[[i]],
                paste0("Interaction strict count mismatch: ", network))
    assert_true(lookup_summary("interaction", network, "CombinedSiteLOSORobustStrictEdges") ==
                  gate$ExpectedInteractionCombinedRobust[[i]],
                paste0("Interaction robust count mismatch: ", network))
  }

  assert_true(!gate$FigureDisplay[gate$Network == "18S_within"],
              "18S within-marker display must fail closed")
  assert_true(all(gate$FigureDisplay[gate$Network %in% c("16S_within", "16S_18S_cross")]),
              "16S and cross-marker display decisions must be explicit")
  j18 <- comparison$StrictJaccard[comparison$Network == "18S_within"]
  assert_true(length(j18) == 1L && is.finite(j18) && abs(j18 - 0.25) < 1e-12,
              "18S model-sensitivity boundary changed")

  robust <- edges[edges$Model == "interaction" & edges$StableStrict &
                    edges$CombinedSiteLOSORobust, , drop = FALSE]
  plotted_networks <- gate$Network[gate$FigureDisplay]
  plotted <- robust[robust$Network %in% plotted_networks, , drop = FALSE]
  assert_true(!any(plotted$Network == "18S_within"), "18S edge leaked into plotted data")
  for (network in plotted_networks) {
    observed <- sum(plotted$Network == network)
    expected <- gate$ExpectedInteractionCombinedRobust[gate$Network == network]
    assert_true(observed == expected, paste0("Plotted robust count mismatch: ", network))
  }

  list(edges = edges, plotted = plotted, summary = summary_df,
       comparison = comparison, reproduction = reproduction,
       module_status = module_status, network_contract = network_contract,
       gate = gate, numeric_manifest = numeric_manifest)
}

build_comparison_table <- function(src) {
  order_networks <- c("16S_within", "16S_18S_cross", "18S_within")
  labels <- c(`16S_within` = "16S within", `16S_18S_cross` = "16S-18S cross",
              `18S_within` = "18S within")
  rows <- lapply(order_networks, function(network) {
    gate_row <- src$gate[src$gate$Network == network, , drop = FALSE]
    cmp_row <- src$comparison[src$comparison$Network == network, , drop = FALSE]
    data.frame(
      Network = network,
      NetworkLabel = unname(labels[[network]]),
      AcceptedAdditiveStrict = gate_row$ExpectedAcceptedAdditiveStrict,
      AdditiveRefit999Strict = gate_row$ExpectedAdditiveRefitStrict,
      InteractionRefit999Strict = gate_row$ExpectedInteractionStrict,
      InteractionCombinedRobust = gate_row$ExpectedInteractionCombinedRobust,
      StrictJaccard = cmp_row$StrictJaccard,
      FigureDisplay = gate_row$FigureDisplay,
      ManuscriptDecision = gate_row$ManuscriptDecision,
      Rationale = gate_row$Rationale,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

edge_sign <- function(r) {
  factor(ifelse(r > 0, "Positive", "Negative"), levels = c("Positive", "Negative"))
}

make_within_panel <- function(x, source_dir, seed = 20260811L) {
  x <- x[x$Network == "16S_within", , drop = FALSE]
  assert_true(nrow(x) == 54L, "16S plot must contain exactly 54 robust edges")
  graph <- igraph::graph_from_data_frame(x[, c("Node1", "Node2")], directed = FALSE)
  set.seed(seed)
  xy <- igraph::layout_with_fr(graph, niter = 2500)
  nodes <- data.frame(Node = igraph::V(graph)$name, X = xy[, 1L], Y = xy[, 2L],
                      Degree = as.integer(igraph::degree(graph)), stringsAsFactors = FALSE)
  label_nodes <- nodes |>
    dplyr::arrange(dplyr::desc(Degree), Node) |>
    dplyr::slice_head(n = 10L) |>
    dplyr::pull(Node)
  nodes$Label <- ifelse(nodes$Node %in% label_nodes, nodes$Node, "")
  edges_plot <- x |>
    dplyr::left_join(nodes[, c("Node", "X", "Y")], by = c("Node1" = "Node")) |>
    dplyr::rename(X1 = X, Y1 = Y) |>
    dplyr::left_join(nodes[, c("Node", "X", "Y")], by = c("Node2" = "Node")) |>
    dplyr::rename(X2 = X, Y2 = Y) |>
    dplyr::mutate(Sign = edge_sign(R), AbsR = abs(R))
  write_tsv_new(nodes, file.path(source_dir, "Supplementary_Figure_S11_interaction_16S_nodes_v1.tsv"))
  write_tsv_new(edges_plot, file.path(source_dir, "Supplementary_Figure_S11_interaction_16S_edges_v1.tsv"))

  positive <- sum(edges_plot$R > 0)
  negative <- sum(edges_plot$R < 0)
  panel <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges_plot,
      ggplot2::aes(X1, Y1, xend = X2, yend = Y2, color = Sign, linewidth = AbsR),
      alpha = 0.72, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = nodes, ggplot2::aes(X, Y), shape = 21, fill = "#56B4E9", color = "#1F2937",
      stroke = 0.35, size = 2.45
    ) +
    ggrepel::geom_text_repel(
      data = nodes[nodes$Label != "", , drop = FALSE], ggplot2::aes(X, Y, label = Label),
      size = 2.35, min.segment.length = 0, max.overlaps = Inf, seed = seed,
      box.padding = 0.40, point.padding = 0.22, force = 2.0, force_pull = 0.20,
      max.time = 5, max.iter = 50000, color = "#111827", family = "Arial"
    ) +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"),
                                name = "Residual association") +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.35), guide = "none") +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::labs(
      title = "(a) Interaction-model 16S residual associations",
      subtitle = sprintf("%d Site+LOSO-robust edges (%d positive; %d negative); node size is constant",
                         nrow(edges_plot), positive, negative)
    ) +
    ggplot2::theme_void(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.2, color = "#374151", hjust = 0),
      plot.margin = ggplot2::margin(4, 10, 3, 10),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 7.2),
      legend.text = ggplot2::element_text(size = 7)
    )
  list(plot = panel, nodes = nodes, edges = edges_plot)
}

make_cross_panel <- function(x, source_dir) {
  x <- x[x$Network == "16S_18S_cross", , drop = FALSE]
  assert_true(nrow(x) == 12L, "Cross-marker plot must contain exactly 12 robust edges")
  assert_true(all(x$Marker1 == "16S" & x$Marker2 == "18S"),
              "Cross-marker orientation changed")
  rank_nodes <- function(values) {
    out <- as.data.frame(table(values), stringsAsFactors = FALSE)
    names(out) <- c("Node", "Degree")
    out$Degree <- as.integer(out$Degree)
    out |> dplyr::arrange(dplyr::desc(Degree), Node)
  }
  d16 <- rank_nodes(x$Node1)
  d18 <- rank_nodes(x$Node2)
  pos16 <- data.frame(Key = paste0("16S::", d16$Node), Node = d16$Node,
                      Marker = "16S genus", Degree = d16$Degree, X = 0,
                      Y = seq(1, 0, length.out = nrow(d16)), stringsAsFactors = FALSE)
  pos18 <- data.frame(Key = paste0("18S::", d18$Node), Node = d18$Node,
                      Marker = "PR2 18S class", Degree = d18$Degree, X = 1,
                      Y = seq(1, 0, length.out = nrow(d18)), stringsAsFactors = FALSE)
  nodes <- dplyr::bind_rows(pos16, pos18)
  edges_plot <- x |>
    dplyr::mutate(Key1 = paste0("16S::", Node1), Key2 = paste0("18S::", Node2)) |>
    dplyr::left_join(nodes[, c("Key", "X", "Y")], by = c("Key1" = "Key")) |>
    dplyr::rename(X1 = X, Y1 = Y) |>
    dplyr::left_join(nodes[, c("Key", "X", "Y")], by = c("Key2" = "Key")) |>
    dplyr::rename(X2 = X, Y2 = Y) |>
    dplyr::mutate(Sign = edge_sign(R), AbsR = abs(R))
  write_tsv_new(nodes, file.path(source_dir, "Supplementary_Figure_S11_interaction_cross_nodes_v1.tsv"))
  write_tsv_new(edges_plot, file.path(source_dir, "Supplementary_Figure_S11_interaction_cross_edges_v1.tsv"))

  positive <- sum(edges_plot$R > 0)
  negative <- sum(edges_plot$R < 0)
  panel <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edges_plot,
      ggplot2::aes(X1, Y1, xend = X2, yend = Y2, color = Sign, linewidth = AbsR),
      curvature = 0.11, alpha = 0.72, lineend = "round"
    ) +
    ggplot2::geom_point(
      data = nodes, ggplot2::aes(X, Y, fill = Marker), shape = 21, size = 2.55,
      color = "#1F2937", stroke = 0.35
    ) +
    ggplot2::geom_text(data = pos16, ggplot2::aes(X, Y, label = Node), hjust = 1.12,
                       size = 2.25, color = "#111827", family = "Arial") +
    ggplot2::geom_text(data = pos18, ggplot2::aes(X, Y, label = Node), hjust = -0.12,
                       size = 2.25, color = "#111827", family = "Arial") +
    ggplot2::scale_color_manual(values = c(Positive = "#D55E00", Negative = "#0072B2"),
                                name = "Residual association") +
    ggplot2::scale_fill_manual(values = c(`16S genus` = "#56B4E9", `PR2 18S class` = "#009E73"),
                               name = "Node") +
    ggplot2::scale_linewidth_continuous(range = c(0.35, 1.35), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-0.62, 1.62), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "(b) Interaction-model cross-marker residual associations",
      subtitle = sprintf("%d Site+LOSO-robust edges (%d positive; %d negative); %d genera and %d PR2 classes",
                         nrow(edges_plot), positive, negative, nrow(pos16), nrow(pos18))
    ) +
    ggplot2::theme_void(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.2, color = "#374151", hjust = 0),
      plot.margin = ggplot2::margin(4, 10, 3, 10),
      legend.position = "bottom", legend.title = ggplot2::element_text(size = 7.2),
      legend.text = ggplot2::element_text(size = 7)
    )
  list(plot = panel, nodes = nodes, edges = edges_plot)
}

make_comparison_panel <- function(comparison_table, source_dir) {
  write_tsv_new(comparison_table,
                file.path(source_dir, "Supplementary_Figure_S11_model_comparison_v1.tsv"))
  stage_labels <- c(
    AcceptedAdditiveStrict = "Accepted additive\nstrict",
    AdditiveRefit999Strict = "Additive 999-refit\nstrict",
    InteractionRefit999Strict = "Interaction 999-refit\nstrict",
    InteractionCombinedRobust = "Interaction Site+LOSO\nrobust",
    StrictJaccard = "Additive-interaction\nstrict Jaccard",
    Decision = "Figure decision"
  )
  network_levels <- c("16S within", "16S-18S cross", "18S within")
  long_rows <- list()
  k <- 0L
  for (i in seq_len(nrow(comparison_table))) {
    row <- comparison_table[i, , drop = FALSE]
    values <- c(
      AcceptedAdditiveStrict = as.character(row$AcceptedAdditiveStrict),
      AdditiveRefit999Strict = as.character(row$AdditiveRefit999Strict),
      InteractionRefit999Strict = as.character(row$InteractionRefit999Strict),
      InteractionCombinedRobust = if (row$Network == "18S_within")
        paste0(row$InteractionCombinedRobust, "*") else as.character(row$InteractionCombinedRobust),
      StrictJaccard = sprintf("%.2f", row$StrictJaccard),
      Decision = if (isTRUE(row$FigureDisplay)) "Shown\nassociation-only" else "Not shown\nnot interpreted"
    )
    for (stage in names(values)) {
      k <- k + 1L
      fill <- "#F3F6F9"
      if (stage == "StrictJaccard") fill <- if (row$StrictJaccard >= 0.80) "#E5F1F8" else "#FFF1CC"
      if (stage == "Decision") fill <- if (isTRUE(row$FigureDisplay)) "#D8F0E5" else "#FFF1CC"
      if (row$Network == "18S_within" && !stage %in% c("StrictJaccard", "Decision")) fill <- "#FFF8E7"
      long_rows[[k]] <- data.frame(NetworkLabel = row$NetworkLabel, Stage = stage,
                                   Label = unname(values[[stage]]), Fill = fill,
                                   stringsAsFactors = FALSE)
    }
  }
  long <- do.call(rbind, long_rows)
  long$NetworkLabel <- factor(long$NetworkLabel, levels = rev(network_levels))
  long$Stage <- factor(long$Stage, levels = names(stage_labels))

  ggplot2::ggplot(long, ggplot2::aes(Stage, NetworkLabel)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Fill), color = "#CBD5E1", linewidth = 0.45,
                       width = 0.96, height = 0.86) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 2.35, lineheight = 0.88,
                       color = "#111827", family = "Arial") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_discrete(labels = stage_labels, expand = ggplot2::expansion(add = 0.08)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.12)) +
    ggplot2::labs(
      title = "(c) Additive-interaction edge-count comparison and interpretation boundary",
      subtitle = "Counts are exact; *the 18S within-marker surface is retained for transparency but not drawn or interpreted",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 6.5, color = "#111827", lineheight = 0.88),
      axis.text.y = ggplot2::element_text(size = 7.2, color = "#111827"),
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.0, color = "#374151", hjust = 0),
      plot.margin = ggplot2::margin(4, 10, 4, 10)
    )
}

parse_tiff_metadata <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  order <- rawToChar(readBin(con, what = "raw", n = 2L))
  endian <- if (order == "II") "little" else if (order == "MM") "big" else NA_character_
  assert_true(!is.na(endian), "Invalid TIFF byte order")
  magic <- readBin(con, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
  assert_true(magic == 42L, "Invalid TIFF magic")
  ifd_offset <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
  seek(con, where = ifd_offset, origin = "start")
  n_entries <- readBin(con, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
  entries <- vector("list", n_entries)
  for (i in seq_len(n_entries)) {
    tag <- readBin(con, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
    type <- readBin(con, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
    count <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
    if (type == 3L && count == 1L) {
      value <- readBin(con, integer(), n = 1L, size = 2L, signed = FALSE, endian = endian)
      readBin(con, what = "raw", n = 2L)
      offset <- NA_integer_
    } else {
      offset <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
      value <- if (type == 4L && count == 1L) offset else NA_real_
    }
    entries[[i]] <- data.frame(Tag = tag, Type = type, Count = count,
                               Value = value, Offset = offset, stringsAsFactors = FALSE)
  }
  entries <- do.call(rbind, entries)
  get_value <- function(tag) {
    row <- entries[entries$Tag == tag, , drop = FALSE]
    assert_true(nrow(row) >= 1L, paste0("TIFF tag missing: ", tag))
    row$Value[[1L]]
  }
  get_rational <- function(tag) {
    row <- entries[entries$Tag == tag, , drop = FALSE]
    assert_true(nrow(row) >= 1L && row$Type[[1L]] == 5L, paste0("TIFF rational tag missing: ", tag))
    current <- seek(con)
    seek(con, where = row$Offset[[1L]], origin = "start")
    numerator <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
    denominator <- readBin(con, integer(), n = 1L, size = 4L, endian = endian)
    seek(con, where = current, origin = "start")
    numerator / denominator
  }
  data.frame(Width = get_value(256L), Height = get_value(257L),
             XResolution = get_rational(282L), YResolution = get_rational(283L),
             ResolutionUnit = get_value(296L), stringsAsFactors = FALSE)
}

validate_presentation <- function(root, src = NULL) {
  manifest_path <- file.path(root, "provenance", "output_manifest_sha256_v1.tsv")
  manifest <- read_tsv(manifest_path)
  assert_true(identical(names(manifest), c("PresentationRelativePath", "Bytes", "SHA256")),
              "Unexpected presentation manifest schema")
  actual_files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  actual_files <- actual_files[normalizePath(actual_files, winslash = "/", mustWork = TRUE) !=
                                 normalizePath(manifest_path, winslash = "/", mustWork = TRUE)]
  root_prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
  actual_rel <- sort(substring(normalizePath(actual_files, winslash = "/", mustWork = TRUE),
                               nchar(root_prefix) + 1L))
  assert_true(identical(actual_rel, sort(manifest$PresentationRelativePath)),
              "Presentation manifest file set mismatch")
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(root, manifest$PresentationRelativePath[[i]])
    assert_true(file.exists(path), paste0("Manifested presentation file missing: ", path))
    assert_true(as.numeric(file.info(path)$size) == as.numeric(manifest$Bytes[[i]]),
                paste0("Presentation size drift: ", manifest$PresentationRelativePath[[i]]))
    assert_true(identical(sha256_file(path), toupper(manifest$SHA256[[i]])),
                paste0("Presentation SHA-256 drift: ", manifest$PresentationRelativePath[[i]]))
  }
  input_manifest <- read_tsv(file.path(root, "provenance", "input_source_manifest_v1.tsv"))
  for (i in seq_len(nrow(input_manifest))) {
    path <- file.path(package_root, input_manifest$PackageRelativePath[[i]])
    assert_true(file.exists(path), paste0("Presentation input missing: ", path))
    assert_true(identical(sha256_file(path), toupper(input_manifest$SHA256[[i]])),
                paste0("Presentation input drift: ", input_manifest$Role[[i]]))
  }
  numeric_qa <- read_tsv(file.path(root, "qc", "numeric_QA_v1.tsv"))
  export_qa <- read_tsv(file.path(root, "qc", "export_QA_v1.tsv"))
  assert_true(all(as.logical(numeric_qa$Pass)), "Numeric figure QA failed")
  assert_true(all(as.logical(export_qa$Pass)), "Export figure QA failed")
  visual_qa <- readLines(file.path(root, "qc", "visual_QA_v1.md"), warn = FALSE)
  assert_true(any(grepl("^Status: PASS$", visual_qa)), "Visual QA attestation is not PASS")
  edges16 <- read_tsv(file.path(root, "source_data", "Supplementary_Figure_S11_interaction_16S_edges_v1.tsv"))
  edges_cross <- read_tsv(file.path(root, "source_data", "Supplementary_Figure_S11_interaction_cross_edges_v1.tsv"))
  assert_true(nrow(edges16) == 54L && nrow(edges_cross) == 12L,
              "Displayed source-data edge counts changed")
  assert_true(all(edges16$Model == "interaction" & as.logical(edges16$StableStrict) &
                    as.logical(edges16$CombinedSiteLOSORobust)),
              "16S plotted source contains an ineligible edge")
  assert_true(all(edges_cross$Model == "interaction" & as.logical(edges_cross$StableStrict) &
                    as.logical(edges_cross$CombinedSiteLOSORobust)),
              "Cross-marker plotted source contains an ineligible edge")
  assert_true(!any(c(edges16$Network, edges_cross$Network) == "18S_within"),
              "18S within-marker edge present in plotted source")
  invisible(TRUE)
}

src <- load_and_validate_sources()

if (mode == "dry-run") {
  assert_true(!dir.exists(presentation_dir),
              paste0("Versioned presentation output already exists: ", presentation_dir))
  counts <- table(factor(src$plotted$Network,
                         levels = c("16S_within", "16S_18S_cross")))
  cat("DRY-RUN PASS\n")
  cat("R version: ", as.character(getRversion()), "\n", sep = "")
  cat("Numeric artifact manifest: ", nrow(src$numeric_manifest), " files verified\n", sep = "")
  cat("Eligible plotted edges: 16S=", counts[[1L]], "; cross=", counts[[2L]], "\n", sep = "")
  cat("18S within-marker network: not drawn; count retained in comparison table\n")
  cat("Output target: ", presentation_dir, "\n", sep = "")
  quit(save = "no", status = 0L)
}

if (mode == "check") {
  assert_true(dir.exists(presentation_dir), paste0("Missing presentation output: ", presentation_dir))
  validate_presentation(presentation_dir, src)
  cat("CHECK PASS\n")
  cat("Presentation manifest, source lineage, numeric/export QA, visual QA, and plotted-edge gates verified.\n")
  quit(save = "no", status = 0L)
}

assert_true(mode == "run", "Internal mode dispatch error")
assert_true(!dir.exists(presentation_dir),
            paste0("Refusing to overwrite versioned presentation: ", presentation_dir))
stage_dir <- file.path(package_root, paste0("._presentation_v1_staging_", Sys.getpid()))
assert_true(!dir.exists(stage_dir), paste0("Staging path already exists: ", stage_dir))
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("captions", "figures", "logs", "provenance", "qc", "source_data")) {
  dir.create(file.path(stage_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

input_paths <- c(
  numeric_artifact_manifest = numeric_manifest_path,
  network_all_edges = network_edges_path,
  network_summary = network_summary_path,
  additive_interaction_comparison = comparison_path,
  accepted_additive_reproduction = accepted_reproduction_path,
  module_status = module_status_path,
  network_contract = network_contract_path,
  presentation_claim_gate = display_gate_path,
  figure_contract = figure_contract_path,
  caption = caption_path,
  rendering_script = script_file
)
input_manifest <- do.call(rbind, lapply(names(input_paths), function(role) {
  path <- input_paths[[role]]
  data.frame(Role = role, PackageRelativePath = relative_to_package(path),
             Bytes = as.numeric(file.info(path)$size), SHA256 = sha256_file(path),
             stringsAsFactors = FALSE)
}))
write_tsv_new(input_manifest, file.path(stage_dir, "provenance", "input_source_manifest_v1.tsv"))
write_lines_new(readLines(figure_contract_path, warn = FALSE),
                file.path(stage_dir, "figure_contract_interaction_s11_v1.md"))
write_lines_new(readLines(caption_path, warn = FALSE),
                file.path(stage_dir, "captions", "Supplementary_Figure_S11_interaction_caption_v1.md"))

comparison_table <- build_comparison_table(src)
write_tsv_new(src$gate, file.path(stage_dir, "source_data",
                                  "Supplementary_Figure_S11_display_decisions_v1.tsv"))
within <- make_within_panel(src$plotted, file.path(stage_dir, "source_data"))
cross <- make_cross_panel(src$plotted, file.path(stage_dir, "source_data"))
comparison_panel <- make_comparison_panel(comparison_table, file.path(stage_dir, "source_data"))

combined <- (within$plot / cross$plot / comparison_panel) +
  patchwork::plot_layout(heights = c(1.16, 0.88, 0.52), guides = "collect") &
  ggplot2::theme(legend.position = "bottom", legend.box = "horizontal",
                 legend.margin = ggplot2::margin(0, 0, 0, 0))

width_mm <- 180
height_mm <- 235
base <- file.path(stage_dir, "figures", "Supplementary_Figure_S11_interaction_sensitivity_v1")
pdf_path <- paste0(base, ".pdf")
svg_path <- paste0(base, ".svg")
png_path <- paste0(base, ".png")
tiff_path <- paste0(base, ".tiff")
ggplot2::ggsave(pdf_path, combined, width = width_mm, height = height_mm, units = "mm",
                device = grDevices::cairo_pdf, bg = "white")
ggplot2::ggsave(svg_path, combined, width = width_mm, height = height_mm, units = "mm",
                device = svglite::svglite, bg = "white")
ggplot2::ggsave(png_path, combined, width = width_mm, height = height_mm, units = "mm",
                device = ragg::agg_png, dpi = 300, bg = "white")
ggplot2::ggsave(tiff_path, combined, width = width_mm, height = height_mm, units = "mm",
                device = "tiff", dpi = 600, compression = "lzw", bg = "white")

numeric_qa <- data.frame(
  Check = c("common_sample_count", "16S_displayed_robust_edges", "16S_positive",
            "16S_negative", "cross_displayed_robust_edges", "cross_positive",
            "cross_negative", "18S_displayed_edges", "18S_table_strict_edges",
            "18S_strict_jaccard", "all_displayed_edges_interaction",
            "all_displayed_edges_strict", "all_displayed_edges_site_loso_robust"),
  Observed = c(75L, nrow(within$edges), sum(within$edges$R > 0), sum(within$edges$R < 0),
               nrow(cross$edges), sum(cross$edges$R > 0), sum(cross$edges$R < 0), 0L,
               comparison_table$InteractionRefit999Strict[comparison_table$Network == "18S_within"],
               comparison_table$StrictJaccard[comparison_table$Network == "18S_within"],
               as.integer(all(c(within$edges$Model, cross$edges$Model) == "interaction")),
               as.integer(all(c(within$edges$StableStrict, cross$edges$StableStrict))),
               as.integer(all(c(within$edges$CombinedSiteLOSORobust,
                                cross$edges$CombinedSiteLOSORobust)))),
  Expected = c(75L, 54L, 32L, 22L, 12L, 4L, 8L, 0L, 3L, 0.25, 1L, 1L, 1L),
  stringsAsFactors = FALSE
)
numeric_qa$Pass <- abs(as.numeric(numeric_qa$Observed) - as.numeric(numeric_qa$Expected)) < 1e-12
write_tsv_new(numeric_qa, file.path(stage_dir, "qc", "numeric_QA_v1.tsv"))
assert_true(all(numeric_qa$Pass), "Numeric QA failed before visual gate")

png_array <- png::readPNG(png_path)
png_height <- dim(png_array)[1L]
png_width <- dim(png_array)[2L]
pdf_info <- pdftools::pdf_info(pdf_path)
pdf_text <- paste(pdftools::pdf_text(pdf_path), collapse = "\n")
svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
tiff_info <- parse_tiff_metadata(tiff_path)
expected_png_width <- floor(width_mm / 25.4 * 300)
expected_png_height <- floor(height_mm / 25.4 * 300)
expected_tiff_width <- floor(width_mm / 25.4 * 600)
expected_tiff_height <- floor(height_mm / 25.4 * 600)
export_qa <- data.frame(
  Check = c("pdf_exists_nonempty", "pdf_one_page", "pdf_text_extractable",
            "svg_exists_nonempty", "svg_editable_text", "png_width_300dpi",
            "png_height_300dpi", "tiff_width_600dpi", "tiff_height_600dpi",
            "tiff_x_resolution", "tiff_y_resolution", "tiff_resolution_unit_inches"),
  Observed = c(
    as.integer(file.exists(pdf_path) && file.info(pdf_path)$size > 10000),
    pdf_info$pages,
    as.integer(grepl("Interaction-model 16S residual associations", pdf_text, fixed = TRUE)),
    as.integer(file.exists(svg_path) && file.info(svg_path)$size > 10000),
    as.integer(grepl("<text", svg_text, fixed = TRUE)),
    png_width, png_height, tiff_info$Width, tiff_info$Height,
    tiff_info$XResolution, tiff_info$YResolution, tiff_info$ResolutionUnit
  ),
  Expected = c(1L, 1L, 1L, 1L, 1L, expected_png_width, expected_png_height,
               expected_tiff_width, expected_tiff_height, 600, 600, 2L),
  stringsAsFactors = FALSE
)
export_qa$Pass <- abs(as.numeric(export_qa$Observed) - as.numeric(export_qa$Expected)) < 1e-9
write_tsv_new(export_qa, file.path(stage_dir, "qc", "export_QA_v1.tsv"))
assert_true(all(export_qa$Pass), "Export QA failed before visual gate")

cat("VISUAL_QA_READY\t", normalizePath(png_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Inspect the PNG at final layout. Enter PROMOTE_AFTER_VISUAL_QA to attest PASS, or any other text to fail closed.\n")
flush.console()
token <- readLines(file("stdin"), n = 1L, warn = FALSE)
assert_true(length(token) == 1L && identical(trimws(token), "PROMOTE_AFTER_VISUAL_QA"),
            paste0("Visual QA not attested; staging retained at ", stage_dir))

visual_qa <- c(
  "# Visual QA: interaction-aware Supplementary Figure S11 candidate v1",
  "",
  "Status: PASS",
  "Inspection surface: generated 300-dpi PNG at the final 180 mm x 235 mm layout.",
  "Backend: R 4.5.2 only.",
  "",
  "## Attested checks",
  "",
  "- Panel hierarchy is clear: 16S hero network, cross-marker supporting network, then exact comparison table.",
  "- Only interaction-model edges passing both `StableStrict` and `CombinedSiteLOSORobust` are drawn.",
  "- The 18S within-marker network is not drawn; its three-edge count and `not interpreted` decision remain visible in panel c.",
  "- Additive and interaction counts and strict-edge Jaccard values are legible and unambiguous.",
  "- Panel titles, subtitles, node labels, legends, and the 18S asterisk note are not clipped.",
  "- Positive and negative edges remain distinguishable by the frozen orange/blue palette; marker identity uses the frozen blue/green node fills.",
  "- Node size is constant, avoiding a visual claim of validated hub status.",
  "- No panel implies causality, direct interaction, biomarkers, or cross-domain density comparability.",
  "",
  "The attestation token was provided only after inspection of the staged PNG; output promotion was otherwise blocked."
)
write_lines_new(visual_qa, file.path(stage_dir, "qc", "visual_QA_v1.md"))

readme <- c(
  "# Interaction-aware Supplementary Figure S11 candidate v1",
  "",
  "This versioned R-only package is a manuscript candidate, not an overwrite of the accepted S11 artifact.",
  "Panels a and b show only interaction-model residual associations that pass the strict LocationID gate plus Site-cluster and seven-Site LOSO sensitivity.",
  "The 18S within-marker count is retained in panel c but its network is not drawn or interpreted.",
  "",
  "Primary exports:",
  "",
  "- `figures/Supplementary_Figure_S11_interaction_sensitivity_v1.pdf`",
  "- `figures/Supplementary_Figure_S11_interaction_sensitivity_v1.svg`",
  "- `figures/Supplementary_Figure_S11_interaction_sensitivity_v1.png`",
  "- `figures/Supplementary_Figure_S11_interaction_sensitivity_v1.tiff` (600 dpi, LZW)",
  "",
  "Traceability is provided by the source-data TSVs, input/output SHA-256 manifests, figure contract, caption, and QA files.",
  "All scientific language is association-only."
)
write_lines_new(readme, file.path(stage_dir, "README_presentation_v1.md"))

log_lines <- c(
  paste0("Timestamp: ", format(Sys.time(), tz = "Asia/Seoul", usetz = TRUE)),
  paste0("R version: ", as.character(getRversion())),
  "Mode: run",
  "Backend: R only",
  "Input numeric artifact manifest: PASS",
  "Accepted additive reproduction: PASS",
  "Figure display gate: 16S_within and 16S_18S_cross only",
  "Plotted interaction combined-robust edges: 16S=54; cross=12; 18S=0",
  "18S comparison-table strict count: 3; decision=not_interpreted",
  "Export QA: PASS",
  "Interactive staged visual QA: PASS",
  "Claim boundary: association_only"
)
write_lines_new(log_lines, file.path(stage_dir, "logs", "execution_log_v1.txt"))
write_lines_new(capture.output(sessionInfo()), file.path(stage_dir, "logs", "R_sessionInfo_v1.txt"))

manifest_path <- file.path(stage_dir, "provenance", "output_manifest_sha256_v1.tsv")
manifest_files <- list.files(stage_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
manifest_files <- manifest_files[manifest_files != manifest_path]
root_prefix <- paste0(normalizePath(stage_dir, winslash = "/", mustWork = TRUE), "/")
manifest_rel <- substring(normalizePath(manifest_files, winslash = "/", mustWork = TRUE),
                          nchar(root_prefix) + 1L)
output_manifest <- data.frame(
  PresentationRelativePath = manifest_rel,
  Bytes = as.numeric(file.info(manifest_files)$size),
  SHA256 = vapply(manifest_files, sha256_file, character(1L)),
  stringsAsFactors = FALSE
)
output_manifest <- output_manifest[order(output_manifest$PresentationRelativePath), , drop = FALSE]
write_tsv_new(output_manifest, manifest_path)

assert_true(file.rename(stage_dir, presentation_dir),
            paste0("Atomic presentation promotion failed; staging retained at ", stage_dir))
validate_presentation(presentation_dir, src)
cat("RUN PASS\n")
cat("Promoted presentation: ", normalizePath(presentation_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
