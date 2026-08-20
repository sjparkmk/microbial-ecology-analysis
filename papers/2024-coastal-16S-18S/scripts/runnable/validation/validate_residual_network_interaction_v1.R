#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)
package_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

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
mode <- sub("^--", "", modes[[1]])
if (mode == "run") stop("The independent validator is read-only; use --check", call. = FALSE)

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
  invisible(TRUE)
}

read_tsv <- function(path) {
  assert_true(file.exists(path), paste0("Missing file: ", path))
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "")
}

assert_true(as.character(getRversion()) == "4.5.2", paste0("R 4.5.2 required; found ", getRversion()))
assert_true(requireNamespace("digest", quietly = TRUE), "R package digest is required")
sha256_file <- function(path) toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))

source_contract <- read_tsv(file.path(package_root, "config", "protected_source_hashes_v1.tsv"))
source_paths <- setNames(file.path(project_root, source_contract$ProjectRelativePath), source_contract$Role)
for (i in seq_len(nrow(source_contract))) {
  path <- source_paths[[source_contract$Role[[i]]]]
  assert_true(file.exists(path), paste0("Protected source missing: ", path))
  assert_true(identical(sha256_file(path), toupper(source_contract$SHA256[[i]])),
              paste0("Protected source hash drift: ", source_contract$Role[[i]]))
}

artifacts <- file.path(package_root, "artifacts")
if (mode == "dry-run") {
  cat("INDEPENDENT VALIDATOR DRY RUN PASS: R and protected source hashes verified; artifacts present=",
      dir.exists(artifacts), "\n", sep = "")
  quit(save = "no", status = 0)
}

assert_true(dir.exists(artifacts), "Promoted artifacts are missing")
manifest <- read_tsv(file.path(artifacts, "provenance", "output_manifest_sha256_v1.tsv"))
for (i in seq_len(nrow(manifest))) {
  path <- file.path(package_root, manifest$PackageRelativePath[[i]])
  assert_true(file.exists(path), paste0("Manifest target missing: ", manifest$PackageRelativePath[[i]]))
  assert_true(identical(sha256_file(path), toupper(manifest$SHA256[[i]])),
              paste0("Manifest hash mismatch: ", manifest$PackageRelativePath[[i]]))
}

load_count <- function(path) {
  x <- read_tsv(path)
  ids <- as.character(x[[1]])
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"; rownames(m) <- ids
  m
}

taxonomy_id_column <- function(x) {
  hit <- intersect(c("OTUID", "ASV_ID", "#OTUID", "FeatureID"), names(x))
  if (length(hit)) hit[[1]] else names(x)[[1]]
}

aggregate_rank <- function(m, tax, rank) {
  id_col <- taxonomy_id_column(tax)
  idx <- match(rownames(m), tax[[id_col]])
  assert_true(!anyNA(idx), paste0("Independent taxonomy linkage failed: ", rank))
  lab <- trimws(as.character(tax[[rank]][idx]))
  bad <- is.na(lab) | lab == "" | tolower(lab) %in% c("na", "unassigned", "unknown", "uncultured", "unclassified")
  lab[bad] <- "Unassigned"
  rowsum(m, lab, reorder = FALSE)
}

design <- read_tsv(source_paths[["Design"]])
assert_true(nrow(design) == 75L && !anyDuplicated(design$SampleID), "Independent design contract failed")
design$Site <- factor(design$Site, levels = c("Dokdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan", "Ulleungdo"))
design$Season <- factor(design$Season, levels = c("Spring", "Autumn"))
design$SampleType <- factor(design$SampleType, levels = c("Sea-sand", "Water"))

prepare_marker <- function(count_role, tax_role, rank, marker, depth) {
  m <- load_count(source_paths[[count_role]])
  assert_true(setequal(colnames(m), design$SampleID), paste0(marker, " SampleID mismatch"))
  m <- m[, design$SampleID, drop = FALSE]
  assert_true(all(colSums(m) == depth), paste0(marker, " depth mismatch"))
  agg <- aggregate_rank(m, read_tsv(source_paths[[tax_role]]), rank)
  rel <- sweep(agg, 2L, colSums(agg), "/")
  keep <- rownames(agg) != "Unassigned" & rowMeans(agg > 0) >= 0.20 & rowMeans(rel) >= 0.0005
  x <- agg[keep, , drop = FALSE]
  if (nrow(x) > 50L) {
    ord <- order(rowMeans(rel[rownames(x), , drop = FALSE]), decreasing = TRUE)
    x <- x[ord[seq_len(50L)], , drop = FALSE]
  }
  z <- log(t(x) + 0.5); clr <- sweep(z, 1L, rowMeans(z), "-")
  list(marker = marker, rank = rank, counts = x, rel = rel, clr = clr,
       feature = data.frame(Marker = marker, Rank = rank, Taxon = rownames(x),
                            Prevalence = rowMeans(x > 0),
                            MeanRelativeAbundance = rowMeans(rel[rownames(x), , drop = FALSE]),
                            stringsAsFactors = FALSE))
}

m16 <- prepare_marker("16S_rarefied_count", "16S_taxonomy", "Genus", "16S", 20000)
m18 <- prepare_marker("18S_rarefied_count", "18S_PR2_taxonomy", "Class", "18S", 50000)
features <- rbind(m16$feature, m18$feature)
saved_features <- read_tsv(file.path(artifacts, "results", "feature_contract_v1.tsv"))
assert_true(identical(paste(features$Marker, features$Taxon), paste(saved_features$Marker, saved_features$Taxon)),
            "Independent feature identity/order mismatch")
assert_true(max(abs(features$Prevalence - saved_features$Prevalence)) < 1e-12, "Independent prevalence mismatch")
assert_true(max(abs(features$MeanRelativeAbundance - saved_features$MeanRelativeAbundance)) < 1e-12,
            "Independent mean-relative-abundance mismatch")

formulas <- list(additive = ~ Site + Season + SampleType,
                 interaction = ~ Site + SampleType * Season)

fit_model <- function(d, c16, c18, formula) {
  dd <- data.frame(Site = droplevels(factor(as.character(d$Site))),
                   Season = droplevels(factor(as.character(d$Season), levels = c("Spring", "Autumn"))),
                   SampleType = droplevels(factor(as.character(d$SampleType), levels = c("Sea-sand", "Water"))))
  mm <- model.matrix(formula, dd); q <- qr(mm)
  if (q$rank != ncol(mm) || nrow(mm) - q$rank <= 2L) return(list(ok = FALSE, rank = q$rank, columns = ncol(mm)))
  r16 <- qr.resid(q, as.matrix(c16)); r18 <- qr.resid(q, as.matrix(c18))
  colnames(r16) <- colnames(c16); colnames(r18) <- colnames(c18)
  list(ok = TRUE, rank = q$rank, columns = ncol(mm), residual16 = r16, residual18 = r18)
}

edge_templates <- function() {
  within <- function(taxa, marker, network) {
    tmp <- matrix(0, length(taxa), length(taxa)); idx <- which(upper.tri(tmp), arr.ind = TRUE)
    data.frame(Network = network, Marker1 = marker, Marker2 = marker,
               Node1 = taxa[idx[, 1]], Node2 = taxa[idx[, 2]], I = idx[, 1], J = idx[, 2], stringsAsFactors = FALSE)
  }
  cross <- function(a, b) {
    idx <- which(matrix(TRUE, length(a), length(b)), arr.ind = TRUE)
    data.frame(Network = "16S_18S_cross", Marker1 = "16S", Marker2 = "18S",
               Node1 = a[idx[, 1]], Node2 = b[idx[, 2]], I = idx[, 1], J = idx[, 2], stringsAsFactors = FALSE)
  }
  list(`16S_within` = within(colnames(m16$clr), "16S", "16S_within"),
       `18S_within` = within(colnames(m18$clr), "18S", "18S_within"),
       `16S_18S_cross` = cross(colnames(m16$clr), colnames(m18$clr)))
}

templates <- edge_templates()
extract_r <- function(fit, network, template) {
  cm <- if (network == "16S_within") cor(fit$residual16) else if (network == "18S_within") cor(fit$residual18) else cor(fit$residual16, fit$residual18)
  cm[cbind(template$I, template$J)]
}
cor_p <- function(r, n) {
  tt <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  2 * pt(abs(tt), df = n - 2, lower.tail = FALSE)
}

saved_edges <- read_tsv(file.path(artifacts, "results", "network_all_edges_v1.tsv"))
full_fits <- list()
max_full_diff <- 0
for (model in names(formulas)) {
  fit <- fit_model(design, m16$clr, m18$clr, formulas[[model]])
  assert_true(fit$ok, paste0("Independent full model unestimable: ", model))
  full_fits[[model]] <- fit
  for (network in names(templates)) {
    template <- templates[[network]]
    r <- extract_r(fit, network, template); p <- cor_p(r, 75L); q <- p.adjust(p, "BH")
    saved <- saved_edges[saved_edges$Model == model & saved_edges$Network == network, , drop = FALSE]
    assert_true(identical(template$Node1, saved$Node1) && identical(template$Node2, saved$Node2),
                paste0("Independent edge order mismatch: ", model, "/", network))
    max_full_diff <- max(max_full_diff, abs(r - saved$R), abs(p - saved$P), abs(q - saved$Q_BH))
  }
}
assert_true(max_full_diff < 1e-10, paste0("Independent full-edge numeric mismatch: ", max_full_diff))

strict_recalc <- saved_edges$Q_BH < 0.01 & abs(saved_edges$R) >= 0.70 &
  saved_edges$LocationSameSignFraction >= 0.95 & saved_edges$LocationValidFraction >= 0.95 &
  ((saved_edges$LocationLower95 > 0 & saved_edges$R > 0) |
     (saved_edges$LocationUpper95 < 0 & saved_edges$R < 0))
assert_true(identical(as.logical(strict_recalc), as.logical(saved_edges$StableStrict)),
            "Independent strict-gate reconstruction mismatch")

rebuild_draw_matrix <- function(path, unit, expected_positions) {
  x <- read_tsv(path)
  assert_true(nrow(x) == 999L * expected_positions, paste0(unit, " draw-ledger row count mismatch"))
  assert_true(all(table(x$Bootstrap) == expected_positions), paste0(unit, " draw size mismatch"))
  x <- x[order(x$Bootstrap, x$DrawPosition), , drop = FALSE]
  matrix(x$ClusterID, nrow = 999L, byrow = TRUE)
}

location_draws <- rebuild_draw_matrix(file.path(artifacts, "source_data", "location_bootstrap_draws_v1.tsv"), "LocationID", 21L)
site_draws <- rebuild_draw_matrix(file.path(artifacts, "source_data", "site_bootstrap_draws_v1.tsv"), "Site", 7L)
assert_true(all(location_draws %in% unique(design$LocationID)), "Unknown LocationID in draw ledger")
assert_true(all(site_draws %in% unique(as.character(design$Site))), "Unknown Site in draw ledger")

diagnostics <- read_tsv(file.path(artifacts, "results", "bootstrap_fit_diagnostics_v1.tsv"))
assert_true(all(table(diagnostics$ClusterUnit, diagnostics$Model) == 999L), "Bootstrap diagnostic count mismatch")
assert_true(all(aggregate(Estimable ~ ClusterUnit + Model, diagnostics, mean)$Estimable >= 0.95),
            "Bootstrap valid fraction below contract")

audit_bootstrap <- function(unit, draws, model) {
  clusters <- if (unit == "LocationID") as.character(design$LocationID) else as.character(design$Site)
  audit_rows <- saved_edges[saved_edges$Model == model & (saved_edges$StableStrict |
    ave(seq_len(nrow(saved_edges)), saved_edges$Model, saved_edges$Network,
        FUN = function(z) seq_along(z)) <= 3L), , drop = FALSE]
  audit_rows <- audit_rows[!duplicated(paste(audit_rows$Network, audit_rows$Node1, audit_rows$Node2)), , drop = FALSE]
  values <- matrix(NA_real_, nrow(audit_rows), 999L)
  feature_index16 <- setNames(seq_len(ncol(m16$clr)), colnames(m16$clr))
  feature_index18 <- setNames(seq_len(ncol(m18$clr)), colnames(m18$clr))
  for (b in seq_len(999L)) {
    idx <- unlist(lapply(draws[b, ], function(z) which(clusters == z)), use.names = FALSE)
    fit <- tryCatch(fit_model(design[idx, , drop = FALSE], m16$clr[idx, , drop = FALSE],
                              m18$clr[idx, , drop = FALSE], formulas[[model]]), error = function(e) list(ok = FALSE))
    if (!isTRUE(fit$ok)) next
    for (i in seq_len(nrow(audit_rows))) {
      edge <- audit_rows[i, ]
      if (edge$Network == "16S_within") {
        values[i, b] <- cor(fit$residual16[, feature_index16[[edge$Node1]]],
                            fit$residual16[, feature_index16[[edge$Node2]]])
      } else if (edge$Network == "18S_within") {
        values[i, b] <- cor(fit$residual18[, feature_index18[[edge$Node1]]],
                            fit$residual18[, feature_index18[[edge$Node2]]])
      } else {
        values[i, b] <- cor(fit$residual16[, feature_index16[[edge$Node1]]],
                            fit$residual18[, feature_index18[[edge$Node2]]])
      }
    }
  }
  prefix <- if (unit == "LocationID") "Location" else "Site"
  lower <- apply(values, 1L, quantile, 0.025, na.rm = TRUE, names = FALSE)
  median <- apply(values, 1L, quantile, 0.5, na.rm = TRUE, names = FALSE)
  upper <- apply(values, 1L, quantile, 0.975, na.rm = TRUE, names = FALSE)
  same <- vapply(seq_len(nrow(values)), function(i) {
    ok <- is.finite(values[i, ]); mean(sign(values[i, ok]) == sign(audit_rows$R[[i]]))
  }, numeric(1L))
  valid <- rowMeans(is.finite(values))
  diffs <- c(abs(lower - audit_rows[[paste0(prefix, "Lower95")]]),
             abs(median - audit_rows[[paste0(prefix, "Median")]]),
             abs(upper - audit_rows[[paste0(prefix, "Upper95")]]),
             abs(same - audit_rows[[paste0(prefix, "SameSignFraction")]]),
             abs(valid - audit_rows[[paste0(prefix, "ValidFraction")]]))
  max(diffs, na.rm = TRUE)
}

bootstrap_audit_diffs <- c()
for (unit in c("LocationID", "Site")) for (model in names(formulas)) {
  draws <- if (unit == "LocationID") location_draws else site_draws
  value <- audit_bootstrap(unit, draws, model)
  bootstrap_audit_diffs <- c(bootstrap_audit_diffs, value)
  cat(sprintf("Independent %s/%s bootstrap audit max difference: %.3g\n", unit, model, value))
}
assert_true(max(bootstrap_audit_diffs) < 1e-10, "Independent bootstrap audit mismatch")

saved_loso <- read_tsv(file.path(artifacts, "source_data", "loso_all_edges_v1.tsv"))
max_loso_diff <- 0
for (site in levels(droplevels(design$Site))) {
  use <- as.character(design$Site) != site
  for (model in names(formulas)) {
    fit <- fit_model(design[use, , drop = FALSE], m16$clr[use, , drop = FALSE], m18$clr[use, , drop = FALSE], formulas[[model]])
    assert_true(fit$ok, paste0("Independent LOSO unestimable: ", model, "/", site))
    for (network in names(templates)) {
      template <- templates[[network]]
      r <- extract_r(fit, network, template); q <- p.adjust(cor_p(r, sum(use)), "BH")
      saved <- saved_loso[saved_loso$Model == model & saved_loso$Network == network & saved_loso$LeftOutSite == site, , drop = FALSE]
      assert_true(identical(template$Node1, saved$Node1) && identical(template$Node2, saved$Node2),
                  "Independent LOSO edge order mismatch")
      max_loso_diff <- max(max_loso_diff, abs(r - saved$R), abs(q - saved$Q_BH))
    }
  }
}
assert_true(max_loso_diff < 1e-10, paste0("Independent LOSO numeric mismatch: ", max_loso_diff))

baseline <- read_tsv(file.path(artifacts, "results", "accepted_additive_baseline_reproduction_v1.tsv"))
expected <- c(`16S_within` = 57L, `18S_within` = 2L, `16S_18S_cross` = 15L)
observed <- setNames(baseline$ReproducedStableEdges, baseline$Network)
assert_true(all(baseline$Status == "PASS") && identical(as.integer(observed[names(expected)]), as.integer(expected)),
            "Independent accepted-baseline status/count check failed")

summary <- read_tsv(file.path(artifacts, "results", "network_summary_v1.tsv"))
for (i in seq_len(nrow(summary))) {
  hit <- saved_edges$Model == summary$Model[[i]] & saved_edges$Network == summary$Network[[i]]
  assert_true(sum(saved_edges$StableStrict[hit]) == summary$StrictEdges[[i]], "Independent summary count mismatch")
  assert_true(sum(saved_edges$CombinedSiteLOSORobust[hit]) == summary$CombinedSiteLOSORobustStrictEdges[[i]],
              "Independent robust-edge summary mismatch")
}

module_status <- read_tsv(file.path(artifacts, "results", "module_status_v1.tsv"))
assert_true(module_status$Status[[1]] == "PASS", "Module did not fail closed")
figure_status <- read_tsv(file.path(artifacts, "results", "figure_generation_status_v1.tsv"))
assert_true(!as.logical(figure_status$Generated[[1]]), "Unexpected unvalidated figure generation")

cat(sprintf(paste0("INDEPENDENT CHECK PASS: %d manifest entries; full-edge max diff %.3g; ",
                   "bootstrap-audit max diff %.3g; LOSO max diff %.3g.\n"),
            nrow(manifest), max_full_diff, max(bootstrap_audit_diffs), max_loso_diff))
