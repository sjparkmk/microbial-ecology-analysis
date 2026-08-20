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

assert_true(as.character(getRversion()) == "4.5.2", paste0("R 4.5.2 required; found ", getRversion()))
assert_true(requireNamespace("digest", quietly = TRUE), "R package digest is required")

sha256_file <- function(path) toupper(digest::digest(file = path, algo = "sha256", serialize = FALSE))

source_contract_path <- file.path(package_root, "config", "protected_source_hashes_v1.tsv")
network_contract_path <- file.path(package_root, "config", "network_contract_v1.tsv")
claim_contract_path <- file.path(package_root, "config", "claim_gate_v1.md")
source_contract <- read_tsv(source_contract_path)
network_contract <- read_tsv(network_contract_path)

source_paths <- setNames(file.path(project_root, source_contract$ProjectRelativePath), source_contract$Role)

validate_sources <- function() {
  rows <- lapply(seq_len(nrow(source_contract)), function(i) {
    path <- source_paths[[source_contract$Role[[i]]]]
    assert_true(file.exists(path), paste0("Protected source missing: ", path))
    observed <- sha256_file(path)
    expected <- toupper(source_contract$SHA256[[i]])
    assert_true(identical(observed, expected), paste0("Protected source SHA-256 drift: ", source_contract$Role[[i]]))
    info <- file.info(path)
    data.frame(Role = source_contract$Role[[i]], ProjectRelativePath = source_contract$ProjectRelativePath[[i]],
               Bytes = as.numeric(info$size), SHA256 = observed, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

load_count_matrix <- function(path) {
  x <- read_tsv(path)
  assert_true(ncol(x) >= 2L, paste0("Count table has fewer than two columns: ", path))
  ids <- as.character(x[[1]])
  assert_true(!anyNA(ids) && !anyDuplicated(ids), paste0("Invalid feature IDs: ", path))
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  assert_true(all(is.finite(mat)) && all(mat >= 0) && all(mat == floor(mat)),
              paste0("Counts must be finite nonnegative integers: ", path))
  mat
}

load_design <- function(path) {
  x <- read_tsv(path)
  need <- c("SampleID", "Site", "SampleType", "Season", "LocationID")
  miss <- setdiff(need, names(x))
  assert_true(!length(miss), paste0("Design fields missing: ", paste(miss, collapse = ", ")))
  assert_true(nrow(x) == 75L && !anyDuplicated(x$SampleID), "Design must contain 75 unique samples")
  site_levels <- c("Dokdo", "Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan", "Ulleungdo")
  x$Site <- factor(x$Site, levels = site_levels)
  x$SampleType <- factor(x$SampleType, levels = c("Sea-sand", "Water"))
  x$Season <- factor(x$Season, levels = c("Spring", "Autumn"))
  assert_true(!anyNA(x$Site) && !anyNA(x$SampleType) && !anyNA(x$Season), "Unexpected design labels")
  assert_true(length(unique(x$LocationID)) == 21L, "Expected 21 fixed LocationID clusters")
  assert_true(length(unique(x$Site)) == 7L, "Expected seven Sites")
  x
}

taxonomy_id_column <- function(x) {
  hits <- intersect(c("OTUID", "ASV_ID", "#OTUID", "FeatureID"), names(x))
  if (length(hits)) hits[[1]] else names(x)[[1]]
}

normalize_tax_label <- function(x) {
  x <- trimws(as.character(x))
  bad <- is.na(x) | x == "" | tolower(x) %in% c("na", "unassigned", "unknown", "uncultured", "unclassified")
  x[bad] <- "Unassigned"
  x
}

aggregate_counts_by_taxon <- function(mat, taxonomy, rank) {
  id_col <- taxonomy_id_column(taxonomy)
  assert_true(rank %in% names(taxonomy), paste0("Taxonomy rank missing: ", rank))
  idx <- match(rownames(mat), taxonomy[[id_col]])
  assert_true(!anyNA(idx), paste0("Taxonomy linkage failure for ", rank))
  labels <- normalize_tax_label(taxonomy[[rank]][idx])
  out <- rowsum(mat, group = labels, reorder = FALSE)
  out[rowSums(out) > 0, , drop = FALSE]
}

relative_abundance <- function(mat) {
  den <- colSums(mat)
  assert_true(all(den > 0), "Zero-depth sample")
  sweep(mat, 2L, den, "/")
}

clr_transform <- function(mat, pseudocount = 0.5) {
  x <- log(t(mat) + pseudocount)
  sweep(x, 1L, rowMeans(x), "-")
}

load_inputs <- function() {
  design <- load_design(source_paths[["Design"]])
  load_marker <- function(marker, count_role, tax_role, rank, depth) {
    mat <- load_count_matrix(source_paths[[count_role]])
    assert_true(setequal(colnames(mat), design$SampleID), paste0(marker, " count/design SampleIDs differ"))
    mat <- mat[, design$SampleID, drop = FALSE]
    assert_true(all(colSums(mat) == depth), paste0(marker, " rarefaction-depth contract failed"))
    taxonomy <- read_tsv(source_paths[[tax_role]])
    agg <- aggregate_counts_by_taxon(mat, taxonomy, rank)
    list(marker = marker, rank = rank, depth = depth, asv = mat, taxonomy = taxonomy, agg = agg)
  }
  markers <- list(
    `16S` = load_marker("16S", "16S_rarefied_count", "16S_taxonomy", "Genus", 20000),
    `18S` = load_marker("18S", "18S_rarefied_count", "18S_PR2_taxonomy", "Class", 50000)
  )
  list(design = design, markers = markers)
}

select_network_features <- function(marker_obj) {
  x <- marker_obj$agg
  rel <- relative_abundance(x)
  keep <- rownames(x) != "Unassigned" & rowMeans(x > 0) >= 0.20 & rowMeans(rel) >= 0.0005
  x <- x[keep, , drop = FALSE]
  if (nrow(x) > 50L) {
    ord <- order(rowMeans(rel[rownames(x), , drop = FALSE]), decreasing = TRUE)
    x <- x[ord[seq_len(50L)], , drop = FALSE]
  }
  clr <- clr_transform(x, pseudocount = 0.5)
  feature_table <- data.frame(
    Marker = marker_obj$marker, Rank = marker_obj$rank, Taxon = rownames(x),
    Prevalence = rowMeans(x > 0),
    MeanRelativeAbundance = rowMeans(rel[rownames(x), , drop = FALSE]),
    SelectionOrder = seq_len(nrow(x)), stringsAsFactors = FALSE
  )
  list(counts = x, clr = clr, features = feature_table)
}

validate_feature_contract <- function(features) {
  observed <- do.call(rbind, lapply(features, `[[`, "features"))
  accepted <- read_tsv(source_paths[["Accepted_retained_taxa"]])
  assert_true(nrow(observed) == 100L && all(table(observed$Marker) == 50L), "Expected 50 taxa per marker")
  key_observed <- paste(observed$Marker, observed$Taxon, sep = "||")
  key_accepted <- paste(accepted$Marker, accepted$Taxon, sep = "||")
  if (!identical(key_observed, key_accepted)) {
    mismatch <- which(key_observed != key_accepted)
    first <- if (length(mismatch)) mismatch[[1]] else NA_integer_
    stop(sprintf("Retained feature identity/order does not reproduce accepted baseline; first mismatch row=%s observed=%s accepted=%s",
                 first, if (is.finite(first)) key_observed[[first]] else "NA",
                 if (is.finite(first)) key_accepted[[first]] else "NA"), call. = FALSE)
  }
  assert_true(max(abs(observed$Prevalence - accepted$Prevalence)) < 1e-12, "Feature prevalence drift")
  assert_true(max(abs(observed$MeanRelativeAbundance - accepted$MeanRelativeAbundance)) < 1e-12,
              "Feature mean-relative-abundance drift")
  observed
}

model_formulas <- list(
  additive = stats::as.formula("~ Site + Season + SampleType"),
  interaction = stats::as.formula("~ Site + SampleType * Season")
)

fit_residual_pair <- function(design, clr16, clr18, formula) {
  result <- tryCatch({
    d <- data.frame(
      Site = droplevels(factor(as.character(design$Site))),
      Season = droplevels(factor(as.character(design$Season), levels = c("Spring", "Autumn"))),
      SampleType = droplevels(factor(as.character(design$SampleType), levels = c("Sea-sand", "Water")))
    )
    mm <- stats::model.matrix(formula, data = d)
    qr_mm <- qr(mm)
    estimable <- qr_mm$rank == ncol(mm) && nrow(mm) - qr_mm$rank > 2L
    if (!estimable) {
      list(ok = FALSE, error = "rank_deficient_or_insufficient_residual_df", matrix = mm,
           rank = qr_mm$rank, columns = ncol(mm), residual_df = nrow(mm) - qr_mm$rank)
    } else {
      r16 <- qr.resid(qr_mm, as.matrix(clr16))
      r18 <- qr.resid(qr_mm, as.matrix(clr18))
      colnames(r16) <- colnames(clr16); colnames(r18) <- colnames(clr18)
      c16 <- suppressWarnings(stats::cor(r16)); c18 <- suppressWarnings(stats::cor(r18))
      ccross <- suppressWarnings(stats::cor(r16, r18))
      finite <- all(is.finite(c16)) && all(is.finite(c18)) && all(is.finite(ccross))
      list(ok = finite, error = if (finite) "" else "nonfinite_correlation", matrix = mm,
           rank = qr_mm$rank, columns = ncol(mm), residual_df = nrow(mm) - qr_mm$rank,
           residual16 = r16, residual18 = r18, cor16 = c16, cor18 = c18, corcross = ccross)
    }
  }, error = function(e) list(ok = FALSE, error = conditionMessage(e), matrix = NULL,
                              rank = NA_integer_, columns = NA_integer_, residual_df = NA_integer_))
  result
}

make_edge_template <- function(taxa1, taxa2 = NULL, marker1, marker2, network) {
  if (is.null(taxa2)) {
    tmp <- matrix(0, length(taxa1), length(taxa1), dimnames = list(taxa1, taxa1))
    idx <- which(upper.tri(tmp), arr.ind = TRUE)
    data.frame(Node1 = taxa1[idx[, 1]], Node2 = taxa1[idx[, 2]],
               Index1 = idx[, 1], Index2 = idx[, 2], Marker1 = marker1, Marker2 = marker2,
               Network = network, stringsAsFactors = FALSE)
  } else {
    idx <- which(matrix(TRUE, length(taxa1), length(taxa2)), arr.ind = TRUE)
    data.frame(Node1 = taxa1[idx[, 1]], Node2 = taxa2[idx[, 2]],
               Index1 = idx[, 1], Index2 = idx[, 2], Marker1 = marker1, Marker2 = marker2,
               Network = network, stringsAsFactors = FALSE)
  }
}

templates_for <- function(features) {
  t16 <- features$`16S`$features$Taxon
  t18 <- features$`18S`$features$Taxon
  list(
    `16S_within` = make_edge_template(t16, marker1 = "16S", marker2 = "16S", network = "16S_within"),
    `18S_within` = make_edge_template(t18, marker1 = "18S", marker2 = "18S", network = "18S_within"),
    `16S_18S_cross` = make_edge_template(t16, t18, marker1 = "16S", marker2 = "18S", network = "16S_18S_cross")
  )
}

cor_vector <- function(fit, network, template) {
  matrix <- switch(network,
    `16S_within` = fit$cor16,
    `18S_within` = fit$cor18,
    `16S_18S_cross` = fit$corcross,
    stop("Unknown network", call. = FALSE)
  )
  matrix[cbind(template$Index1, template$Index2)]
}

cor_p <- function(r, n) {
  tt <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
  2 * stats::pt(abs(tt), df = n - 2, lower.tail = FALSE)
}

full_edge_stats <- function(fit, n, model, templates) {
  rows <- lapply(names(templates), function(network) {
    template <- templates[[network]]
    r <- cor_vector(fit, network, template)
    p <- cor_p(r, n)
    out <- template[, c("Node1", "Node2", "Marker1", "Marker2", "Network")]
    out$Model <- model
    out$R <- r
    out$P <- p
    out$Q_BH <- stats::p.adjust(p, method = "BH")
    out$EdgeKey <- paste(out$Network, out$Node1, out$Node2, sep = "||")
    out
  })
  do.call(rbind, rows)
}

legacy_bootstrap <- function(fit, design, template, network, seed, nboot = 199L) {
  clusters <- as.character(design$LocationID)
  levels <- unique(clusters)
  set.seed(seed)
  out <- matrix(NA_real_, nrow(template), nboot)
  for (b in seq_len(nboot)) {
    draw <- sample(levels, length(levels), replace = TRUE)
    idx <- unlist(lapply(draw, function(z) which(clusters == z)), use.names = FALSE)
    if (network == "16S_within") matrix <- stats::cor(fit$residual16[idx, , drop = FALSE])
    else if (network == "18S_within") matrix <- stats::cor(fit$residual18[idx, , drop = FALSE])
    else matrix <- stats::cor(fit$residual16[idx, , drop = FALSE], fit$residual18[idx, , drop = FALSE])
    out[, b] <- matrix[cbind(template$Index1, template$Index2)]
  }
  out
}

summarize_bootstrap <- function(values, full_r, prefix) {
  valid_n <- rowSums(is.finite(values))
  qfun <- function(x, p) if (sum(is.finite(x)) > 0L) stats::quantile(x, p, na.rm = TRUE, names = FALSE) else NA_real_
  lower <- apply(values, 1L, qfun, p = 0.025)
  median <- apply(values, 1L, qfun, p = 0.5)
  upper <- apply(values, 1L, qfun, p = 0.975)
  same <- vapply(seq_len(nrow(values)), function(i) {
    ok <- is.finite(values[i, ])
    if (!any(ok)) return(NA_real_)
    mean(sign(values[i, ok]) == sign(full_r[[i]]))
  }, numeric(1L))
  out <- data.frame(valid_n, valid_n / ncol(values), lower, median, upper, same)
  names(out) <- paste0(prefix, c("ValidN", "ValidFraction", "Lower95", "Median", "Upper95", "SameSignFraction"))
  out
}

reproduce_accepted_baseline <- function(full_additive, design, templates) {
  seeds <- c(`16S_within` = 130001L, `18S_within` = 130002L, `16S_18S_cross` = 130003L)
  accepted <- read_tsv(source_paths[["Accepted_all_edges"]])
  accepted$EdgeKey <- paste(accepted$Network, accepted$Node1, accepted$Node2, sep = "||")
  rows <- list(); edge_rows <- list()
  for (network in names(templates)) {
    template <- templates[[network]]
    full <- full_additive[full_additive$Network == network, , drop = FALSE]
    boot <- legacy_bootstrap(full_fit_cache$additive, design, template, network, seeds[[network]], 199L)
    bs <- summarize_bootstrap(boot, full$R, "Legacy")
    calc <- cbind(full, bs)
    calc$ScreenPass <- calc$Q_BH < 0.05 & abs(calc$R) >= 0.40 & calc$LegacySameSignFraction >= 0.80 &
      ((calc$LegacyLower95 > 0 & calc$R > 0) | (calc$LegacyUpper95 < 0 & calc$R < 0))
    calc$Stable <- calc$Q_BH < 0.01 & abs(calc$R) >= 0.70 & calc$LegacySameSignFraction >= 0.95 &
      ((calc$LegacyLower95 > 0 & calc$R > 0) | (calc$LegacyUpper95 < 0 & calc$R < 0))
    acc <- accepted[match(calc$EdgeKey, accepted$EdgeKey), , drop = FALSE]
    assert_true(!anyNA(acc$EdgeKey), paste0("Accepted edge linkage failed: ", network))
    numeric_pairs <- list(R = c("R", "R"), P = c("P", "P"), Q = c("Q_BH", "Q_BH"),
                          Lower = c("LegacyLower95", "BootstrapLower95"),
                          Median = c("LegacyMedian", "BootstrapMedian"),
                          Upper = c("LegacyUpper95", "BootstrapUpper95"),
                          SameSign = c("LegacySameSignFraction", "SameSignFraction"))
    diffs <- vapply(numeric_pairs, function(z) max(abs(calc[[z[[1]]]] - acc[[z[[2]]]]), na.rm = TRUE), numeric(1L))
    stable_mismatch <- sum(as.logical(calc$Stable) != as.logical(acc$Stable))
    screen_mismatch <- sum(as.logical(calc$ScreenPass) != as.logical(acc$ScreenPass))
    rows[[network]] <- data.frame(
      Network = network, TestedEdges = nrow(calc), AcceptedStableEdges = sum(as.logical(acc$Stable)),
      ReproducedStableEdges = sum(calc$Stable), StableIdentityMismatches = stable_mismatch,
      ScreenIdentityMismatches = screen_mismatch, MaxAbsRDiff = diffs[["R"]],
      MaxAbsPDiff = diffs[["P"]], MaxAbsQDiff = diffs[["Q"]],
      MaxAbsLowerDiff = diffs[["Lower"]], MaxAbsMedianDiff = diffs[["Median"]],
      MaxAbsUpperDiff = diffs[["Upper"]], MaxAbsSameSignDiff = diffs[["SameSign"]],
      Status = if (stable_mismatch == 0L && screen_mismatch == 0L && max(diffs) < 1e-10) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
    edge_rows[[network]] <- calc
  }
  summary <- do.call(rbind, rows)
  assert_true(all(summary$Status == "PASS"), "Accepted additive baseline did not reproduce")
  expected <- c(`16S_within` = 57L, `18S_within` = 2L, `16S_18S_cross` = 15L)
  got <- setNames(summary$ReproducedStableEdges, summary$Network)
  assert_true(identical(as.integer(got[names(expected)]), as.integer(expected)), "Accepted strict-edge counts did not reproduce")
  list(summary = summary, edges = do.call(rbind, edge_rows))
}

make_draw_matrix <- function(levels, nboot, seed) {
  set.seed(seed)
  matrix(sample(levels, nboot * length(levels), replace = TRUE), nrow = nboot, byrow = TRUE)
}

draw_ledger <- function(draws, unit) {
  data.frame(Bootstrap = rep(seq_len(nrow(draws)), each = ncol(draws)),
             DrawPosition = rep(seq_len(ncol(draws)), times = nrow(draws)),
             ClusterUnit = unit, ClusterID = as.vector(t(draws)), stringsAsFactors = FALSE)
}

run_refit_bootstrap <- function(design, features, templates, unit, draws, log_fun) {
  clusters <- if (unit == "LocationID") as.character(design$LocationID) else as.character(design$Site)
  nboot <- nrow(draws)
  stores <- lapply(names(model_formulas), function(model) {
    setNames(lapply(templates, function(template) matrix(NA_real_, nrow(template), nboot)), names(templates))
  })
  names(stores) <- names(model_formulas)
  diagnostics <- vector("list", nboot * length(model_formulas)); k <- 0L
  for (b in seq_len(nboot)) {
    idx <- unlist(lapply(draws[b, ], function(z) which(clusters == z)), use.names = FALSE)
    db <- design[idx, , drop = FALSE]
    c16 <- features$`16S`$clr[idx, , drop = FALSE]
    c18 <- features$`18S`$clr[idx, , drop = FALSE]
    for (model in names(model_formulas)) {
      k <- k + 1L
      fit <- fit_residual_pair(db, c16, c18, model_formulas[[model]])
      diagnostics[[k]] <- data.frame(
        ClusterUnit = unit, Bootstrap = b, Model = model, Samples = length(idx),
        UniqueSites = length(unique(as.character(db$Site))),
        UniqueLocations = length(unique(as.character(db$LocationID))),
        DesignColumns = fit$columns, DesignRank = fit$rank, ResidualDf = fit$residual_df,
        Estimable = isTRUE(fit$ok), Error = fit$error, stringsAsFactors = FALSE
      )
      if (isTRUE(fit$ok)) {
        for (network in names(templates)) {
          stores[[model]][[network]][, b] <- cor_vector(fit, network, templates[[network]])
        }
      }
    }
    if (b %% 100L == 0L || b == nboot) log_fun(sprintf("%s refit bootstrap %d/%d", unit, b, nboot))
  }
  list(values = stores, diagnostics = do.call(rbind, diagnostics))
}

run_loso <- function(design, features, templates, full_edges) {
  sites <- levels(droplevels(design$Site))
  stores <- lapply(names(model_formulas), function(model) {
    setNames(lapply(templates, function(template) list(
      R = matrix(NA_real_, nrow(template), length(sites), dimnames = list(NULL, sites)),
      Q = matrix(NA_real_, nrow(template), length(sites), dimnames = list(NULL, sites)))), names(templates))
  })
  names(stores) <- names(model_formulas)
  diagnostics <- list(); k <- 0L
  for (s in sites) {
    use <- as.character(design$Site) != s
    for (model in names(model_formulas)) {
      fit <- fit_residual_pair(design[use, , drop = FALSE], features$`16S`$clr[use, , drop = FALSE],
                               features$`18S`$clr[use, , drop = FALSE], model_formulas[[model]])
      k <- k + 1L
      diagnostics[[k]] <- data.frame(LeftOutSite = s, Model = model, Samples = sum(use),
        DesignColumns = fit$columns, DesignRank = fit$rank, ResidualDf = fit$residual_df,
        Estimable = isTRUE(fit$ok), Error = fit$error, stringsAsFactors = FALSE)
      if (isTRUE(fit$ok)) {
        for (network in names(templates)) {
          r <- cor_vector(fit, network, templates[[network]])
          q <- stats::p.adjust(cor_p(r, sum(use)), method = "BH")
          stores[[model]][[network]]$R[, s] <- r
          stores[[model]][[network]]$Q[, s] <- q
        }
      }
    }
  }
  long_rows <- list(); summary_rows <- list(); k <- 0L
  for (model in names(model_formulas)) for (network in names(templates)) {
    k <- k + 1L
    template <- templates[[network]]
    rmat <- stores[[model]][[network]]$R
    qmat <- stores[[model]][[network]]$Q
    full <- full_edges[full_edges$Model == model & full_edges$Network == network, , drop = FALSE]
    long_rows[[k]] <- data.frame(
      Model = model, Network = network,
      Node1 = rep(template$Node1, times = length(sites)), Node2 = rep(template$Node2, times = length(sites)),
      Marker1 = rep(template$Marker1, times = length(sites)), Marker2 = rep(template$Marker2, times = length(sites)),
      LeftOutSite = rep(sites, each = nrow(template)),
      R = as.vector(rmat), Q_BH = as.vector(qmat),
      CoreMagnitudeQPass = as.vector(abs(rmat) >= 0.70 & qmat < 0.01), stringsAsFactors = FALSE)
    valid_n <- rowSums(is.finite(rmat))
    same <- vapply(seq_len(nrow(rmat)), function(i) {
      ok <- is.finite(rmat[i, ])
      if (!any(ok)) NA_real_ else mean(sign(rmat[i, ok]) == sign(full$R[[i]]))
    }, numeric(1L))
    core_fraction <- rowMeans(abs(rmat) >= 0.70 & qmat < 0.01, na.rm = TRUE)
    min_abs <- apply(abs(rmat), 1L, function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_)
    median_r <- apply(rmat, 1L, function(x) if (any(is.finite(x))) median(x, na.rm = TRUE) else NA_real_)
    max_delta <- vapply(seq_len(nrow(rmat)), function(i) {
      x <- rmat[i, ]; if (!any(is.finite(x))) NA_real_ else max(abs(x - full$R[[i]]), na.rm = TRUE)
    }, numeric(1L))
    summary_rows[[k]] <- data.frame(
      Model = model, Network = network, Node1 = template$Node1, Node2 = template$Node2,
      LOSOValidN = valid_n, LOSOValidFraction = valid_n / length(sites),
      LOSOSameSignFraction = same, LOSOCoreMagnitudeQFraction = core_fraction,
      LOSOMinAbsR = min_abs, LOSOMedianR = median_r, LOSOMaxAbsDeltaR = max_delta,
      stringsAsFactors = FALSE)
  }
  list(long = do.call(rbind, long_rows), summary = do.call(rbind, summary_rows),
       diagnostics = do.call(rbind, diagnostics))
}

edge_degree_correlation <- function(a, b, template) {
  nodes <- unique(c(paste(template$Marker1, template$Node1, sep = "::"),
                    paste(template$Marker2, template$Node2, sep = "::")))
  degree_vec <- function(edges) {
    out <- setNames(rep(0, length(nodes)), nodes)
    if (nrow(edges)) {
      ends <- c(paste(edges$Marker1, edges$Node1, sep = "::"), paste(edges$Marker2, edges$Node2, sep = "::"))
      tab <- table(ends)
      out[names(tab)] <- as.numeric(tab)
    }
    out
  }
  da <- degree_vec(a); db <- degree_vec(b)
  if (stats::sd(da) == 0 || stats::sd(db) == 0) return(NA_real_)
  suppressWarnings(stats::cor(da, db, method = "spearman"))
}

compare_models <- function(all_edges, templates) {
  rows <- lapply(names(templates), function(network) {
    a <- all_edges[all_edges$Model == "additive" & all_edges$Network == network, , drop = FALSE]
    b <- all_edges[all_edges$Model == "interaction" & all_edges$Network == network, , drop = FALSE]
    assert_true(identical(a$EdgeKey, b$EdgeKey), paste0("Model edge order drift: ", network))
    sa <- a[a$StableStrict, , drop = FALSE]; sb <- b[b$StableStrict, , drop = FALSE]
    overlap <- intersect(sa$EdgeKey, sb$EdgeKey); union <- union(sa$EdgeKey, sb$EdgeKey)
    overlap_sign <- if (length(overlap)) {
      aa <- sa$R[match(overlap, sa$EdgeKey)]; bb <- sb$R[match(overlap, sb$EdgeKey)]
      mean(sign(aa) == sign(bb))
    } else NA_real_
    data.frame(
      Network = network, AdditiveStrictEdges = nrow(sa), InteractionStrictEdges = nrow(sb),
      StrictOverlap = length(overlap), StrictUnion = length(union),
      StrictJaccard = if (length(union)) length(overlap) / length(union) else NA_real_,
      OverlapSignConcordance = overlap_sign,
      AllEdgeRPearson = stats::cor(a$R, b$R, method = "pearson"),
      AllEdgeRSpearman = stats::cor(a$R, b$R, method = "spearman"),
      MedianAbsRChange = stats::median(abs(a$R - b$R)),
      MaxAbsRChange = max(abs(a$R - b$R)),
      StrictDegreeRankCorrelation = edge_degree_correlation(sa, sb, templates[[network]]),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

safe_atomic_stage <- function(stage) {
  expected_parent <- normalizePath(package_root, winslash = "/", mustWork = TRUE)
  stage_parent <- normalizePath(dirname(stage), winslash = "/", mustWork = TRUE)
  assert_true(identical(expected_parent, stage_parent), "Unsafe staging parent")
  assert_true(grepl("^\\._artifacts_staging_[0-9]+$", basename(stage)), "Unsafe staging basename")
}

verify_output_manifest <- function() {
  artifacts <- file.path(package_root, "artifacts")
  manifest_path <- file.path(artifacts, "provenance", "output_manifest_sha256_v1.tsv")
  manifest <- read_tsv(manifest_path)
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(package_root, manifest$PackageRelativePath[[i]])
    assert_true(file.exists(path), paste0("Manifest file missing: ", manifest$PackageRelativePath[[i]]))
    assert_true(identical(sha256_file(path), toupper(manifest$SHA256[[i]])),
                paste0("Output manifest hash mismatch: ", manifest$PackageRelativePath[[i]]))
  }
  manifest
}

internal_check <- function() {
  validate_sources()
  manifest <- verify_output_manifest()
  artifacts <- file.path(package_root, "artifacts")
  baseline <- read_tsv(file.path(artifacts, "results", "accepted_additive_baseline_reproduction_v1.tsv"))
  assert_true(all(baseline$Status == "PASS"), "Baseline reproduction check failed")
  edges <- read_tsv(file.path(artifacts, "results", "network_all_edges_v1.tsv"))
  summary <- read_tsv(file.path(artifacts, "results", "network_summary_v1.tsv"))
  diagnostics <- read_tsv(file.path(artifacts, "results", "bootstrap_fit_diagnostics_v1.tsv"))
  loso_diag <- read_tsv(file.path(artifacts, "results", "loso_fit_diagnostics_v1.tsv"))
  assert_true(nrow(edges) == 9900L, "Expected 9,900 model-network edge rows")
  edge_counts <- table(edges$Model, edges$Network)
  assert_true(all(edge_counts[, "16S_within"] == 1225L) && all(edge_counts[, "18S_within"] == 1225L) &&
                all(edge_counts[, "16S_18S_cross"] == 2500L), "Edge family size mismatch")
  strict_recalc <- edges$Q_BH < 0.01 & abs(edges$R) >= 0.70 &
    edges$LocationSameSignFraction >= 0.95 & edges$LocationValidFraction >= 0.95 &
    ((edges$LocationLower95 > 0 & edges$R > 0) | (edges$LocationUpper95 < 0 & edges$R < 0))
  assert_true(identical(as.logical(edges$StableStrict), as.logical(strict_recalc)), "Strict-edge gate mismatch")
  for (i in seq_len(nrow(summary))) {
    hit <- edges$Model == summary$Model[[i]] & edges$Network == summary$Network[[i]]
    assert_true(sum(edges$StableStrict[hit]) == summary$StrictEdges[[i]], "Network summary strict count mismatch")
  }
  assert_true(all(table(diagnostics$ClusterUnit, diagnostics$Model) == 999L), "Bootstrap diagnostic count mismatch")
  assert_true(nrow(loso_diag) == 14L && all(table(loso_diag$Model) == 7L), "LOSO diagnostic count mismatch")
  status <- read_tsv(file.path(artifacts, "results", "module_status_v1.tsv"))
  assert_true(status$Status[[1]] == "PASS", "Module status is not PASS")
  cat(sprintf("CHECK PASS: %d manifested files; %d edges; baseline and gates verified.\n", nrow(manifest), nrow(edges)))
  invisible(TRUE)
}

if (mode == "check") {
  assert_true(dir.exists(file.path(package_root, "artifacts")), "No promoted artifacts to check")
  internal_check()
  quit(save = "no", status = 0)
}

source_hash_start <- validate_sources()
inputs <- load_inputs()
features <- lapply(inputs$markers, select_network_features)
feature_contract <- validate_feature_contract(features)
templates <- templates_for(features)

full_fit_cache <- lapply(names(model_formulas), function(model) {
  fit <- fit_residual_pair(inputs$design, features$`16S`$clr, features$`18S`$clr, model_formulas[[model]])
  assert_true(isTRUE(fit$ok), paste0("Full ", model, " model is not estimable: ", fit$error))
  fit
})
names(full_fit_cache) <- names(model_formulas)

full_edge_rows <- lapply(names(model_formulas), function(model) {
  full_edge_stats(full_fit_cache[[model]], nrow(inputs$design), model, templates)
})
full_edges <- do.call(rbind, full_edge_rows)

if (mode == "dry-run") {
  assert_true(!dir.exists(file.path(package_root, "artifacts")), "artifacts/ already exists; run would refuse overwrite")
  accepted_summary <- read_tsv(source_paths[["Accepted_network_summary"]])
  assert_true(setequal(accepted_summary$StableEdges, c(57L, 2L, 15L)), "Accepted baseline summary changed")
  cat("DRY RUN PASS\n")
  cat("  R:", as.character(getRversion()), "\n")
  cat("  Samples/Sites/LocationIDs:", nrow(inputs$design), length(unique(inputs$design$Site)),
      length(unique(inputs$design$LocationID)), "\n")
  cat("  Retained taxa: 50 16S genera and 50 PR2 18S classes\n")
  cat("  Full additive and interaction design matrices are estimable\n")
  cat("  Protected source SHA-256 values match the frozen contract\n")
  quit(save = "no", status = 0)
}

artifacts_final <- file.path(package_root, "artifacts")
assert_true(!file.exists(artifacts_final), "Refusing overwrite of existing artifacts/")
stage <- file.path(package_root, paste0("._artifacts_staging_", Sys.getpid()))
assert_true(!file.exists(stage), "Staging path already exists")
dir.create(file.path(stage, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "provenance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "source_data"), recursive = TRUE, showWarnings = FALSE)
safe_atomic_stage(stage)
promoted <- FALSE
on.exit({
  if (!promoted && dir.exists(stage)) {
    safe_atomic_stage(stage)
    unlink(stage, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

log_lines <- character()
log_message <- function(text) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text)
  cat(line, "\n")
  log_lines <<- c(log_lines, line)
}

log_message("Starting accepted additive-baseline reproduction")
baseline <- reproduce_accepted_baseline(full_edges[full_edges$Model == "additive", , drop = FALSE],
                                        inputs$design, templates)
log_message("Accepted additive baseline reproduced")

location_levels <- unique(as.character(inputs$design$LocationID))
site_levels <- unique(as.character(inputs$design$Site))
location_draws <- make_draw_matrix(location_levels, 999L, 2026081401L)
site_draws <- make_draw_matrix(site_levels, 999L, 2026081402L)
assert_true(nrow(unique(as.data.frame(location_draws))) > 990L, "Unexpectedly duplicated LocationID draw realizations")
assert_true(nrow(unique(as.data.frame(site_draws))) > 900L, "Unexpectedly duplicated Site draw realizations")

log_message("Starting 999 LocationID-cluster refit bootstraps")
location_boot <- run_refit_bootstrap(inputs$design, features, templates, "LocationID", location_draws, log_message)
log_message("Starting 999 Site-cluster refit sensitivity bootstraps")
site_boot <- run_refit_bootstrap(inputs$design, features, templates, "Site", site_draws, log_message)

bootstrap_diagnostics <- rbind(location_boot$diagnostics, site_boot$diagnostics)
valid_fraction <- aggregate(Estimable ~ ClusterUnit + Model, bootstrap_diagnostics, mean)
assert_true(all(valid_fraction$Estimable >= 0.95), "A refit-bootstrap model had <95% estimable fits")

combined_edges <- list(); k <- 0L
for (model in names(model_formulas)) for (network in names(templates)) {
  k <- k + 1L
  full <- full_edges[full_edges$Model == model & full_edges$Network == network, , drop = FALSE]
  loc <- summarize_bootstrap(location_boot$values[[model]][[network]], full$R, "Location")
  site <- summarize_bootstrap(site_boot$values[[model]][[network]], full$R, "Site")
  out <- cbind(full, loc, site)
  out$StableStrict <- out$Q_BH < 0.01 & abs(out$R) >= 0.70 &
    out$LocationSameSignFraction >= 0.95 & out$LocationValidFraction >= 0.95 &
    ((out$LocationLower95 > 0 & out$R > 0) | (out$LocationUpper95 < 0 & out$R < 0))
  out$SiteSensitivityPass <- out$Q_BH < 0.01 & abs(out$R) >= 0.70 &
    out$SiteSameSignFraction >= 0.95 & out$SiteValidFraction >= 0.95 &
    ((out$SiteLower95 > 0 & out$R > 0) | (out$SiteUpper95 < 0 & out$R < 0))
  combined_edges[[k]] <- out
}
all_edges <- do.call(rbind, combined_edges)

log_message("Starting seven leave-one-Site-out residual-model/network refits")
loso <- run_loso(inputs$design, features, templates, full_edges)
loso_key <- paste(loso$summary$Model, loso$summary$Network, loso$summary$Node1, loso$summary$Node2, sep = "||")
edge_loso_key <- paste(all_edges$Model, all_edges$Network, all_edges$Node1, all_edges$Node2, sep = "||")
idx_loso <- match(edge_loso_key, loso_key)
assert_true(!anyNA(idx_loso), "LOSO summary linkage failed")
all_edges <- cbind(all_edges, loso$summary[idx_loso, setdiff(names(loso$summary), c("Model", "Network", "Node1", "Node2")), drop = FALSE])
all_edges$CombinedSiteLOSORobust <- all_edges$StableStrict & all_edges$SiteSensitivityPass &
  all_edges$LOSOValidFraction == 1 & all_edges$LOSOSameSignFraction >= (6 / 7)

network_summary_rows <- list(); k <- 0L
for (model in names(model_formulas)) for (network in names(templates)) {
  k <- k + 1L
  x <- all_edges[all_edges$Model == model & all_edges$Network == network, , drop = FALSE]
  strict <- x$StableStrict
  robust_n <- sum(x$CombinedSiteLOSORobust)
  network_summary_rows[[k]] <- data.frame(
    Model = model, Network = network, Marker1 = x$Marker1[[1]], Marker2 = x$Marker2[[1]],
    TestedEdges = nrow(x), BHSignificantEdges = sum(x$Q_BH < 0.01),
    MagnitudeEdges = sum(abs(x$R) >= 0.70), StrictEdges = sum(strict),
    StrictPositive = sum(strict & x$R > 0), StrictNegative = sum(strict & x$R < 0),
    StrictDensity = sum(strict) / nrow(x), SiteSensitivityPassStrictEdges = sum(strict & x$SiteSensitivityPass),
    CombinedSiteLOSORobustStrictEdges = robust_n,
    CombinedRobustFraction = if (sum(strict)) robust_n / sum(strict) else NA_real_,
    LocationBootstrapValidFraction = valid_fraction$Estimable[valid_fraction$ClusterUnit == "LocationID" & valid_fraction$Model == model],
    SiteBootstrapValidFraction = valid_fraction$Estimable[valid_fraction$ClusterUnit == "Site" & valid_fraction$Model == model],
    stringsAsFactors = FALSE)
}
network_summary <- do.call(rbind, network_summary_rows)
model_comparison <- compare_models(all_edges, templates)

baseline_counts <- setNames(baseline$summary$ReproducedStableEdges, baseline$summary$Network)
claim_rows <- lapply(c("16S_within", "18S_within", "16S_18S_cross"), function(network) {
  x <- network_summary[network_summary$Model == "interaction" & network_summary$Network == network, , drop = FALSE]
  reproducible <- all(baseline$summary$Status == "PASS")
  fit_ok <- x$LocationBootstrapValidFraction >= 0.95 && x$SiteBootstrapValidFraction >= 0.95
  count_ok <- x$StrictEdges >= 3L && x$StrictEdges <= 60L
  density_ok <- x$StrictDensity <= 0.05
  stability_ok <- is.finite(x$CombinedRobustFraction) && x$CombinedRobustFraction >= 0.80
  pass <- reproducible && fit_ok && count_ok && density_ok && stability_ok
  status <- if (pass) "retain_SI_association_only" else if (network == "18S_within" || x$StrictEdges < 3L) "not_informative" else "table_only_or_remove"
  data.frame(
    Network = network, AcceptedBaselineStrictEdges = baseline_counts[[network]],
    InteractionStrictEdges = x$StrictEdges, InteractionStrictDensity = x$StrictDensity,
    CombinedRobustFraction = x$CombinedRobustFraction,
    BaselineReproductionPass = reproducible, BootstrapFitPass = fit_ok,
    EdgeCountPass = count_ok, DensityPass = density_ok, SiteLOSOStabilityPass = stability_ok,
    Status = status,
    AllowedClaim = if (pass) "hypothesis-generating residual association in Supplementary Information only" else "no positive network claim",
    ProhibitedClaims = "direct interaction; hub; bridge; causality; biomarker; cross-marker density comparison",
    stringsAsFactors = FALSE)
})
claim_gate <- do.call(rbind, claim_rows)

module_status <- data.frame(
  Module = "residual_network_interaction_v1", Status = "PASS",
  BaselineReproduction = "PASS", FullModelsEstimable = TRUE,
  LocationBootstrapMinimumValidFraction = min(valid_fraction$Estimable[valid_fraction$ClusterUnit == "LocationID"]),
  SiteBootstrapMinimumValidFraction = min(valid_fraction$Estimable[valid_fraction$ClusterUnit == "Site"]),
  ClaimBoundary = "association_only", stringsAsFactors = FALSE)

model_diag <- do.call(rbind, lapply(names(full_fit_cache), function(model) {
  fit <- full_fit_cache[[model]]
  data.frame(Model = model, Formula = paste(deparse(model_formulas[[model]]), collapse = " "),
             Samples = nrow(inputs$design), DesignColumns = fit$columns, DesignRank = fit$rank,
             ResidualDf = fit$residual_df, Estimable = fit$ok, Error = fit$error,
             stringsAsFactors = FALSE)
}))

write_tsv_new(feature_contract, file.path(stage, "results", "feature_contract_v1.tsv"))
write_tsv_new(model_diag, file.path(stage, "results", "full_model_diagnostics_v1.tsv"))
write_tsv_new(baseline$summary, file.path(stage, "results", "accepted_additive_baseline_reproduction_v1.tsv"))
write_tsv_new(baseline$edges, file.path(stage, "source_data", "accepted_additive_baseline_reproduced_edges_v1.tsv"))
write_tsv_new(all_edges, file.path(stage, "results", "network_all_edges_v1.tsv"))
write_tsv_new(all_edges[all_edges$StableStrict, , drop = FALSE], file.path(stage, "results", "network_strict_edges_v1.tsv"))
write_tsv_new(network_summary, file.path(stage, "results", "network_summary_v1.tsv"))
write_tsv_new(model_comparison, file.path(stage, "results", "additive_interaction_comparison_v1.tsv"))
write_tsv_new(claim_gate, file.path(stage, "results", "network_claim_gate_v1.tsv"))
write_tsv_new(module_status, file.path(stage, "results", "module_status_v1.tsv"))
write_tsv_new(bootstrap_diagnostics, file.path(stage, "results", "bootstrap_fit_diagnostics_v1.tsv"))
write_tsv_new(draw_ledger(location_draws, "LocationID"), file.path(stage, "source_data", "location_bootstrap_draws_v1.tsv"))
write_tsv_new(draw_ledger(site_draws, "Site"), file.path(stage, "source_data", "site_bootstrap_draws_v1.tsv"))
write_tsv_new(loso$long, file.path(stage, "source_data", "loso_all_edges_v1.tsv"))
write_tsv_new(loso$summary, file.path(stage, "results", "loso_edge_stability_v1.tsv"))
write_tsv_new(loso$diagnostics, file.path(stage, "results", "loso_fit_diagnostics_v1.tsv"))

figure_supported <- all(claim_gate$Status[match(c("16S_within", "16S_18S_cross"), claim_gate$Network)] == "retain_SI_association_only")
figure_status <- data.frame(Figure = "interaction_residual_networks", Generated = FALSE,
                            Reason = if (figure_supported) "claim_gate_supported_but_figure_not_required_for_numeric_sensitivity" else "claim_gate_did_not_support_replacement_figure",
                            stringsAsFactors = FALSE)
write_tsv_new(figure_status, file.path(stage, "results", "figure_generation_status_v1.tsv"))

source_hash_end <- validate_sources()
source_manifest <- merge(source_hash_start, source_hash_end, by = c("Role", "ProjectRelativePath", "Bytes"),
                         suffixes = c("_Start", "_End"), sort = FALSE)
source_manifest$Unchanged <- source_manifest$SHA256_Start == source_manifest$SHA256_End
assert_true(all(source_manifest$Unchanged), "Protected source changed during analysis")
write_tsv_new(source_manifest, file.path(stage, "provenance", "protected_source_manifest_v1.tsv"))

write_lines_new(capture.output(sessionInfo()), file.path(stage, "logs", "R_sessionInfo_v1.txt"))
log_message("Analysis and claim gates completed")
write_lines_new(log_lines, file.path(stage, "logs", "execution_log_v1.txt"))

result_readme <- c(
  "# Residual-association interaction sensitivity results v1", "",
  "The protected accepted additive baseline was reproduced before any sensitivity interpretation.", "",
  "The additive model is `CLR ~ Site + Season + SampleType`; the interaction model is `CLR ~ Site + SampleType * Season`.",
  "Both use the same 75 samples, 50 16S genera, 50 PR2-primary 18S classes, pseudocount 0.5, CLR transformation, network-wide BH families, and frozen strict edge gate.", "",
  "LocationID and Site cluster bootstraps each used 999 saved draws and refitted the residual model in every draw. Seven leave-one-Site-out model/network refits were also completed.", "",
  "Use `network_claim_gate_v1.tsv` before any manuscript recommendation. All edges remain hypothesis-generating residual associations and cannot establish direct biological interaction, hubs, bridges, causality, or biomarkers.", "",
  paste0("Figure generated: no. Reason: ", figure_status$Reason[[1]], ".")
)
write_lines_new(result_readme, file.path(stage, "README_results_v1.md"))

package_files <- c(file.path(package_root, "README.md"), source_contract_path, network_contract_path,
                   claim_contract_path, script_file, file.path(script_dir, "validate_residual_network_interaction_v1.R"))
assert_true(all(file.exists(package_files)), "Package source file missing before manifest creation")
stage_files <- list.files(stage, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
stage_files <- stage_files[file.info(stage_files)$isdir %in% FALSE]
stage_rel <- substring(normalizePath(stage_files, winslash = "/", mustWork = TRUE), nchar(normalizePath(stage, winslash = "/", mustWork = TRUE)) + 2L)
package_rel <- substring(normalizePath(package_files, winslash = "/", mustWork = TRUE), nchar(normalizePath(package_root, winslash = "/", mustWork = TRUE)) + 2L)
output_manifest <- rbind(
  data.frame(PackageRelativePath = file.path("artifacts", stage_rel), Bytes = file.info(stage_files)$size,
             SHA256 = vapply(stage_files, sha256_file, character(1L)), stringsAsFactors = FALSE),
  data.frame(PackageRelativePath = package_rel, Bytes = file.info(package_files)$size,
             SHA256 = vapply(package_files, sha256_file, character(1L)), stringsAsFactors = FALSE)
)
output_manifest$PackageRelativePath <- gsub("\\\\", "/", output_manifest$PackageRelativePath)
write_tsv_new(output_manifest, file.path(stage, "provenance", "output_manifest_sha256_v1.tsv"))

if (!file.rename(stage, artifacts_final)) stop("Atomic artifacts promotion failed", call. = FALSE)
promoted <- TRUE
cat("RUN PASS: promoted immutable artifacts to ", normalizePath(artifacts_final, winslash = "/"), "\n", sep = "")
internal_check()
