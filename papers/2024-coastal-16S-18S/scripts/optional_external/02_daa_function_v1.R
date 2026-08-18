#!/usr/bin/env Rscript

# Original-pipeline-aligned differential-abundance and functional module.
#
# Statistical calculations are rebuilt on the curated 75-sample inputs.  The
# original EasyAmplicon scripts define module order and visual grammar only;
# legacy numeric results are never read here.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(ALDEx2)
  library(ANCOMBC)
  library(ggrepel)
})

SCRIPT_VERSION <- "1.2.0"
SEED <- 1L
ANCOMBC_WORKERS <- 6L

args <- commandArgs(trailingOnly = TRUE)
mode_flags <- intersect(args, c("--dry-run", "--run", "--check"))
if (length(mode_flags) != 1L) {
  stop("Choose exactly one mode: --dry-run, --run, or --check", call. = FALSE)
}
mode <- sub("^--", "", mode_flags[[1]])

arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}

project_root <- normalizePath(arg_value("--project-root", getwd()), winslash = "/", mustWork = TRUE)
output_root <- arg_value(
  "--output-root",
  file.path(project_root, "revision", "re_analysis", "12_original_pipeline_reanalysis_v1", "daa_function_v2")
)
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
script_dir <- file.path(project_root, "revision", "re_analysis", "original_style_pipeline_v1")
source(file.path(script_dir, "common_original_style.R"), local = TRUE)

paths <- list(
  design = file.path(project_root, "revision/re_analysis/00_design_v4/master_sample_metadata_v4.tsv"),
  count16 = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_asv_counts_analysis_75_v1.tsv"),
  tax16 = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_taxonomy_analysis_75_v1.tsv"),
  count18 = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_asv_counts_pr2_analysis_v1.tsv"),
  tax18 = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_taxonomy_pr2_primary_analysis_v1.tsv"),
  flags18 = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/qc/18s_asv_flags_ledger_v2.tsv"),
  fap_db = "/your/data/FAPROTAX.txt", # insert your path
  fap_upstream = "/your/data/collapse_table.py", # insert your path
  fap_wrapper = file.path(script_dir, "faprotax_compat_runner_v1.py")
)

required_outputs <- c(
  "results/daa_all_methods_v2.tsv",
  "results/daa_consensus_taxa_v2.tsv",
  "results/daa_surface_summary_v2.tsv",
  "results/faprotax_function_counts_v2.tsv",
  "results/faprotax_mapping_coverage_v2.tsv",
  "results/faprotax_contrasts_v2.tsv",
  "results/faprotax_18s_not_applicable_v2.tsv",
  "qc/module_status_daa_function_v2.tsv",
  "provenance/input_sha256_v2.tsv",
  "provenance/analysis_contract_v2.tsv",
  "provenance/r_session_info_v2.txt"
)

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("R package digest is required", call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

load_inputs <- function() {
  invisible(lapply(paths[c("design", "count16", "tax16", "count18", "tax18", "flags18", "fap_db", "fap_upstream", "fap_wrapper")], function(p) {
    assert_true(file.exists(p), paste0("Missing input: ", p))
  }))
  design <- load_design(paths$design)
  c16 <- align_counts_design(load_count_matrix(paths$count16)$matrix, design)$matrix
  c18 <- align_counts_design(load_count_matrix(paths$count18)$matrix, design)$matrix
  t16 <- read_tsv_strict(paths$tax16)
  t18 <- read_tsv_strict(paths$tax18)
  f18 <- read_tsv_strict(paths$flags18)
  assert_true(setequal(rownames(c16), t16[[taxonomy_id_column(t16)]]), "16S count/taxonomy feature sets differ")
  assert_true(setequal(rownames(c18), t18[[taxonomy_id_column(t18)]]), "18S count/taxonomy feature sets differ")
  fkeep <- f18[tolower(as.character(f18$analysis_keep)) == "true", , drop = FALSE]
  assert_true(setequal(rownames(c18), fkeep$ASV_ID), "18S count/flag analysis feature sets differ")
  list(design = design, c16 = c16, c18 = c18, t16 = t16, t18 = t18, f18 = fkeep)
}

aggregate_dataset <- function(mat, taxonomy, rank, dataset, marker, omit_ids = character()) {
  if (length(omit_ids)) mat <- mat[!rownames(mat) %in% omit_ids, , drop = FALSE]
  agg <- aggregate_counts_by_taxon(mat, taxonomy, rank)
  list(matrix = agg, dataset = dataset, marker = marker, rank = rank)
}

prepare_datasets <- function(inp) {
  metazoa <- inp$f18$ASV_ID[tolower(as.character(inp$f18$flag_pr2_metazoa)) == "true"]
  list(
    aggregate_dataset(inp$c16, inp$t16, "Genus", "16S_genus", "16S"),
    aggregate_dataset(inp$c18, inp$t18, "Class", "18S_PR2_class_full", "18S"),
    aggregate_dataset(inp$c18, inp$t18, "Class", "18S_PR2_class_non_metazoan", "18S", metazoa)
  )
}

surface_specs <- function(design) {
  complete_sites <- c("Goseong", "Pohang", "Pyoseon", "Sehwa", "Seongsan")
  list(
    list(
      id = "A72_Substrate", samples = as.character(design$SampleID[design$PairID_Substrate != ""]),
      condition = "SampleType", baseline = "Sea-sand", target = "Water",
      pair = "PairID_Substrate", fix = "PairID_Substrate + SampleType",
      random = NULL, term_tokens = c("SampleType", "Water"),
      original_reference = "ancombc2_multigroup_script.R; ANCOMBC2_plots_and_tables_windows_v6.R"
    ),
    list(
      id = "A36_SeaSand_Season", samples = as.character(design$SampleID[design$SampleType == "Sea-sand" & design$PairID_Season != ""]),
      condition = "Season", baseline = "Spring", target = "Autumn",
      pair = "LocationID", fix = "LocationID + Season", random = NULL,
      term_tokens = c("Season", "Autumn"),
      original_reference = "ancombc2_multigroup_script.R; multigroup_compare_volcano_up11_taxonomy_fix.R"
    ),
    list(
      id = "A30_Water_Season", samples = as.character(design$SampleID[design$SampleType == "Water" & design$PairID_Season != ""]),
      condition = "Season", baseline = "Spring", target = "Autumn",
      pair = "LocationID", fix = "LocationID + Season", random = NULL,
      term_tokens = c("Season", "Autumn"),
      original_reference = "ancombc2_multigroup_script.R; multigroup_compare_volcano_up11_taxonomy_fix.R"
    ),
    list(
      id = "A60_Interaction", samples = as.character(design$SampleID[as.character(design$Site) %in% complete_sites]),
      condition = "interaction", baseline = "Sea-sand_Spring", target = "Water_Autumn interaction",
      pair = "LocationID", fix = "LocationID + SampleType * Season", random = NULL,
      term_tokens = c("SampleType", "Water", "Season", "Autumn"),
      original_reference = "ancombc2_multigroup_script.R; ANCOMBC2_plots_and_tables_windows_v6.R"
    )
  )
}

set_surface_factors <- function(md, spec) {
  md$Site <- factor(as.character(md$Site))
  md$SiteSeasonBlock <- factor(as.character(md$SiteSeasonBlock))
  md$LocationID <- factor(as.character(md$LocationID))
  md$PairID_Substrate <- factor(as.character(md$PairID_Substrate))
  md$SampleType <- factor(as.character(md$SampleType), levels = c("Sea-sand", "Water"))
  md$Season <- factor(as.character(md$Season), levels = c("Spring", "Autumn"))
  rownames(md) <- md$SampleID
  md
}

prevalence_filter <- function(mat, min_fraction = 0.10) {
  keep <- rowMeans(mat > 0) >= min_fraction & rowSums(mat) > 0
  mat[keep, , drop = FALSE]
}

find_term_column <- function(nms, prefix, tokens) {
  cand <- nms[startsWith(nms, prefix)]
  for (tok in tokens) cand <- cand[grepl(tok, cand, fixed = TRUE)]
  if (length(cand) != 1L) return(NA_character_)
  cand[[1]]
}

run_ancombc <- function(mat, md, spec) {
  set.seed(SEED)
  fit <- tryCatch(
    ANCOMBC::ancombc2(
      data = mat, taxa_are_rows = TRUE, meta_data = md,
      fix_formula = spec$fix, rand_formula = spec$random,
      p_adj_method = "BH", pseudo = 0, pseudo_sens = FALSE,
      prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
      group = NULL, struc_zero = FALSE, neg_lb = FALSE,
      alpha = 0.05, n_cl = ANCOMBC_WORKERS, verbose = FALSE,
      global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(ok = FALSE, message = conditionMessage(fit), table = data.frame()))
  }
  res <- as.data.frame(fit$res, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"taxon" %in% names(res)) res$taxon <- rownames(res)
  lc <- find_term_column(names(res), "lfc_", spec$term_tokens)
  pc <- find_term_column(names(res), "p_", spec$term_tokens)
  qc <- find_term_column(names(res), "q_", spec$term_tokens)
  dc <- find_term_column(names(res), "diff_", spec$term_tokens)
  if (anyNA(c(lc, pc, qc))) {
    return(list(ok = FALSE, message = paste0("Target ANCOM-BC2 coefficient not uniquely found. Columns: ", paste(names(res), collapse = ",")), table = res))
  }
  out <- data.frame(
    Taxon = as.character(res$taxon), Method = "ANCOMBC2",
    Effect = as.numeric(res[[lc]]), PValue = as.numeric(res[[pc]]),
    QValue = as.numeric(res[[qc]]), PassedSensitivity = TRUE,
    MethodDetail = paste0("fix=", spec$fix, "; random=", spec$random %||% "none", "; coefficient=", lc),
    stringsAsFactors = FALSE
  )
  if (!is.na(dc)) out$DeclaredDifferential <- as.logical(res[[dc]]) else out$DeclaredDifferential <- out$QValue < 0.05
  # Official ANCOM-BC2 pseudo-count sensitivity is expensive on the complete
  # taxon set.  Run it only for primary BH-significant taxa; non-significant
  # rows cannot become candidates and therefore do not require the extra fit.
  sig_taxa <- out$Taxon[is.finite(out$QValue) & out$QValue < 0.05]
  if (length(sig_taxa)) {
    sens <- tryCatch(
      ANCOMBC::ancombc2(
        data = mat[sig_taxa, , drop = FALSE], taxa_are_rows = TRUE, meta_data = md,
        fix_formula = spec$fix, rand_formula = spec$random,
        p_adj_method = "BH", pseudo = 0, pseudo_sens = TRUE,
        prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
        group = NULL, struc_zero = FALSE, neg_lb = FALSE,
        alpha = 0.05, n_cl = ANCOMBC_WORKERS, verbose = FALSE,
        global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE
      ),
      error = function(e) e
    )
    out$PassedSensitivity[out$Taxon %in% sig_taxa] <- FALSE
    if (!inherits(sens, "error")) {
      sr <- as.data.frame(sens$res, check.names = FALSE, stringsAsFactors = FALSE)
      if (!"taxon" %in% names(sr)) sr$taxon <- rownames(sr)
      ssc <- find_term_column(names(sr), "passed_ss_", spec$term_tokens)
      if (!is.na(ssc)) {
        ii <- match(out$Taxon, as.character(sr$taxon))
        hit <- out$Taxon %in% sig_taxa & !is.na(ii)
        out$PassedSensitivity[hit] <- as.logical(sr[[ssc]][ii[hit]])
      }
    }
  }
  list(ok = TRUE, message = "ok", table = out)
}

run_deseq2 <- function(mat, md, spec) {
  if (spec$id == "A60_Interaction") {
    design_formula <- ~ LocationID + SampleType + Season + SampleType:Season
    target_tokens <- c("SampleTypeWater", "SeasonAutumn")
  } else if (spec$condition == "SampleType") {
    design_formula <- ~ PairID_Substrate + SampleType
    target_tokens <- c("SampleType", "Water")
  } else {
    design_formula <- ~ LocationID + Season
    target_tokens <- c("Season", "Autumn")
  }
  fit <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(round(mat), colData = md, design = design_formula)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    rn <- DESeq2::resultsNames(dds)
    cand <- rn
    for (tok in target_tokens) cand <- cand[grepl(tok, cand, fixed = TRUE)]
    if (length(cand) != 1L) stop("DESeq2 target coefficient is not unique: ", paste(rn, collapse = ","))
    rr <- as.data.frame(DESeq2::results(dds, name = cand[[1]], alpha = 0.05), stringsAsFactors = FALSE)
    rr$Taxon <- rownames(rr)
    list(ok = TRUE, coefficient = cand[[1]], res = rr)
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (!fit$ok) return(list(ok = FALSE, message = fit$message, table = data.frame()))
  out <- data.frame(
    Taxon = fit$res$Taxon, Method = "DESeq2", Effect = fit$res$log2FoldChange,
    PValue = fit$res$pvalue, QValue = fit$res$padj, PassedSensitivity = TRUE,
    DeclaredDifferential = is.finite(fit$res$padj) & fit$res$padj < 0.05,
    MethodDetail = paste0("paired fixed-effect design; coefficient=", fit$coefficient),
    stringsAsFactors = FALSE
  )
  list(ok = TRUE, message = "ok", table = out)
}

paired_order <- function(md, spec) {
  p <- as.character(md[[spec$pair]])
  cnd <- as.character(md[[spec$condition]])
  tab <- table(p, cnd)
  assert_true(ncol(tab) == 2L && all(tab == 1L), paste0(spec$id, " must contain one sample per pair and condition"))
  ord <- order(p, match(cnd, c(spec$baseline, spec$target)))
  md[ord, , drop = FALSE]
}

run_aldex2 <- function(mat, md, spec) {
  if (spec$id == "A60_Interaction") {
    return(list(ok = FALSE, message = "ALDEx2 has no registered paired interaction estimator", table = data.frame()))
  }
  # ALDEx2 is a confirmatory sensitivity surface.  Restrict it to taxa present
  # in at least 20% of this contrast so Monte-Carlo effort is spent on taxa
  # with enough paired information; ANCOM-BC2 and DESeq2 retain the 10% gate.
  mat <- mat[rowMeans(mat > 0) >= 0.20, , drop = FALSE]
  if (nrow(mat) > 150L) {
    mat <- mat[head(order(rowSums(mat), decreasing = TRUE), 150L), , drop = FALSE]
  }
  if (nrow(mat) < 2L) return(list(ok = FALSE, message = "Fewer than two taxa passed the ALDEx2 20% prevalence sensitivity gate", table = data.frame()))
  md <- paired_order(md, spec)
  mat <- mat[, md$SampleID, drop = FALSE]
  # ALDEx2 1.42.0 accepts a character condition vector here. Passing a factor
  # triggers an internal `round.factor` error before the paired test begins.
  conds <- as.character(md[[spec$condition]])
  fit <- tryCatch({
    set.seed(SEED)
    clr <- ALDEx2::aldex.clr(mat, conds, mc.samples = 32, denom = "all", verbose = FALSE)
    tt <- as.data.frame(ALDEx2::aldex.ttest(clr, paired.test = TRUE, hist.plot = FALSE, verbose = FALSE))
    ef <- as.data.frame(ALDEx2::aldex.effect(clr, verbose = FALSE, paired.test = TRUE))
    tt$Taxon <- rownames(tt); ef$Taxon <- rownames(ef)
    list(ok = TRUE, tt = tt, ef = ef)
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  if (!fit$ok) return(list(ok = FALSE, message = fit$message, table = data.frame()))
  qcol <- intersect(c("we.eBH", "wi.eBH"), names(fit$tt))[[1]]
  pcol <- intersect(c("we.ep", "wi.ep"), names(fit$tt))[[1]]
  ecol <- intersect(c("effect", "diff.btw"), names(fit$ef))[[1]]
  z <- merge(fit$tt[, c("Taxon", pcol, qcol), drop = FALSE], fit$ef[, c("Taxon", ecol), drop = FALSE], by = "Taxon")
  out <- data.frame(
    Taxon = z$Taxon, Method = "ALDEx2", Effect = as.numeric(z[[ecol]]),
    PValue = as.numeric(z[[pcol]]), QValue = as.numeric(z[[qcol]]),
    PassedSensitivity = TRUE,
    DeclaredDifferential = is.finite(z[[qcol]]) & z[[qcol]] < 0.05,
    MethodDetail = paste0("paired Welch/Wilcoxon Monte Carlo sensitivity; prevalence>=20%; top150 by total count if needed; mc.samples=32; effect=", ecol),
    stringsAsFactors = FALSE
  )
  list(ok = TRUE, message = "ok", table = out)
}

simple_effect_and_loso <- function(mat, md, spec) {
  x <- log(mat + 0.5)
  x <- sweep(x, 2, colMeans(x), "-")
  get_effect <- function(keep_sites = NULL) {
    mm <- md
    if (!is.null(keep_sites)) mm <- mm[as.character(mm$Site) %in% keep_sites, , drop = FALSE]
    xx <- x[, mm$SampleID, drop = FALSE]
    if (spec$id == "A60_Interaction") {
      keys <- unique(as.character(mm$LocationID))
      vals <- sapply(keys, function(k) {
        q <- mm[as.character(mm$LocationID) == k, , drop = FALSE]
        if (nrow(q) != 4L) return(rep(NA_real_, nrow(xx)))
        at <- function(tp, ss) {
          sid <- q$SampleID[as.character(q$SampleType) == tp & as.character(q$Season) == ss]
          if (length(sid) != 1L) return(rep(NA_real_, nrow(xx)))
          xx[, sid]
        }
        (at("Water", "Autumn") - at("Water", "Spring")) - (at("Sea-sand", "Autumn") - at("Sea-sand", "Spring"))
      })
    } else {
      keys <- unique(as.character(mm[[spec$pair]]))
      vals <- sapply(keys, function(k) {
        q <- mm[as.character(mm[[spec$pair]]) == k, , drop = FALSE]
        a <- q$SampleID[as.character(q[[spec$condition]]) == spec$target]
        b <- q$SampleID[as.character(q[[spec$condition]]) == spec$baseline]
        if (length(a) != 1L || length(b) != 1L) return(rep(NA_real_, nrow(xx)))
        xx[, a] - xx[, b]
      })
    }
    if (is.null(dim(vals))) vals <- matrix(vals, nrow = nrow(xx))
    apply(vals, 1, stats::median, na.rm = TRUE)
  }
  full <- get_effect()
  sites <- sort(unique(as.character(md$Site)))
  loo <- sapply(sites, function(s) get_effect(setdiff(sites, s)))
  if (is.null(dim(loo))) loo <- matrix(loo, nrow = nrow(mat))
  stable <- rowMeans(sign(loo) == sign(full), na.rm = TRUE)
  data.frame(Taxon = rownames(mat), SimpleCLREffect = full, LOSODirectionStability = stable, stringsAsFactors = FALSE)
}

orient_aldex <- function(tab, simple) {
  if (!nrow(tab)) return(tab)
  z <- merge(tab[, c("Taxon", "Effect")], simple[, c("Taxon", "SimpleCLREffect")], by = "Taxon")
  ok <- is.finite(z$Effect) & is.finite(z$SimpleCLREffect) & z$Effect != 0 & z$SimpleCLREffect != 0
  if (sum(ok) >= 10L && mean(sign(z$Effect[ok]) == sign(z$SimpleCLREffect[ok])) < 0.5) tab$Effect <- -tab$Effect
  tab
}

consensus_table <- function(methods, simple, marker, dataset, rank, spec) {
  if (!nrow(methods)) return(data.frame())
  methods$Significant <- is.finite(methods$QValue) & methods$QValue < 0.05 & methods$PassedSensitivity
  taxa <- sort(unique(methods$Taxon))
  rows <- lapply(taxa, function(tx) {
    z <- methods[methods$Taxon == tx, , drop = FALSE]
    sig <- z[z$Significant, , drop = FALSE]
    primary <- sig$Method == "ANCOMBC2"
    ndir <- length(unique(sign(sig$Effect[is.finite(sig$Effect) & sig$Effect != 0])))
    support <- nrow(sig)
    primary_sig <- any(primary)
    cls <- if (support == 0L) "no_signal" else if (ndir > 1L) "discordant" else if (primary_sig && support >= 3L) "robust_3method" else if (primary_sig && support >= 2L) "supported_2method" else if (primary_sig) "primary_only" else "sensitivity_only"
    data.frame(
      Marker = marker, Dataset = dataset, Rank = rank, SurfaceID = spec$id,
      Taxon = tx, MethodsAvailable = paste(sort(unique(z$Method)), collapse = ";"),
      SignificantMethods = paste(sort(unique(sig$Method)), collapse = ";"),
      NMethodsSignificant = support, ConsensusClass = cls,
      ConsensusDirection = if (support && ndir == 1L) ifelse(sign(sig$Effect[[1]]) > 0, paste0("higher_in_", gsub(" ", "_", spec$target)), paste0("higher_in_", gsub(" ", "_", spec$baseline))) else "none_or_discordant",
      MinQ = if (support) min(sig$QValue, na.rm = TRUE) else min(z$QValue, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows) %>% left_join(simple, by = "Taxon")
  out$StableDirection <- is.finite(out$LOSODirectionStability) & out$LOSODirectionStability >= 0.80
  out$ManuscriptStatus <- ifelse(
    out$ConsensusClass == "robust_3method" & out$StableDirection & out$MinQ < 0.01, "main_candidate",
    ifelse(out$ConsensusClass %in% c("robust_3method", "supported_2method") & out$StableDirection, "supplement_candidate", "not_informative")
  )
  out
}

run_daa <- function(inp, stage) {
  datasets <- prepare_datasets(inp)
  specs <- surface_specs(inp$design)
  all_methods <- list(); all_consensus <- list(); summaries <- list(); status <- list()
  k <- 0L
  for (ds in datasets) for (spec in specs) {
    # The non-metazoan 18S lane is a predeclared sensitivity for the dominant
    # substrate contrast, not a duplicate four-surface discovery analysis.
    if (ds$dataset == "18S_PR2_class_non_metazoan" && spec$id != "A72_Substrate") next
    k <- k + 1L
    md <- inp$design[match(spec$samples, inp$design$SampleID), , drop = FALSE]
    md <- set_surface_factors(md, spec)
    mat <- ds$matrix[, md$SampleID, drop = FALSE]
    before <- nrow(mat)
    mat <- prevalence_filter(mat, 0.10)
    simple <- simple_effect_and_loso(mat, md, spec)
    ac <- run_ancombc(mat, md, spec)
    de <- run_deseq2(mat, md, spec)
    al <- run_aldex2(mat, md, spec)
    if (al$ok) al$table <- orient_aldex(al$table, simple)
    mt <- bind_rows(
      if (ac$ok) ac$table else NULL,
      if (de$ok) de$table else NULL,
      if (al$ok) al$table else NULL
    )
    if (nrow(mt)) {
      mt$Marker <- ds$marker; mt$Dataset <- ds$dataset; mt$Rank <- ds$rank; mt$SurfaceID <- spec$id
      mt <- mt %>% left_join(simple, by = "Taxon")
      mt$Significant <- is.finite(mt$QValue) & mt$QValue < 0.05 & mt$PassedSensitivity
      all_methods[[k]] <- mt
      cc <- consensus_table(mt, simple, ds$marker, ds$dataset, ds$rank, spec)
      all_consensus[[k]] <- cc
      ncand <- sum(cc$ManuscriptStatus %in% c("main_candidate", "supplement_candidate"))
      plot_df <- mt[mt$Method == "ANCOMBC2", , drop = FALSE] %>% left_join(cc[, c("Taxon", "ConsensusClass", "ManuscriptStatus")], by = "Taxon")
      plot_df$NegLog10Q <- -log10(pmax(plot_df$QValue, 1e-300))
      plot_df$Label <- ""
      lab <- head(plot_df$Taxon[order(plot_df$QValue, -abs(plot_df$Effect), na.last = NA)], 10L)
      plot_df$Label[plot_df$Taxon %in% lab & plot_df$QValue < 0.10] <- plot_df$Taxon[plot_df$Taxon %in% lab & plot_df$QValue < 0.10]
      p <- ggplot(plot_df, aes(Effect, NegLog10Q, color = ManuscriptStatus)) +
        geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, color = "grey45") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
        geom_point(size = 1.35, alpha = 0.82) +
        ggrepel::geom_text_repel(aes(label = Label), size = 2, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, show.legend = FALSE) +
        scale_color_manual(values = c(main_candidate = "#B2182B", supplement_candidate = "#EF8A62", not_informative = "#9E9E9E"), drop = FALSE) +
        labs(title = paste(ds$dataset, spec$id), subtitle = "ANCOM-BC2 primary; colors require cross-method and leave-one-site-out support", x = "ANCOM-BC2 log fold change", y = expression(-log[10](BH~q)), color = "Result status") +
        theme_easyamplicon_original(8)
      fig_base <- file.path(stage, "figures", paste0("DA_", ds$dataset, "_", spec$id, "_volcano_v2"))
      save_plot_bundle(p, fig_base, 150, 110)
    } else {
      ncand <- 0L
    }
    messages <- paste(c(if (!ac$ok) paste0("ANCOMBC2: ", ac$message), if (!de$ok) paste0("DESeq2: ", de$message), if (!al$ok) paste0("ALDEx2: ", al$message)), collapse = " | ")
    summaries[[k]] <- data.frame(
      Marker = ds$marker, Dataset = ds$dataset, Rank = ds$rank, SurfaceID = spec$id,
      NSamples = ncol(mat), NIndependentPairsOrLocations = length(unique(as.character(md[[spec$pair]]))),
      TaxaBeforeFilter = before, TaxaAfterPrevalence10Pct = nrow(mat),
      ANCOMBC2_OK = ac$ok, DESeq2_OK = de$ok, ALDEx2_OK = al$ok,
      NCandidateTaxa = ncand, Notes = messages, stringsAsFactors = FALSE
    )
    status[[k]] <- data.frame(
      ModuleID = paste0("DA_", ds$dataset, "_", spec$id), RunStatus = if (ac$ok && de$ok) "completed" else "partial",
      NSamples = ncol(mat), NIndependentUnits = length(unique(as.character(md[[spec$pair]]))),
      Gate = "ANCOM-BC2 BH q<0.05 + passed_ss; >=1 sensitivity same direction; LOSO direction >=80%",
      GateResult = if (ncand > 0L) "pass_candidate_taxa_present" else "no_taxon_passed",
      ManuscriptRole = if (ncand > 0L) "candidate_taxa_table_and_plot" else "not_reportable",
      Reason = if (ncand > 0L) paste(ncand, "taxa passed") else paste0("No robust stable taxon. ", messages),
      OriginalScriptReference = spec$original_reference, stringsAsFactors = FALSE
    )
  }
  methods <- bind_rows(all_methods)
  consensus <- bind_rows(all_consensus)
  summaries <- bind_rows(summaries)
  status <- bind_rows(status)
  write_tsv_atomic(methods, file.path(stage, "results/daa_all_methods_v2.tsv"))
  write_tsv_atomic(consensus, file.path(stage, "results/daa_consensus_taxa_v2.tsv"))
  write_tsv_atomic(summaries, file.path(stage, "results/daa_surface_summary_v2.tsv"))
  list(status = status, methods = methods, consensus = consensus, summary = summaries)
}

build_faprotax_input <- function(inp, path) {
  tx <- inp$t16
  idc <- taxonomy_id_column(tx)
  ranks <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), names(tx))
  z <- tx[match(rownames(inp$c16), tx[[idc]]), , drop = FALSE]
  lin <- apply(z[, ranks, drop = FALSE], 1, function(v) paste(normalize_tax_label(v), collapse = ";"))
  # Row names need not be unique according to collapse_table.py; keeping ASV-level
  # rows preserves exact mapped-read coverage.
  out <- data.frame(taxonomy = lin, inp$c16, check.names = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(out, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
  invisible(path)
}

run_faprotax <- function(inp, stage) {
  raw_dir <- file.path(stage, "source_data/faprotax")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  input <- file.path(raw_dir, "16s_taxonomy_count_input_v2.tsv")
  output <- file.path(raw_dir, "16s_faprotax_collapsed_raw_v2.tsv")
  report <- file.path(raw_dir, "16s_faprotax_assignment_report_v2.txt")
  mapping <- file.path(raw_dir, "16s_faprotax_taxonomy_to_groups_v2.tsv")
  build_faprotax_input(inp, input)
  py <- Sys.which("python")
  assert_true(nzchar(py), "Python executable not found for protected FAPROTAX runner")
  cmd_args <- c(
    paths$fap_wrapper, "-i", input, "-g", paths$fap_db, "-o", output, "-r", report,
    "--out_groups2records_table_dense", mapping,
    "--row_names_are_in_column", "taxonomy", "--column_names_are_in", "first_data_line",
    "--group_members_defined_as", "words", "--non_numeric", "ignore",
    "--group_leftovers_as", "unassigned", "--omit_unrepresented_groups",
    "--output_format_collapsed", "classical"
  )
  log_out <- file.path(stage, "logs/faprotax_stdout_v2.log")
  log_err <- file.path(stage, "logs/faprotax_stderr_v2.log")
  dir.create(dirname(log_out), recursive = TRUE, showWarnings = FALSE)
  rc <- system2(py, args = shQuote(cmd_args), stdout = log_out, stderr = log_err)
  assert_true(identical(as.integer(rc), 0L) && file.exists(output), paste0("FAPROTAX failed; see ", log_err))
  fun <- read.delim(output, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  assert_true(ncol(fun) == 76L, "Unexpected FAPROTAX output schema")
  names(fun)[1] <- "Function"
  sample_cols <- inp$design$SampleID
  assert_true(setequal(names(fun)[-1], sample_cols), "FAPROTAX sample set differs")
  fun <- fun[, c("Function", sample_cols), drop = FALSE]
  for (j in sample_cols) fun[[j]] <- as.numeric(fun[[j]])
  un_idx <- which(tolower(fun$Function) == "unassigned")
  assert_true(length(un_idx) == 1L, "FAPROTAX unassigned row missing or duplicated")
  source_depth <- colSums(inp$c16)
  unassigned <- as.numeric(fun[un_idx, sample_cols, drop = TRUE])
  coverage <- data.frame(
    SampleID = sample_cols, SourceReads = as.numeric(source_depth[sample_cols]),
    FAPROTAXUnassignedReads = unassigned,
    FAPROTAXMappedReads = as.numeric(source_depth[sample_cols]) - unassigned,
    MappedFraction = (as.numeric(source_depth[sample_cols]) - unassigned) / as.numeric(source_depth[sample_cols]),
    stringsAsFactors = FALSE
  ) %>% left_join(as.data.frame(inp$design)[, c("SampleID", "Site", "GeologicalHabitat", "SampleType", "Season", "LocationID", "PairID_Substrate", "PairID_Season")], by = "SampleID")
  fun_assigned <- fun[-un_idx, , drop = FALSE]
  write_tsv_atomic(fun, file.path(stage, "results/faprotax_function_counts_v2.tsv"))
  write_tsv_atomic(coverage, file.path(stage, "results/faprotax_mapping_coverage_v2.tsv"))
  contrasts <- faprotax_contrasts(fun_assigned, inp$design)
  write_tsv_atomic(contrasts, file.path(stage, "results/faprotax_contrasts_v2.tsv"))
  pass_coverage <- stats::median(coverage$MappedFraction) >= 0.30 && min(coverage$MappedFraction) >= 0.20
  cand <- contrasts$QValue < 0.05 & contrasts$LOSODirectionStability >= 0.80
  cand[is.na(cand)] <- FALSE
  if (nrow(fun_assigned)) {
    totals <- rowSums(as.matrix(fun_assigned[, sample_cols, drop = FALSE]))
    keep <- head(order(totals, decreasing = TRUE), 15L)
    plot_fun <- fun_assigned[keep, , drop = FALSE]
    long <- pivot_longer(plot_fun, all_of(sample_cols), names_to = "SampleID", values_to = "Count") %>%
      left_join(as.data.frame(inp$design)[, c("SampleID", "Site", "SampleType", "Season")], by = "SampleID") %>%
      group_by(Function, Site, SampleType, Season) %>% summarise(MeanCount = mean(Count), .groups = "drop")
    long$Function <- factor(long$Function, levels = rev(plot_fun$Function))
    p <- ggplot(long, aes(MeanCount, Function, color = SampleType, shape = Season)) +
      geom_point(alpha = 0.82, size = 1.8) +
      scale_x_log10(labels = scales::label_number()) +
      scale_color_manual(values = sampletype_palette) +
      scale_shape_manual(values = c(Spring = 16, Autumn = 17)) +
      facet_wrap(~Site, ncol = 4, scales = "free_x") +
      labs(title = "16S taxonomy-inferred functional groups (FAPROTAX 1.2.7)", subtitle = "Exploratory only; points are cell mean inferred counts and functions can overlap", x = "Mean inferred count (log10 scale)", y = NULL) +
      theme_easyamplicon_original(7)
    save_plot_bundle(p, file.path(stage, "figures/FAPROTAX_16S_top15_original_bubble_v2"), 220, 150)
  }
  status <- data.frame(
    ModuleID = "FAPROTAX_16S", RunStatus = "completed", NSamples = 75L, NIndependentUnits = 21L,
    Gate = "median mapped reads >=30%; every sample >=20%; paired BH q<0.05; LOSO direction >=80%",
    GateResult = if (pass_coverage && any(cand)) "pass_candidate_functions_present" else if (!pass_coverage) "failed_mapping_coverage" else "no_function_passed_inference_gate",
    ManuscriptRole = if (pass_coverage && any(cand)) "supplementary_exploratory" else "coverage_table_only",
    Reason = paste0("median mapped fraction=", sprintf("%.3f", median(coverage$MappedFraction)), "; minimum=", sprintf("%.3f", min(coverage$MappedFraction)), "; stable BH-significant functions=", sum(cand)),
    OriginalScriptReference = "function_data_process.R ordering/visual reference; prediction replaced by protected FAPROTAX 1.2.7",
    stringsAsFactors = FALSE
  )
  not18 <- data.frame(
    Marker = "18S", RequestedModule = "functional_prediction", RunStatus = "not_run_incompatible_marker_or_database",
    ManuscriptRole = "not_applicable",
    Reason = "The original 18S Greengenes/PICRUSt/BugBase lane is bacterial and invalid for mixed marine PR2 18S; PR2 v5.1.1 local trait columns contain no usable annotations.",
    stringsAsFactors = FALSE
  )
  write_tsv_atomic(not18, file.path(stage, "results/faprotax_18s_not_applicable_v2.tsv"))
  list(status = bind_rows(status, data.frame(
    ModuleID = "FUNCTION_18S", RunStatus = "not_run_incompatible_marker_or_database", NSamples = 75L, NIndependentUnits = 21L,
    Gate = "marker-compatible validated functional backend required", GateResult = "not_applicable",
    ManuscriptRole = "not_applicable", Reason = not18$Reason,
    OriginalScriptReference = "Original 18S Greengenes/PICRUSt/BugBase block rejected", stringsAsFactors = FALSE
  )), coverage = coverage, contrasts = contrasts)
}

faprotax_contrasts <- function(fun, design) {
  mat <- as.matrix(fun[, design$SampleID, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- fun$Function
  # Normalize overlapping function counts only for within-FAPROTAX comparisons.
  mat <- sweep(mat + 0.5, 2, colSums(mat + 0.5), "/")
  logp <- log(mat)
  specs <- surface_specs(design)[1:3]
  rows <- list()
  for (i in seq_along(specs)) {
    sp <- specs[[i]]
    md <- design[match(sp$samples, design$SampleID), , drop = FALSE]
    md <- paired_order(md, sp)
    pair_ids <- unique(as.character(md[[sp$pair]]))
    diffs <- sapply(pair_ids, function(k) {
      q <- md[as.character(md[[sp$pair]]) == k, , drop = FALSE]
      a <- q$SampleID[as.character(q[[sp$condition]]) == sp$target]
      b <- q$SampleID[as.character(q[[sp$condition]]) == sp$baseline]
      logp[, a] - logp[, b]
    })
    if (is.null(dim(diffs))) diffs <- matrix(diffs, nrow = nrow(logp))
    p <- apply(diffs, 1, function(v) {
      v <- v[is.finite(v)]
      if (length(unique(v)) < 2L) return(1)
      suppressWarnings(stats::wilcox.test(v, mu = 0, exact = FALSE)$p.value)
    })
    effect <- apply(diffs, 1, median, na.rm = TRUE)
    sites <- sort(unique(as.character(md$Site)))
    loo <- sapply(sites, function(s) {
      keep_pairs <- pair_ids[!pair_ids %in% as.character(md[[sp$pair]][as.character(md$Site) == s])]
      if (!length(keep_pairs)) return(rep(NA_real_, nrow(logp)))
      apply(diffs[, pair_ids %in% keep_pairs, drop = FALSE], 1, median, na.rm = TRUE)
    })
    stability <- rowMeans(sign(loo) == sign(effect), na.rm = TRUE)
    rows[[i]] <- data.frame(
      SurfaceID = sp$id, Function = rownames(logp), EffectLogProportion = effect,
      PValue = p, QValue = safe_bh(p), LOSODirectionStability = stability,
      Direction = ifelse(effect > 0, paste0("higher_in_", gsub(" ", "_", sp$target)), paste0("higher_in_", gsub(" ", "_", sp$baseline))),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

write_provenance <- function(stage) {
  input_paths <- unlist(paths, use.names = TRUE)
  inp <- data.frame(
    InputID = names(input_paths), Path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
    Bytes = as.numeric(file.info(input_paths)$size), SHA256 = vapply(input_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write_tsv_atomic(inp, file.path(stage, "provenance/input_sha256_v2.tsv"))
  contract <- data.frame(
    Component = c("input_counts", "display_depths", "DA_primary", "DA_sensitivities", "DA_parallelism", "DA_pairing", "function_16S", "function_18S", "visual_grammar"),
    Value = c(
      "Unrarefied 75-sample analysis-tier counts", "16S 20,000 and 18S 50,000 only for non-DA displays",
      "ANCOM-BC2 2.12.1; prevalence>=10%; BH q<0.05; pseudo-sensitivity required",
      "DESeq2 1.50.2 and paired ALDEx2 1.42.0 where estimable",
      paste0("ANCOM-BC2 n_cl=", ANCOMBC_WORKERS, " on a 24-logical-core Windows host; inferential settings unchanged"),
      "PairID_Substrate/LocationID fixed or random blocking; leave-one-Site-out direction stability",
      "Protected local FAPROTAX 1.2.7 via in-memory NumPy compatibility wrapper",
      "Not applicable: no valid mixed-marine-18S local functional backend",
      "Original EasyAmplicon ANCOM-BC2 volcano/dot and function bubble grammar; canonical palettes"
    ), stringsAsFactors = FALSE
  )
  write_tsv_atomic(contract, file.path(stage, "provenance/analysis_contract_v2.tsv"))
  dir.create(file.path(stage, "provenance"), recursive = TRUE, showWarnings = FALSE)
  capture.output(sessionInfo(), file = file.path(stage, "provenance/r_session_info_v2.txt"))
}

write_readme <- function(stage, daa, fap) {
  txt <- c(
    "# Original-style DA and functional analysis v2",
    "",
    "This package uses the curated 75-sample 16S and PR2-primary 18S inputs. The original EasyAmplicon scripts supplied analysis order and figure grammar; no legacy 18S numerical result was reused.",
    "",
    "## Differential abundance",
    "",
    "ANCOM-BC2 is primary. DESeq2 and paired ALDEx2 are sensitivity methods where the contrast is estimable. Candidate taxa require BH control, ANCOM-BC2 pseudo-count sensitivity, agreement with at least one sensitivity method, and >=80% leave-one-Site-out direction stability. Full method tables are retained even when no taxon passes.",
    "",
    "## Functional inference",
    "",
    paste0("16S FAPROTAX median mapped fraction: ", sprintf("%.3f", median(fap$coverage$MappedFraction)), "; minimum: ", sprintf("%.3f", min(fap$coverage$MappedFraction)), "."),
    "FAPROTAX outputs are taxonomy-inferred ecological functional groups, not measured genes or pathways. Functions overlap, so totals are not metagenome abundances.",
    "18S functional prediction is explicitly not applicable because the legacy Greengenes/PICRUSt/BugBase workflow is bacterial and no validated local mixed-marine-18S trait backend is available.",
    "",
    "## Original-script references",
    "",
    "- ancombc2_multigroup_script.R",
    "- ANCOMBC2_plots_and_tables_windows_v6.R",
    "- multigroup_compare_volcano_up11_taxonomy_fix.R",
    "- function_data_process.R (order/visual reference only; its hard-coded example was not executed)",
    "- FAPROTAX_1.2.7/collapse_table.py and FAPROTAX.txt (protected source)",
    "",
    paste0("Script version: ", SCRIPT_VERSION)
  )
  path <- file.path(stage, "README_daa_function_v2.md")
  writeLines(txt, path, useBytes = TRUE)
}

validate_output <- function(root) {
  miss <- required_outputs[!file.exists(file.path(root, required_outputs))]
  assert_true(!length(miss), paste0("Missing outputs: ", paste(miss, collapse = ", ")))
  st <- read_tsv_strict(file.path(root, "qc/module_status_daa_function_v2.tsv"))
  assert_true(all(c("ModuleID", "RunStatus", "GateResult", "ManuscriptRole", "Reason") %in% names(st)), "Module status schema invalid")
  assert_true(any(st$ModuleID == "FAPROTAX_16S") && any(st$ModuleID == "FUNCTION_18S"), "Functional status rows missing")
  sm <- read_tsv_strict(file.path(root, "results/daa_surface_summary_v2.tsv"))
  assert_true(nrow(sm) == 9L, "Expected 9 DA dataset/surface rows")
  assert_true(all(tolower(as.character(sm$ANCOMBC2_OK)) == "true"), "ANCOM-BC2 did not complete for every DA surface")
  assert_true(all(tolower(as.character(sm$DESeq2_OK)) == "true"), "DESeq2 did not complete for every DA surface")
  interaction <- sm$SurfaceID == "A60_Interaction"
  assert_true(all(tolower(as.character(sm$ALDEx2_OK[!interaction])) == "true"), "Paired ALDEx2 did not complete for every non-interaction surface")
  assert_true(all(tolower(as.character(sm$ALDEx2_OK[interaction])) == "false"), "Interaction surfaces must remain explicitly non-estimable by paired ALDEx2")
  assert_true(!any(grepl("round.*factor", as.character(sm$Notes), ignore.case = TRUE)), "ALDEx2 factor-input regression detected")
  mt <- read_tsv_strict(file.path(root, "results/daa_all_methods_v2.tsv"))
  assert_true(any(mt$Method == "ALDEx2"), "ALDEx2 result rows are missing")
  invisible(TRUE)
}

inp <- load_inputs()
if (mode == "dry-run") {
  dsets <- prepare_datasets(inp)
  specs <- surface_specs(inp$design)
  cat("DRY_RUN_OK\n")
  cat("datasets=", paste(vapply(dsets, `[[`, character(1), "dataset"), collapse = ","), "\n", sep = "")
  cat("surfaces=", paste(vapply(specs, `[[`, character(1), "id"), collapse = ","), "\n", sep = "")
  cat("samples=75; raw16_total=", sum(inp$c16), "; raw18_total=", sum(inp$c18), "\n", sep = "")
  cat("output_absent=", !dir.exists(output_root), "\n", sep = "")
  quit(status = 0)
}

if (mode == "check") {
  assert_true(dir.exists(output_root), paste0("Output package missing: ", output_root))
  validate_output(output_root)
  cat("CHECK_OK\n")
  quit(status = 0)
}

assert_true(!dir.exists(output_root), paste0("Refusing output overwrite: ", output_root))
stage <- paste0(output_root, ".staging-", Sys.getpid())
assert_true(!dir.exists(stage), paste0("Staging path exists: ", stage))
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
ok <- FALSE
on.exit(if (!ok && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

daa <- run_daa(inp, stage)
fap <- run_faprotax(inp, stage)
status <- bind_rows(daa$status, fap$status)
write_tsv_atomic(status, file.path(stage, "qc/module_status_daa_function_v2.tsv"))
write_provenance(stage)
write_readme(stage, daa, fap)
validate_output(stage)
assert_true(file.rename(stage, output_root), "Atomic package promotion failed")
ok <- TRUE
cat("RUN_OK\n")
