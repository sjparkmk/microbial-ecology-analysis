#!/usr/bin/env Rscript

# Advanced ecology analyses aligned to the original EasyAmplicon workflow.
# Numerical inputs are restricted to the curated 75-sample revision package.
# Original scripts/results are visual/sequence references only and are never edited.

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript so the project root can be resolved", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
script_dir <- dirname(script_file)
project_root <- normalizePath("/your/data/", winslash = "/", mustWork = TRUE) # insert your path
source(file.path(script_dir, "common_original_style.R"), local = TRUE)

suppressPackageStartupMessages({
  library(ranger)
  library(pROC)
  library(NST)
  library(igraph)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811)

cli <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% cli

paths <- list(
  design = file.path(project_root, "revision/re_analysis/00_design_v4/master_sample_metadata_v4.tsv"),
  count16 = file.path(project_root, "revision/re_analysis/04_input_16s_rarefied_20000_v1/results/16s_asv_counts_rarefied_20000_analysis_v1.tsv"),
  tax16 = file.path(project_root, "revision/re_analysis/01_input_16s_75_v1/16s_taxonomy_analysis_75_v1.tsv"),
  count18 = file.path(project_root, "revision/re_analysis/03_18s_rarefied_50000_v1/results/18s_asv_counts_rarefied_50000_analysis_v1.tsv"),
  tax18 = file.path(project_root, "revision/re_analysis/02_18s_unoise_pr2_v1/results/18s_taxonomy_pr2_primary_analysis_v1.tsv"),
  out_parent = file.path(project_root, "revision/re_analysis/12_original_pipeline_reanalysis_v1"),
  out = file.path(project_root, "revision/re_analysis/12_original_pipeline_reanalysis_v1/advanced_ecology_v2")
)

required_packages <- c("ranger", "pROC", "NST", "igraph", "ggrepel", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
assert_true(!length(missing_packages), paste0("Missing R packages: ", paste(missing_packages, collapse = ", ")))

design <- load_design(paths$design)
load_marker <- function(marker, count_path, tax_path, rank, depth) {
  mat <- load_count_matrix(count_path)$matrix
  mat <- align_counts_design(mat, design)$matrix
  assert_true(all(colSums(mat) == depth), paste0(marker, " rarefaction depth contract failed"))
  tax <- read_tsv_strict(tax_path)
  agg <- aggregate_counts_by_taxon(mat, tax, rank)
  list(marker = marker, rank = rank, depth = depth, asv = mat, tax = tax, agg = agg)
}

markers <- list(
  `16S` = load_marker("16S", paths$count16, paths$tax16, "Genus", 20000),
  `18S` = load_marker("18S", paths$count18, paths$tax18, "Class", 50000)
)

cat("Advanced ecology input contract\n")
cat("  design:", nrow(design), "samples,", nlevels(droplevels(design$Site)), "sites\n")
for (z in markers) {
  cat(" ", z$marker, ":", nrow(z$asv), "ASVs;", nrow(z$agg), z$rank, "taxa; depth", z$depth, "\n")
}
cat("  output:", paths$out, "\n")
if (dry_run) {
  cat("DRY RUN PASS: inputs, sample alignment, depths, taxonomy linkage, and packages validated.\n")
  quit(save = "no", status = 0)
}

assert_true(!dir.exists(paths$out), paste0("Refusing overwrite of existing output directory: ", paths$out))
dir.create(paths$out_parent, recursive = TRUE, showWarnings = FALSE)
stage <- file.path(paths$out_parent, paste0("._advanced_ecology_v2_staging_", Sys.getpid()))
assert_true(!file.exists(stage), paste0("Staging path already exists: ", stage))
dir.create(file.path(stage, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "logs"), recursive = TRUE, showWarnings = FALSE)
promoted <- FALSE
on.exit({
  if (!promoted && dir.exists(stage)) unlink(stage, recursive = TRUE, force = TRUE)
}, add = TRUE)

as_v2_table <- function(name) sub("_v1\\.tsv$", "_v2.tsv", name)
as_v2_figure <- function(name) sub("_v1$", "_v2", name)
as_v2_log <- function(name) sub("_v1\\.txt$", "_v2.txt", name)
table_path <- function(name) file.path(stage, "tables", as_v2_table(name))
figure_path <- function(name) file.path(stage, "figures", as_v2_figure(name))
log_path <- function(name) file.path(stage, "logs", as_v2_log(name))

status_rows <- list()
add_status <- function(module, marker, surface, status, evidence, limitation, manuscript_use) {
  status_rows[[length(status_rows) + 1L]] <<- data.frame(
    Module = module, Marker = marker, Surface = surface, Status = status,
    Evidence = evidence, Limitation = limitation, ManuscriptUse = manuscript_use,
    stringsAsFactors = FALSE
  )
}

# ------------------------------ Random forest ------------------------------

rf_surface <- function(marker_obj, task, nperm = 99L, seed = 1L) {
  if (task == "SampleType") {
    use <- !is.na(design$PairID_Substrate) & design$PairID_Substrate != ""
    pair_col <- "PairID_Substrate"
    lev <- c("Sea-sand", "Water")
  } else if (task == "Season") {
    use <- !is.na(design$PairID_Season) & design$PairID_Season != ""
    pair_col <- "PairID_Season"
    lev <- c("Spring", "Autumn")
  } else stop("Unknown RF task", call. = FALSE)
  d <- droplevels(design[use, , drop = FALSE])
  y <- factor(as.character(d[[task]]), levels = lev)
  pair <- as.character(d[[pair_col]])
  pair_check <- split(as.character(y), pair)
  assert_true(all(lengths(pair_check) == 2L), paste0(task, " RF requires two samples per pair"))
  assert_true(all(vapply(pair_check, function(v) setequal(v, lev), logical(1))),
              paste0(task, " RF pair labels are not balanced"))
  counts <- marker_obj$agg[, d$SampleID, drop = FALSE]
  counts <- counts[rownames(counts) != "Unassigned", , drop = FALSE]
  sites <- droplevels(d$Site)
  folds <- levels(sites)

  fit_loso <- function(y_fit, trees, importance = FALSE, seed_offset = 0L) {
    pred_rows <- list()
    imp_rows <- list()
    for (ii in seq_along(folds)) {
      held <- folds[[ii]]
      tr <- sites != held
      te <- !tr
      rel_tr <- sweep(counts[, tr, drop = FALSE], 2, colSums(counts[, tr, drop = FALSE]), "/")
      keep <- rowMeans(counts[, tr, drop = FALSE] > 0) >= 0.10 & rowMeans(rel_tr) >= 0.0001
      train_counts <- counts[keep, tr, drop = FALSE]
      test_counts <- counts[keep, te, drop = FALSE]
      assert_true(nrow(train_counts) >= 5L, paste0(marker_obj$marker, " ", task, " RF has <5 training taxa"))
      xtr <- clr_transform(train_counts, pseudocount = 0.5)
      xte <- clr_transform(test_counts, pseudocount = 0.5)
      colnames(xtr) <- make.unique(make.names(colnames(xtr)))
      colnames(xte) <- colnames(xtr)
      set.seed(seed + seed_offset + ii)
      fit <- ranger::ranger(
        x = as.data.frame(xtr), y = droplevels(y_fit[tr]), probability = TRUE,
        num.trees = trees, mtry = max(1L, floor(sqrt(ncol(xtr)))), min.node.size = 2L,
        importance = if (importance) "permutation" else "none", seed = seed + seed_offset + ii,
        num.threads = 1L
      )
      pp <- predict(fit, data = as.data.frame(xte))$predictions
      if (is.null(dim(pp))) {
        pp <- matrix(pp, ncol = length(lev), dimnames = list(NULL, lev))
      }
      miss <- setdiff(lev, colnames(pp))
      assert_true(!length(miss), paste0("RF probability columns missing: ", paste(miss, collapse = ", ")))
      pred <- lev[max.col(pp[, lev, drop = FALSE], ties.method = "first")]
      pred_rows[[ii]] <- data.frame(
        SampleID = d$SampleID[te], Site = as.character(d$Site[te]), PairID = pair[te],
        Observed = as.character(y_fit[te]), Predicted = pred,
        ProbabilityPositive = pp[, lev[[2]]], Fold = held, stringsAsFactors = FALSE
      )
      if (importance) {
        raw_names <- rownames(train_counts)
        imp <- fit$variable.importance
        # ranger uses syntactic feature names; retain taxa through positional matching.
        imp_rows[[ii]] <- data.frame(
          Taxon = raw_names[match(names(imp), colnames(xtr))], Fold = held,
          Importance = as.numeric(imp), stringsAsFactors = FALSE
        )
      }
    }
    list(pred = bind_rows(pred_rows), importance = bind_rows(imp_rows))
  }

  metric <- function(pred) {
    yy <- factor(pred$Observed, levels = lev)
    ph <- factor(pred$Predicted, levels = lev)
    sens <- vapply(lev, function(cl) mean(ph[yy == cl] == cl), numeric(1))
    auc <- suppressWarnings(as.numeric(pROC::auc(pROC::roc(
      response = as.integer(yy == lev[[2]]), predictor = pred$ProbabilityPositive,
      levels = c(0, 1), direction = "<", quiet = TRUE
    ))))
    c(BalancedAccuracy = mean(sens), AUC = auc,
      SensitivityLevel1 = sens[[1]], SensitivityLevel2 = sens[[2]])
  }

  obs_fit <- fit_loso(y, trees = 500L, importance = TRUE)
  obs_metric <- metric(obs_fit$pred)
  set.seed(seed + 100000L)
  null_rows <- vector("list", nperm)
  for (bb in seq_len(nperm)) {
    yp <- y
    swap <- stats::setNames(stats::runif(length(unique(pair))) < 0.5, unique(pair))
    for (p in names(swap)[swap]) {
      idx <- which(pair == p)
      yp[idx] <- rev(yp[idx])
    }
    perm_fit <- fit_loso(yp, trees = 150L, importance = FALSE, seed_offset = bb * 1000L)
    mm <- metric(perm_fit$pred)
    null_rows[[bb]] <- data.frame(Permutation = bb, t(mm), stringsAsFactors = FALSE)
  }
  null <- bind_rows(null_rows)
  p_perm <- (1 + sum(null$BalancedAccuracy >= obs_metric[["BalancedAccuracy"]])) / (nperm + 1)
  pred <- obs_fit$pred %>% mutate(
    Marker = marker_obj$marker, Rank = marker_obj$rank, Task = task,
    Correct = Observed == Predicted
  )
  imp <- obs_fit$importance %>%
    filter(!is.na(Taxon)) %>%
    group_by(Taxon) %>%
    summarise(FoldsSelected = n(), MeanImportance = mean(Importance),
              MedianImportance = median(Importance), PositiveFoldFraction = mean(Importance > 0),
              .groups = "drop") %>%
    arrange(desc(MeanImportance), Taxon) %>%
    mutate(RankImportance = row_number(), Marker = marker_obj$marker, Rank = marker_obj$rank, Task = task)
  metrics <- data.frame(
    Marker = marker_obj$marker, Rank = marker_obj$rank, Task = task,
    N = nrow(d), Sites = length(folds), Pairs = length(unique(pair)),
    BalancedAccuracy = obs_metric[["BalancedAccuracy"]], AUC = obs_metric[["AUC"]],
    SensitivityLevel1 = obs_metric[["SensitivityLevel1"]],
    SensitivityLevel2 = obs_metric[["SensitivityLevel2"]],
    Level1 = lev[[1]], Level2 = lev[[2]], Permutations = nperm,
    PermutationP_BalancedAccuracy = p_perm, stringsAsFactors = FALSE
  )
  null$Marker <- marker_obj$marker
  null$Task <- task
  informative <- is.finite(metrics$AUC) && metrics$BalancedAccuracy >= 0.65 &&
    metrics$AUC >= 0.65 && metrics$PermutationP_BalancedAccuracy <= 0.05
  status <- if (informative) "supplement_candidate" else "not_informative"
  add_status(
    "RF", marker_obj$marker, task, status,
    sprintf("LOSO-site balanced accuracy %.3f; AUC %.3f; paired-label permutation p=%.3f",
            metrics$BalancedAccuracy, metrics$AUC, metrics$PermutationP_BalancedAccuracy),
    "Taxa were filtered inside each training fold; labels were swapped only within matched field pairs. Site count remains small.",
    if (informative) "Supplementary predictive evidence only; do not interpret importance as causation."
    else "Do not cite as discriminatory evidence."
  )
  list(metrics = metrics, pred = pred, null = null, importance = imp)
}

cat("Running site-blocked RF surfaces...\n")
rf_results <- list()
seed_counter <- 0L
for (mn in names(markers)) for (task in c("SampleType", "Season")) {
  seed_counter <- seed_counter + 1L
  cat("  RF", mn, task, "\n")
  rf_results[[paste(mn, task, sep = "_")]] <- rf_surface(markers[[mn]], task, nperm = 99L,
                                                         seed = 12000L + seed_counter * 100000L)
}
rf_metrics <- bind_rows(lapply(rf_results, `[[`, "metrics"))
rf_predictions <- bind_rows(lapply(rf_results, `[[`, "pred"))
rf_null <- bind_rows(lapply(rf_results, `[[`, "null"))
rf_importance <- bind_rows(lapply(rf_results, `[[`, "importance"))
write_tsv_atomic(rf_metrics, table_path("RF_site_blocked_metrics_v1.tsv"))
write_tsv_atomic(rf_predictions, table_path("RF_site_blocked_predictions_v1.tsv"))
write_tsv_atomic(rf_null, table_path("RF_paired_label_permutation_null_v1.tsv"))
write_tsv_atomic(rf_importance, table_path("RF_fold_stable_feature_importance_v1.tsv"))

rf_plot_data <- rf_null %>%
  select(Marker, Task, Permutation, BalancedAccuracy) %>%
  mutate(Type = "Paired-label null")
rf_obs_plot <- rf_metrics %>% transmute(Marker, Task, BalancedAccuracy, Type = "Observed LOSO")
p_rf <- ggplot(rf_plot_data, aes(x = Type, y = BalancedAccuracy)) +
  geom_violin(fill = "grey88", color = "grey35", linewidth = 0.35, width = 0.8) +
  geom_jitter(width = 0.10, height = 0, alpha = 0.18, size = 0.65, color = "grey35") +
  geom_point(data = rf_obs_plot, aes(x = Type, y = BalancedAccuracy),
             inherit.aes = FALSE, shape = 21, size = 2.8, stroke = 0.55, fill = "#D55E00") +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35) +
  facet_grid(Task ~ Marker) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Balanced accuracy", title = "Site-blocked random-forest performance",
       subtitle = "Observed leave-one-site-out performance against 99 within-pair label permutations",
       caption = "Prediction is evaluated on held-out sites; 0.5 denotes chance for the balanced binary surfaces.") +
  theme_project_original_style(8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot_bundle(p_rf, figure_path("Fig_Adv_RF_site_blocked_performance_v1"), 178, 128)

rf_qualifying <- rf_metrics %>%
  filter(BalancedAccuracy >= 0.65, AUC >= 0.65, PermutationP_BalancedAccuracy <= 0.05) %>%
  select(Marker, Task)
imp_plot <- rf_importance %>% inner_join(rf_qualifying, by = c("Marker", "Task")) %>%
  group_by(Marker, Task) %>% slice_max(MeanImportance, n = 10, with_ties = FALSE) %>%
  ungroup() %>% mutate(Facet = paste(Marker, Task, sep = " | "),
                       Label = reorder(paste(Facet, Taxon, sep = "___"), MeanImportance))
p_imp <- ggplot(imp_plot, aes(x = MeanImportance, y = Label, fill = Marker)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  facet_wrap(~Facet, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_manual(values = c(`16S` = "#0072B2", `18S` = "#009E73")) +
  labs(x = "Mean permutation importance across selected LOSO folds", y = NULL,
       title = "Fold-stable random-forest feature importance",
       caption = "Importance is predictive, not causal; feature filtering was repeated within every training fold.") +
  theme_project_original_style(7.5) + theme(legend.position = "none")
save_plot_bundle(p_imp, figure_path("Fig_Adv_RF_feature_importance_v1"), 178, 190)

# --------------------------- Sloan neutral model ---------------------------

ncm_fit <- function(counts) {
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  n_samples <- ncol(counts)
  N <- mean(colSums(counts))
  mean_abund <- rowSums(counts) / sum(counts)
  occurrence <- rowMeans(counts > 0)
  pclip <- pmin(pmax(mean_abund, 1e-12), 1 - 1e-12)
  predict_occ <- function(m) {
    stats::pbeta(1 / N, N * m * pclip, N * m * (1 - pclip), lower.tail = FALSE)
  }
  objective <- function(log_m) sum((occurrence - predict_occ(exp(log_m)))^2)
  opt <- stats::optimize(objective, interval = log(c(1e-8, 100)))
  m <- exp(opt$minimum)
  pred <- pmin(pmax(predict_occ(m), 0), 1)
  sst <- sum((occurrence - mean(occurrence))^2)
  r2 <- if (sst > 0) 1 - sum((occurrence - pred)^2) / sst else NA_real_
  lower <- stats::qbinom(0.025, size = n_samples, prob = pred) / n_samples
  upper <- stats::qbinom(0.975, size = n_samples, prob = pred) / n_samples
  classification <- ifelse(occurrence > upper, "Above prediction",
                           ifelse(occurrence < lower, "Below prediction", "Within prediction"))
  data <- data.frame(
    Taxon = rownames(counts), MeanRelativeAbundance = mean_abund,
    OccurrenceFrequency = occurrence, PredictedOccurrence = pred,
    Lower95 = lower, Upper95 = upper, Classification = classification,
    stringsAsFactors = FALSE
  )
  list(m = m, R2 = r2, N = N, n = n_samples, data = data)
}

ncm_surface <- function(marker_obj, sample_type, nboot = 199L, seed = 1L) {
  use <- !is.na(design$PairID_Substrate) & design$PairID_Substrate != "" &
    as.character(design$SampleType) == sample_type
  d <- droplevels(design[use, , drop = FALSE])
  x <- marker_obj$asv[, d$SampleID, drop = FALSE]
  full <- ncm_fit(x)
  sites <- unique(as.character(d$Site))
  loso <- lapply(sites, function(s) {
    fit <- ncm_fit(x[, as.character(d$Site) != s, drop = FALSE])
    data.frame(LeftOutSite = s, m = fit$m, R2 = fit$R2, N = fit$N, Samples = fit$n)
  }) %>% bind_rows()
  set.seed(seed)
  boot <- vector("list", nboot)
  for (bb in seq_len(nboot)) {
    draw <- sample(sites, length(sites), replace = TRUE)
    idx <- unlist(lapply(draw, function(s) which(as.character(d$Site) == s)), use.names = FALSE)
    fit <- ncm_fit(x[, idx, drop = FALSE])
    boot[[bb]] <- data.frame(Bootstrap = bb, m = fit$m, R2 = fit$R2)
  }
  boot <- bind_rows(boot)
  summary <- data.frame(
    Marker = marker_obj$marker, SampleType = sample_type, Rank = "ASV",
    Samples = nrow(d), Sites = length(sites), Features = nrow(full$data),
    N = full$N, m = full$m, R2 = full$R2,
    Bootstrap_m_Lower95 = quantile(boot$m, 0.025, na.rm = TRUE),
    Bootstrap_m_Median = median(boot$m, na.rm = TRUE),
    Bootstrap_m_Upper95 = quantile(boot$m, 0.975, na.rm = TRUE),
    Bootstrap_R2_Lower95 = quantile(boot$R2, 0.025, na.rm = TRUE),
    Bootstrap_R2_Median = median(boot$R2, na.rm = TRUE),
    Bootstrap_R2_Upper95 = quantile(boot$R2, 0.975, na.rm = TRUE),
    BootstrapPositiveR2Fraction = mean(boot$R2 > 0, na.rm = TRUE),
    LOSOPositiveR2Fraction = mean(loso$R2 > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  informative <- is.finite(full$R2) && full$R2 >= 0.10 &&
    summary$Bootstrap_R2_Lower95 > 0 && summary$BootstrapPositiveR2Fraction >= 0.95 &&
    summary$LOSOPositiveR2Fraction >= 0.80
  status <- if (informative) "supplement_candidate" else "not_informative"
  add_status(
    "NCM", marker_obj$marker, paste0("A72_", sample_type), status,
    sprintf("Sloan ASV model R2=%.3f; m=%.4g; positive-R2 bootstrap=%.2f; LOSO=%.2f",
            full$R2, full$m, summary$BootstrapPositiveR2Fraction, summary$LOSOPositiveR2Fraction),
    "Site-cluster bootstrap and leave-one-site-out assess stability, but the model remains phenomenological.",
    if (informative) "Supplementary community-assembly description; avoid deterministic/stochastic causal wording."
    else "Do not use for a neutral-assembly claim."
  )
  full$data$Marker <- marker_obj$marker
  full$data$SampleType <- sample_type
  loso$Marker <- marker_obj$marker
  loso$SampleType <- sample_type
  boot$Marker <- marker_obj$marker
  boot$SampleType <- sample_type
  list(summary = summary, data = full$data, loso = loso, boot = boot)
}

cat("Running Sloan neutral-community models...\n")
ncm_results <- list()
ii <- 0L
for (mn in names(markers)) for (st in c("Sea-sand", "Water")) {
  ii <- ii + 1L
  cat("  NCM", mn, st, "\n")
  ncm_results[[paste(mn, st, sep = "_")]] <- ncm_surface(markers[[mn]], st, 199L, 73000L + ii)
}
ncm_summary <- bind_rows(lapply(ncm_results, `[[`, "summary"))
ncm_data <- bind_rows(lapply(ncm_results, `[[`, "data"))
ncm_loso <- bind_rows(lapply(ncm_results, `[[`, "loso"))
ncm_boot <- bind_rows(lapply(ncm_results, `[[`, "boot"))
write_tsv_atomic(ncm_summary, table_path("NCM_A72_summary_v1.tsv"))
write_tsv_atomic(ncm_data, table_path("NCM_A72_taxon_fit_source_v1.tsv"))
write_tsv_atomic(ncm_loso, table_path("NCM_A72_leave_one_site_out_v1.tsv"))
write_tsv_atomic(ncm_boot, table_path("NCM_A72_site_cluster_bootstrap_v1.tsv"))

ncm_curve <- ncm_data %>%
  group_by(Marker, SampleType) %>% arrange(MeanRelativeAbundance, .by_group = TRUE) %>%
  select(Marker, SampleType, MeanRelativeAbundance, PredictedOccurrence, Lower95, Upper95) %>%
  distinct() %>% ungroup()
p_ncm <- ggplot(ncm_data, aes(x = MeanRelativeAbundance, y = OccurrenceFrequency, color = Classification)) +
  geom_point(alpha = 0.28, size = 0.55) +
  geom_ribbon(data = ncm_curve,
              aes(x = MeanRelativeAbundance, y = PredictedOccurrence, ymin = Lower95, ymax = Upper95,
                  group = interaction(Marker, SampleType)),
              inherit.aes = FALSE, fill = "grey75", alpha = 0.45) +
  geom_line(data = ncm_curve,
            aes(x = MeanRelativeAbundance, y = PredictedOccurrence,
                group = interaction(Marker, SampleType)), inherit.aes = FALSE,
            linewidth = 0.55, color = "black") +
  facet_grid(SampleType ~ Marker, scales = "free_x") +
  scale_x_log10(labels = scales::label_scientific()) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c(`Above prediction` = "#D55E00", `Within prediction` = "#7F7F7F",
                                `Below prediction` = "#0072B2")) +
  labs(x = "Mean relative abundance (log10)", y = "Occurrence frequency",
       title = "Sloan neutral-community model on matched A72 samples",
       subtitle = "ASV-level fits shown separately by marker and sample type",
       caption = "Stability was evaluated by site-cluster bootstrap and leave-one-site-out refitting.",
       color = "Observed position") +
  theme_project_original_style(8) + theme(legend.position = "bottom")
save_plot_bundle(p_ncm, figure_path("Fig_Adv_NCM_A72_v1"), 178, 145)

# -------------------- Taxonomic normalized stochasticity -------------------

nst_surface <- function(marker_obj, task, seed = 1L, rand_null = 199L, rand_boot = 999L) {
  if (task == "SampleType") {
    use <- !is.na(design$PairID_Substrate) & design$PairID_Substrate != ""
    group_var <- "SampleType"
    surface <- "A72_SampleType"
  } else {
    use <- !is.na(design$PairID_Season) & design$PairID_Season != ""
    group_var <- "Season"
    surface <- "A66_Season"
  }
  d <- droplevels(design[use, , drop = FALSE])
  x <- marker_obj$agg[, d$SampleID, drop = FALSE]
  rel <- sweep(x, 2, colSums(x), "/")
  keep <- rownames(x) != "Unassigned" & rowMeans(x > 0) >= 0.20 & rowMeans(rel) >= 0.00001
  x <- x[keep, , drop = FALSE]
  assert_true(nrow(x) >= 20L, paste0(marker_obj$marker, " ", task, " NST has <20 retained taxa"))
  group <- data.frame(Group = droplevels(d[[group_var]]), row.names = d$SampleID)
  meta_group <- data.frame(Site = droplevels(d$Site), row.names = d$SampleID)

  run_one <- function(method, abundance, seed_offset) {
    set.seed(seed + seed_offset)
    progress <- capture.output({
      fit <- NST::tNST(
        comm = t(x), group = group, meta.group = meta_group,
        dist.method = method, abundance.weighted = abundance,
        rand = rand_null, output.rand = TRUE, nworker = 1L,
        null.model = "PF", between.group = FALSE
      )
    })
    set.seed(seed + seed_offset + 50000L)
    boot <- NST::nst.boot(fit, group = group, rand = rand_boot, trace = FALSE,
                          out.detail = TRUE, between.group = FALSE, nworker = 1L)
    idx <- fit$index.grp
    nst_col <- grep("^NST\\.i", names(idx), value = TRUE)[1]
    out <- data.frame(
      Marker = marker_obj$marker, Rank = marker_obj$rank, Surface = surface,
      Group = as.character(idx$group), PairwiseComparisons = idx$size,
      Method = if (abundance) paste0(method, "_abundance") else paste0(method, "_incidence"),
      NST = idx[[nst_col]], RetainedTaxa = nrow(x), Samples = nrow(d),
      Sites = nlevels(droplevels(d$Site)), NullRandomizations = rand_null,
      stringsAsFactors = FALSE
    )
    bs <- boot$summary %>% filter(Index == "NST") %>%
      transmute(Group = as.character(Group), BootstrapMean = mean, BootstrapSD = stdev,
                BootstrapMin = Min, BootstrapQ25 = Quantile25, BootstrapMedian = Median,
                BootstrapQ75 = Quantile75, BootstrapMax = Max)
    out <- left_join(out, bs, by = "Group")
    list(summary = out, compare = boot$compare, progress = progress)
  }

  primary <- run_one("bray", TRUE, 0L)
  sensitivity <- run_one("jaccard", FALSE, 10000L)
  out <- bind_rows(primary$summary, sensitivity$summary)
  directions <- out %>% group_by(Method) %>% arrange(Group, .by_group = TRUE) %>%
    summarise(Direction = sign(diff(NST)), Difference = diff(NST), .groups = "drop")
  concordant <- nrow(directions) == 2L && length(unique(directions$Direction[directions$Direction != 0])) == 1L
  material <- all(abs(directions$Difference) >= 0.05)
  status <- if (concordant && material) "descriptive_only" else "not_informative"
  add_status(
    "tNST", marker_obj$marker, surface, status,
    paste0("PF null model; Bray-abundance and Jaccard-incidence direction concordant=", concordant,
           "; both absolute group differences >=0.05=", material),
    "Package bootstrap resamples pairwise values rather than independent sites; site is used as metacommunity, so comparisons remain descriptive.",
    if (status == "descriptive_only") "Supplementary descriptive surface only; no p-value or causal assembly claim."
    else "Do not use for a group-level assembly claim."
  )
  list(summary = out, directions = directions,
       compare = bind_rows(transform(primary$compare, Method = "bray_abundance"),
                           transform(sensitivity$compare, Method = "jaccard_incidence")))
}

cat("Running tNST primary and sensitivity surfaces...\n")
nst_results <- list()
ii <- 0L
for (mn in names(markers)) for (task in c("SampleType", "Season")) {
  ii <- ii + 1L
  cat("  tNST", mn, task, "\n")
  nst_results[[paste(mn, task, sep = "_")]] <- nst_surface(markers[[mn]], task,
                                                            seed = 91000L + ii * 1000L)
}
nst_summary <- bind_rows(lapply(nst_results, `[[`, "summary"))
nst_directions <- bind_rows(Map(function(x, nm) transform(x$directions, Analysis = nm),
                                nst_results, names(nst_results)))
nst_compare <- bind_rows(Map(function(x, nm) transform(x$compare, Analysis = nm),
                             nst_results, names(nst_results)))
write_tsv_atomic(nst_summary, table_path("tNST_summary_and_bootstrap_v1.tsv"))
write_tsv_atomic(nst_directions, table_path("tNST_method_sensitivity_direction_v1.tsv"))
write_tsv_atomic(nst_compare, table_path("tNST_package_pairwise_bootstrap_comparison_descriptive_v1.tsv"))

p_nst <- ggplot(nst_summary, aes(x = Group, y = NST, color = Method, group = Method)) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35, color = "grey40") +
  geom_errorbar(aes(ymin = BootstrapQ25, ymax = BootstrapQ75), width = 0.12,
                position = position_dodge(width = 0.35), linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.35), size = 2.2) +
  facet_wrap(~Surface + Marker, scales = "free_x", ncol = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c(bray_abundance = "#D55E00", jaccard_incidence = "#0072B2")) +
  labs(x = NULL, y = "Taxonomic normalized stochasticity ratio (tNST)",
       title = "Taxonomic normalized stochasticity under paired analysis surfaces",
       subtitle = "Points are observed group estimates; bars are package bootstrap interquartile ranges",
       caption = "Site defines the metacommunity. Bootstrap intervals are descriptive because pairwise turnover values are not independent.",
       color = "Sensitivity surface") +
  theme_project_original_style(8) + theme(legend.position = "bottom")
save_plot_bundle(p_nst, figure_path("Fig_Adv_tNST_sensitivity_v1"), 178, 155)

# ------------------------- Residual association networks -------------------

network_features <- function(marker_obj, max_taxa = 50L) {
  x <- marker_obj$agg
  rel <- sweep(x, 2, colSums(x), "/")
  keep <- rownames(x) != "Unassigned" & rowMeans(x > 0) >= 0.20 & rowMeans(rel) >= 0.0005
  x <- x[keep, , drop = FALSE]
  if (nrow(x) > max_taxa) {
    ord <- order(rowMeans(rel[rownames(x), , drop = FALSE]), decreasing = TRUE)
    x <- x[ord[seq_len(max_taxa)], , drop = FALSE]
  }
  clr <- clr_transform(x, pseudocount = 0.5)
  design_df <- data.frame(Site = droplevels(design$Site), Season = droplevels(design$Season),
                          SampleType = droplevels(design$SampleType))
  mm <- model.matrix(~ Site + Season + SampleType, data = design_df)
  qr_mm <- qr(mm)
  resid <- apply(clr, 2, function(v) qr.resid(qr_mm, v))
  if (is.vector(resid)) resid <- matrix(resid, ncol = 1L)
  colnames(resid) <- colnames(clr)
  rownames(resid) <- rownames(clr)
  list(counts = x, residuals = resid,
       filter = data.frame(Marker = marker_obj$marker, Rank = marker_obj$rank,
                           Taxon = rownames(x), Prevalence = rowMeans(x > 0),
                           MeanRelativeAbundance = rowMeans(rel[rownames(x), , drop = FALSE])))
}

cor_p <- function(r, n) {
  rr <- pmin(pmax(r, -0.999999), 0.999999)
  tt <- rr * sqrt((n - 2) / pmax(1e-12, 1 - rr^2))
  2 * stats::pt(abs(tt), df = n - 2, lower.tail = FALSE)
}

cluster_boot_cor <- function(a, b = NULL, nboot = 199L, seed = 1L) {
  clusters <- as.character(design$LocationID)
  lev <- unique(clusters)
  set.seed(seed)
  if (is.null(b)) {
    nr <- ncol(a)
    vals <- array(NA_real_, dim = c(nr, nr, nboot))
  } else {
    vals <- array(NA_real_, dim = c(ncol(a), ncol(b), nboot))
  }
  for (bb in seq_len(nboot)) {
    draw <- sample(lev, length(lev), replace = TRUE)
    idx <- unlist(lapply(draw, function(z) which(clusters == z)), use.names = FALSE)
    if (is.null(b)) vals[, , bb] <- stats::cor(a[idx, , drop = FALSE])
    else vals[, , bb] <- stats::cor(a[idx, , drop = FALSE], b[idx, , drop = FALSE])
  }
  vals
}

edge_table <- function(a, b = NULL, marker_a, marker_b, nboot = 199L, seed = 1L) {
  if (is.null(b)) {
    rmat <- stats::cor(a)
    idx <- which(upper.tri(rmat), arr.ind = TRUE)
    edges <- data.frame(Node1 = colnames(a)[idx[, 1]], Node2 = colnames(a)[idx[, 2]],
                        R = rmat[idx], stringsAsFactors = FALSE)
    boot <- cluster_boot_cor(a, NULL, nboot, seed)
    boot_mat <- vapply(seq_len(nboot), function(k) boot[, , k][idx], numeric(nrow(idx)))
  } else {
    rmat <- stats::cor(a, b)
    idx <- which(matrix(TRUE, nrow(rmat), ncol(rmat)), arr.ind = TRUE)
    edges <- data.frame(Node1 = rownames(rmat)[idx[, 1]], Node2 = colnames(rmat)[idx[, 2]],
                        R = rmat[idx], stringsAsFactors = FALSE)
    boot <- cluster_boot_cor(a, b, nboot, seed)
    boot_mat <- vapply(seq_len(nboot), function(k) boot[, , k][idx], numeric(nrow(idx)))
  }
  edges$Marker1 <- marker_a
  edges$Marker2 <- marker_b
  edges$P <- cor_p(edges$R, nrow(a))
  edges$Q_BH <- safe_bh(edges$P)
  edges$BootstrapLower95 <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  edges$BootstrapMedian <- apply(boot_mat, 1, median, na.rm = TRUE)
  edges$BootstrapUpper95 <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  edges$SameSignFraction <- rowMeans(sign(boot_mat) == sign(edges$R), na.rm = TRUE)
  edges$ScreenPass <- edges$Q_BH < 0.05 & abs(edges$R) >= 0.40 & edges$SameSignFraction >= 0.80 &
    ((edges$BootstrapLower95 > 0 & edges$R > 0) | (edges$BootstrapUpper95 < 0 & edges$R < 0))
  edges$Stable <- edges$Q_BH < 0.01 & abs(edges$R) >= 0.70 & edges$SameSignFraction >= 0.95 &
    ((edges$BootstrapLower95 > 0 & edges$R > 0) | (edges$BootstrapUpper95 < 0 & edges$R < 0))
  edges
}

cat("Running residual association networks...\n")
net_features <- lapply(markers, network_features)
filter_table <- bind_rows(lapply(net_features, `[[`, "filter"))
write_tsv_atomic(filter_table, table_path("Network_retained_taxa_v1.tsv"))

edge16 <- edge_table(net_features$`16S`$residuals, marker_a = "16S", marker_b = "16S",
                     nboot = 199L, seed = 130001L)
edge18 <- edge_table(net_features$`18S`$residuals, marker_a = "18S", marker_b = "18S",
                     nboot = 199L, seed = 130002L)
edge_cross <- edge_table(net_features$`16S`$residuals, net_features$`18S`$residuals,
                         marker_a = "16S", marker_b = "18S", nboot = 199L, seed = 130003L)
edges <- bind_rows(transform(edge16, Network = "16S_within"),
                   transform(edge18, Network = "18S_within"),
                   transform(edge_cross, Network = "16S_18S_cross"))
write_tsv_atomic(edges, table_path("Network_all_tested_edges_v1.tsv"))
write_tsv_atomic(edges %>% filter(Stable), table_path("Network_stable_edges_v1.tsv"))

network_summary <- edges %>% group_by(Network, Marker1, Marker2) %>%
  summarise(TestedEdges = n(), BHSignificantEdges = sum(Q_BH < 0.05),
            ScreenPassEdges = sum(ScreenPass),
            StableEdges = sum(Stable), StablePositive = sum(Stable & R > 0),
            StableNegative = sum(Stable & R < 0),
            StableDensity = StableEdges / TestedEdges, .groups = "drop")
write_tsv_atomic(network_summary, table_path("Network_summary_v1.tsv"))

for (i in seq_len(nrow(network_summary))) {
  rr <- network_summary[i, ]
  informative <- rr$StableEdges >= 3L && rr$StableEdges <= 60L && rr$StableDensity <= 0.05
  status <- if (informative) "supplement_candidate" else "not_informative"
  add_status(
    "Network", paste(rr$Marker1, rr$Marker2, sep = "-"), rr$Network, status,
    sprintf("%d/%d core edges passed BH q<0.01, |r|>=0.70, same-sign>=0.95, and cluster-bootstrap CI exclusion of zero (screen-pass=%d)",
            rr$StableEdges, rr$TestedEdges, rr$ScreenPassEdges),
    "Associations are CLR-residual correlations after Site+Season+SampleType adjustment; they are not direct interactions.",
    if (informative) "Supplementary hypothesis-generating network only."
    else "Do not present as an ecological network result."
  )
}

plot_network <- function(edge_df, network_name, seed = 1L) {
  stable <- edge_df %>% filter(Network == network_name, Stable)
  if (!nrow(stable)) {
    diag <- edge_df %>% filter(Network == network_name)
    return(ggplot(diag, aes(x = R, y = -log10(pmax(Q_BH, 1e-300)))) +
      geom_point(alpha = 0.35, size = 0.75, color = "grey35") +
      geom_vline(xintercept = c(-0.4, 0.4), linetype = 2, linewidth = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
      labs(x = "Residual correlation", y = expression(-log[10](BH~q)),
           title = paste0(network_name, ": no stable edge"),
           subtitle = "Diagnostic surface retained to document the explicit negative result") +
      theme_project_original_style(8))
  }
  vertices <- unique(c(stable$Node1, stable$Node2))
  g <- graph_from_data_frame(stable[, c("Node1", "Node2")], directed = FALSE,
                             vertices = data.frame(name = vertices))
  set.seed(seed)
  xy <- layout_with_fr(g, niter = 1000)
  nodes <- data.frame(Node = V(g)$name, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
  nodes$Marker <- ifelse(nodes$Node %in% stable$Node1[stable$Marker1 == "16S"] |
                           nodes$Node %in% stable$Node2[stable$Marker2 == "16S"], "16S", "18S")
  # Disambiguate cross-domain labels while preserving readable labels.
  if (network_name == "16S_18S_cross") {
    stable$Key1 <- paste0("16S::", stable$Node1)
    stable$Key2 <- paste0("18S::", stable$Node2)
    g <- graph_from_data_frame(stable[, c("Key1", "Key2")], directed = FALSE)
    set.seed(seed)
    xy <- layout_as_bipartite(g, types = grepl("^18S::", V(g)$name))
    nodes <- data.frame(Key = V(g)$name, x = xy[, 1], y = xy[, 2], stringsAsFactors = FALSE)
    nodes$Marker <- ifelse(grepl("^16S::", nodes$Key), "16S", "18S")
    nodes$Node <- sub("^[^:]+::", "", nodes$Key)
    ed <- stable %>% left_join(nodes %>% select(Key, x, y), by = c("Key1" = "Key")) %>%
      rename(x1 = x, y1 = y) %>% left_join(nodes %>% select(Key, x, y), by = c("Key2" = "Key")) %>%
      rename(x2 = x, y2 = y)
  } else {
    ed <- stable %>% left_join(nodes, by = c("Node1" = "Node")) %>%
      rename(x1 = x, y1 = y) %>% left_join(nodes, by = c("Node2" = "Node"), suffix = c("", ".two")) %>%
      rename(x2 = x, y2 = y)
  }
  node_key <- if ("Key" %in% names(nodes)) nodes$Key else nodes$Node
  nodes$Degree <- as.numeric(igraph::degree(g)[node_key])
  label_n <- min(15L, nrow(nodes))
  label_idx <- order(nodes$Degree, decreasing = TRUE)[seq_len(label_n)]
  nodes$Label <- ""
  nodes$Label[label_idx] <- nodes$Node[label_idx]
  ggplot() +
    geom_segment(data = ed, aes(x = x1, y = y1, xend = x2, yend = y2,
                                color = R > 0, linewidth = abs(R)), alpha = 0.60) +
    geom_point(data = nodes, aes(x = x, y = y, fill = Marker), shape = 21, size = 2.6,
               color = "black", stroke = 0.35) +
    ggrepel::geom_text_repel(data = nodes %>% filter(Label != ""),
                             aes(x = x, y = y, label = Label), size = 2.2,
                             max.overlaps = Inf, min.segment.length = 0, seed = seed) +
    scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"),
                       labels = c(`TRUE` = "Positive", `FALSE` = "Negative"), name = "Association") +
    scale_fill_manual(values = c(`16S` = "#56B4E9", `18S` = "#009E73")) +
    scale_linewidth_continuous(range = c(0.3, 1.2), guide = "none") +
    coord_equal(clip = "off") +
    labs(x = NULL, y = NULL, title = network_name,
         subtitle = paste(nrow(stable), "stable residual associations")) +
    theme_project_original_style(8) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "bottom")
}

p_net16 <- plot_network(edges, "16S_within", 140001L)
p_net18 <- plot_network(edges, "18S_within", 140002L)
p_netcross <- plot_network(edges, "16S_18S_cross", 140003L)
p_network <- ((p_net16 | p_net18) / p_netcross) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Design-adjusted residual association networks",
    subtitle = "CLR residuals adjusted for Site + Season + SampleType.\nCore edges require BH, effect-size, and LocationID bootstrap support.",
    caption = "Edges denote robust residual associations, not direct biological interactions.",
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9), plot.caption = element_text(size = 7))
  ) & theme(legend.position = "bottom")
save_plot_bundle(p_network, figure_path("Fig_Adv_residual_association_networks_v1"), 178, 190)

# ------------------------- Gating and provenance ---------------------------

status <- bind_rows(status_rows)
status_levels <- c("main_candidate", "supplement_candidate", "descriptive_only",
                   "not_informative", "not_applicable", "failed_qc")
assert_true(all(status$Status %in% status_levels), "Invalid status label generated")
write_tsv_atomic(status, table_path("advanced_ecology_status_gate_v1.tsv"))

parameters <- data.frame(
  Module = c("RF", "RF", "NCM", "NCM", "tNST", "tNST", "Network", "Network", "Network"),
  Parameter = c("evaluation", "null_test", "model", "stability", "primary_and_sensitivity",
                "bootstrap_caveat", "transform_and_adjustment", "edge_gate", "bootstrap_unit"),
  Value = c(
    "leave-one-Site-out on matched A72 SampleType and A66 Season surfaces; filtering within training fold",
    "99 label permutations by swapping labels within matched pair; 500 observed and 150 null trees",
    "Sloan beta-distribution occurrence model at ASV level, separately for Sea-sand and Water in A72",
    "199 Site-cluster bootstraps and leave-one-Site-out refits",
    "Bray abundance-weighted plus Jaccard incidence; PF null; 199 randomizations; Site metacommunity",
    "999 package bootstraps summarize pairwise values and are descriptive, not independent-site inference",
    "Genus (16S) or PR2 Class (18S), prevalence>=20%, mean RA>=0.05%, top 50; CLR residualized by Site+Season+SampleType",
    "screen: BH q<0.05, |r|>=0.40, same-sign>=0.80; plotted core: BH q<0.01, |r|>=0.70, same-sign>=0.95, percentile CI excludes zero",
    "199 resamples of the 21 fixed LocationID clusters"
  ), stringsAsFactors = FALSE
)
write_tsv_atomic(parameters, table_path("advanced_ecology_parameters_v1.tsv"))

figure_contracts <- data.frame(
  Figure = c("Fig_Adv_RF_site_blocked_performance_v1", "Fig_Adv_RF_feature_importance_v1",
             "Fig_Adv_NCM_A72_v1", "Fig_Adv_tNST_sensitivity_v1",
             "Fig_Adv_residual_association_networks_v1"),
  DefendedSentence = c(
    "Site-blocked predictive performance is compared with a matched-pair label null.",
    "Predictive taxa are shown only when repeatedly selected across LOSO training folds.",
    "Neutral-model fit and deviations are displayed with site-level stability diagnostics in source tables.",
    "tNST direction is checked across abundance and incidence definitions without treating pairwise turnover as independent.",
    "Only design-adjusted associations stable to fixed-location cluster resampling are drawn."
  ),
  Backend = "R/ggplot2",
  EditableVector = "PDF+SVG",
  Raster = "PNG 600 dpi",
  ColorContract = c("Observed orange; null grey", "16S blue; 18S green",
                    "Above orange; within grey; below blue", "Bray orange; Jaccard blue",
                    "Positive orange; negative blue; 16S light blue; 18S green"),
  stringsAsFactors = FALSE
)
figure_contracts$Figure <- as_v2_figure(figure_contracts$Figure)
write_tsv_atomic(figure_contracts, table_path("advanced_ecology_figure_contracts_v1.tsv"))

input_manifest <- data.frame(
  Role = c("Design", "16S rarefied count", "16S taxonomy", "18S rarefied count", "18S PR2 taxonomy"),
  Path = unlist(paths[c("design", "count16", "tax16", "count18", "tax18")], use.names = FALSE),
  MD5 = unname(tools::md5sum(unlist(paths[c("design", "count16", "tax16", "count18", "tax18")], use.names = FALSE))),
  stringsAsFactors = FALSE
)
write_tsv_atomic(input_manifest, table_path("advanced_ecology_input_manifest_v1.tsv"))

supersession <- data.frame(
  Package = c("advanced_ecology_v1", "advanced_ecology_v2"),
  Decision = c("qa_rejected", "accepted"),
  ManuscriptUse = c("not_manuscript_use", "use_status_gate"),
  Reason = c(
    "Superseded after visual QA: network screen was over-dense and tNST facets retained empty cross-surface categories.",
    "Accepted after stricter core-network gate, corrected tNST facets, conservative NCM bootstrap gate, and figure QA."
  ),
  ReplacementPath = c("advanced_ecology_v2", "advanced_ecology_v2"),
  stringsAsFactors = FALSE
)
write_tsv_atomic(supersession, table_path("advanced_ecology_supersession_ledger_v1.tsv"))

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, log_path("R_sessionInfo_v1.txt"), useBytes = TRUE)

readme <- c(
  "# Advanced ecology reanalysis v2",
  "",
  "This package uses only the curated 75-sample revision inputs. The legacy EasyAmplicon scripts informed module order and visual grammar; their numerical assumptions were not reused.",
  "",
  "## Interpretation gate",
  "",
  "Use `tables/advanced_ecology_status_gate_v1.tsv` before citing any result. `supplement_candidate` denotes a stable, limited supplementary surface; `descriptive_only` is not confirmatory; `not_informative` must not be promoted to the manuscript as supporting evidence.",
  "",
  "RF is evaluated by leave-one-site-out prediction and matched-pair label permutations. NCM is checked by Site-cluster bootstrap and leave-one-Site-out refitting. tNST includes two dissimilarity definitions but its package bootstrap operates on pairwise values and is therefore descriptive. Networks are CLR residual associations adjusted for Site, Season, and SampleType and require BH plus LocationID-cluster-bootstrap support.",
  "",
  "All figures are exported as editable PDF/SVG and 600-dpi PNG with matching source TSVs."
)
writeLines(readme, file.path(stage, "README_advanced_ecology_v2.md"), useBytes = TRUE)

# Final integrity checks before the single directory promotion.
expected_tables <- c(
  "RF_site_blocked_metrics_v1.tsv", "RF_site_blocked_predictions_v1.tsv",
  "RF_paired_label_permutation_null_v1.tsv", "RF_fold_stable_feature_importance_v1.tsv",
  "NCM_A72_summary_v1.tsv", "NCM_A72_taxon_fit_source_v1.tsv",
  "NCM_A72_leave_one_site_out_v1.tsv", "NCM_A72_site_cluster_bootstrap_v1.tsv",
  "tNST_summary_and_bootstrap_v1.tsv", "tNST_method_sensitivity_direction_v1.tsv",
  "Network_all_tested_edges_v1.tsv", "Network_stable_edges_v1.tsv",
  "Network_summary_v1.tsv", "advanced_ecology_status_gate_v1.tsv"
)
assert_true(all(file.exists(file.path(stage, "tables", as_v2_table(expected_tables)))), "Expected output table missing")
expected_figures <- paste0(rep(c(
  "Fig_Adv_RF_site_blocked_performance_v1", "Fig_Adv_RF_feature_importance_v1",
  "Fig_Adv_NCM_A72_v1", "Fig_Adv_tNST_sensitivity_v1",
  "Fig_Adv_residual_association_networks_v1"), each = 3L), rep(c(".pdf", ".svg", ".png"), 5L))
# Build names explicitly because rep ordering above is easy to misread.
expected_figures <- unlist(lapply(as_v2_figure(c(
  "Fig_Adv_RF_site_blocked_performance_v1", "Fig_Adv_RF_feature_importance_v1",
  "Fig_Adv_NCM_A72_v1", "Fig_Adv_tNST_sensitivity_v1",
  "Fig_Adv_residual_association_networks_v1")), function(z) paste0(z, c(".pdf", ".svg", ".png"))))
assert_true(all(file.exists(file.path(stage, "figures", expected_figures))), "Expected figure bundle missing")
assert_true(all(file.info(file.path(stage, "figures", expected_figures))$size > 0), "Zero-byte figure detected")

if (!file.rename(stage, paths$out)) stop("Atomic promotion of advanced ecology output failed", call. = FALSE)
promoted <- TRUE
cat("ADVANCED ECOLOGY COMPLETE:", paths$out, "\n")
print(status)
