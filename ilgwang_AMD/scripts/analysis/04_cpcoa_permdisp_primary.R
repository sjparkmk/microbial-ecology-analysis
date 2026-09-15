suppressPackageStartupMessages(library(vegan))
args <- commandArgs(trailingOnly = TRUE)
input_root <- if (length(args) >= 1L) args[[1]] else "/input/youre/data/"
output_root <- if (length(args) >= 2L) args[[2]] else "/output/your/data"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
seqtab <- read.csv(file.path(input_root, "output_data", "seqtab_final_long.csv"), check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(file.path(input_root, "input_data", "sample_sheet.tsv"), stringsAsFactors = FALSE)
samples <- metadata$analysis_name
stopifnot(length(samples) == 7L, all(samples %in% colnames(seqtab)))
counts <- t(as.matrix(seqtab[, samples, drop = FALSE]))
storage.mode(counts) <- "integer"
rownames(counts) <- samples
stopifnot(min(rowSums(counts)) == 88397L)
set.seed(100)
rarefied <- rrarefy(counts, sample = 88397L)
bray <- vegdist(rarefied, method = "bray")
sqrt_bray <- sqrt(bray)
model_data <- data.frame(Type = factor(metadata$Type[match(rownames(rarefied), metadata$analysis_name)]))
rownames(model_data) <- rownames(rarefied)
fit <- capscale(sqrt_bray ~ Type, data = model_data)
r2 <- RsquareAdj(fit)
set.seed(100)
model_test <- anova.cca(fit, permutations = 999, by = NULL)
dispersion <- betadisper(sqrt_bray, model_data$Type, type = "median", bias.adjust = FALSE)
set.seed(100)
dispersion_test <- permutest(dispersion, permutations = 999)
site <- as.data.frame(scores(fit, display = "sites"))
site$Sample <- rownames(site)
site$Type <- model_data[site$Sample, "Type"]
site$Group <- metadata$Group[match(site$Sample, metadata$analysis_name)]
site <- site[, c("Sample", "Type", "Group", "CAP1", "CAP2")]
summary <- data.frame(
  R2 = unname(r2$r.squared),
  adjusted_R2 = unname(r2$adj.r.squared),
  pseudo_F = unname(model_test$F[1]),
  model_df = unname(model_test$Df[1]),
  residual_df = unname(model_test$Df[2]),
  model_P = unname(model_test$`Pr(>F)`[1]),
  CPCoA_requested_permutations = 999L,
  CPCoA_nonidentity_sample_permutation_space = factorial(nrow(rarefied)) - 1L,
  PERMDISP_F = unname(dispersion_test$tab[1, "F"]),
  PERMDISP_P = unname(dispersion_test$tab[1, "Pr(>F)"]),
  PERMDISP_requested_permutations = 999L,
  rarefaction_depth = 88397L,
  rarefaction_seed = 100L
)
write.csv(summary, file.path(output_root, "cpcoa_permdisp_primary.csv"), row.names = FALSE)
write.table(site, file.path(output_root, "cpcoa_site_scores_primary.tsv"), sep = "	", row.names = FALSE, quote = FALSE)
write.table(as.matrix(bray), file.path(output_root, "bray_curtis_rarefied_matrix.tsv"), sep = "	", quote = FALSE, col.names = NA)
writeLines(capture.output(sessionInfo()), file.path(output_root, "cpcoa_session_info.txt"))
print(summary)
