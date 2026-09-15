suppressMessages(library(vegan))
set.seed(100)
input_root <- "/input/youre/data/"
output_root <- "/output/your/data"
results_dir <- file.path(output_root, "results")

seqtab.final <- readRDS(file.path(results_dir, "seqtab_final.rds"))
sample_sheet <- read.delim(file.path(input_root, "metadata", "sample_sheet.tsv"), stringsAsFactors = FALSE)
m <- match(rownames(seqtab.final), sample_sheet$analysis_name)
Type <- factor(sample_sheet$Type_label[m])

run_one <- function(mat, sqrt_dist) {
  bc <- vegdist(mat, method = "bray")
  if (sqrt_dist) bc <- sqrt(bc)
  cap <- capscale(bc ~ Type)
  r2 <- RsquareAdj(cap)
  ev <- cap$CA$eig
  neg_inertia <- if (!is.null(ev)) sum(ev[ev < 0]) else NA
  set.seed(100)
  pt <- anova.cca(cap, permutations = 999, by = NULL)
  data.frame(R2 = r2$r.squared, adj_R2 = r2$adj.r.squared,
             negative_eigenvalue_inertia = neg_inertia,
             F_model = pt$F[1], p_model_999perm = pt$`Pr(>F)`[1])
}

raw <- seqtab.final
relabund <- sweep(raw, 1, rowSums(raw), "/")
min_depth <- min(rowSums(raw))
set.seed(100)
rar <- rrarefy(raw, min_depth)

out <- rbind(
  cbind(input = "unrarefied_counts", sqrt_distance = FALSE, run_one(raw, FALSE)),
  cbind(input = "unrarefied_counts", sqrt_distance = TRUE,  run_one(raw, TRUE)),
  cbind(input = "relative_abundance", sqrt_distance = FALSE, run_one(relabund, FALSE)),
  cbind(input = "relative_abundance", sqrt_distance = TRUE,  run_one(relabund, TRUE)),
  cbind(input = paste0("rarefied", min_depth), sqrt_distance = FALSE, run_one(rar, FALSE)),
  cbind(input = paste0("rarefied", min_depth), sqrt_distance = TRUE,  run_one(rar, TRUE))
)
print(out)
write.csv(out, file.path(results_dir, "cpcoa_distance_sensitivity_sonnet_own_data.csv"), row.names = FALSE)
cat("\nMin depth used for rarefaction:", min_depth, "\n")
