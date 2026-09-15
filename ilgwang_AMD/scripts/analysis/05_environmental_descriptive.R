input_root <- "/input/youre/data/"
output_root <- "/output/your/data"
results_dir <- file.path(output_root, "results")

env <- read.delim(file.path(input_root, "metadata", "env_group_level.tsv"), stringsAsFactors = FALSE)
sample_sheet <- read.delim(file.path(input_root, "metadata", "sample_sheet.tsv"), stringsAsFactors = FALSE)
alpha <- read.csv(file.path(results_dir, "alpha_diversity_per_sample.csv"), stringsAsFactors = FALSE)
site_scores <- read.csv(file.path(results_dir, "cpcoa_site_scores.csv"), stringsAsFactors = FALSE)

cat("=== Group-level environmental measurements (n=3 groups; descriptive only) ===\n")
print(env)


merged <- merge(alpha[, c("sample","Group","Type_label","Observed","Shannon","Pielou_evenness")],
                 env, by = "Group", all.x = TRUE)
merged <- merge(merged, site_scores[, c("sample","CAP1","CAP2")], by = "sample", all.x = TRUE)
merged <- merged[order(merged$Group, merged$sample), ]

write.csv(merged, file.path(results_dir, "environmental_descriptive_by_sample.csv"), row.names = FALSE)
cat("\n=== Sample-level community metrics juxtaposed with Group-level env (descriptive) ===\n")
print(merged)


env_vars <- env[, c("Temperature","pH","DO","EC","ORP")]
cor_mat <- round(cor(env_vars, method = "pearson"), 2)
write.csv(cor_mat, file.path(results_dir, "environmental_variable_correlation_descriptive.csv"))
cat("\n=== Descriptive Pearson correlation among env variables (Group level, n=3; NO p-values reported) ===\n")
print(cor_mat)

cat("\nNOTE: No Mantel test, no correlation p-value, and no formal environment-community\n")
cat("association test is performed or reported by this script, by design. Only n=3\n")
cat("independent environmental measurements exist for 7 community samples.\n")

cat("\nDONE\n")
