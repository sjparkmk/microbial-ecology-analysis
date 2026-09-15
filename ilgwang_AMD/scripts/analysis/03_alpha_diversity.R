suppressMessages(library(vegan))

set.seed(100)
input_root <- "/input/youre/data/"
output_root <- "/output/your/data"
results_dir <- file.path(output_root, "results")

seqtab.final <- readRDS(file.path(results_dir, "seqtab_final.rds"))
sample_sheet <- read.delim(file.path(input_root, "metadata", "sample_sheet.tsv"), stringsAsFactors = FALSE)

depths <- rowSums(seqtab.final)
cat("Sequencing depths (post length-filter):\n"); print(depths)
min_depth <- min(depths)
cat("Rarefying to minimum depth:", min_depth, "\n")

otu_rare <- rrarefy(seqtab.final, sample = min_depth)
storage.mode(otu_rare) <- "integer"

est <- estimateR(otu_rare)
shannon <- diversity(otu_rare, index = "shannon")
simpson <- diversity(otu_rare, index = "simpson")
invsimpson <- diversity(otu_rare, index = "invsimpson")
obs_rich <- specnumber(otu_rare)
pielou <- shannon / log(obs_rich)

alpha <- data.frame(
  sample = rownames(otu_rare),
  depth_rarefied = min_depth,
  Observed = obs_rich,
  Chao1 = est["S.chao1", ],
  Chao1_se = est["se.chao1", ],
  ACE = est["S.ACE", ],
  ACE_se = est["se.ACE", ],
  Shannon = shannon,
  Simpson = simpson,
  InvSimpson = invsimpson,
  Pielou_evenness = pielou,
  row.names = NULL
)

m <- match(alpha$sample, sample_sheet$analysis_name)
alpha$Group <- sample_sheet$Group[m]
alpha$Type_label <- sample_sheet$Type_label[m]
alpha <- alpha[, c("sample", "Group", "Type_label", setdiff(names(alpha), c("sample","Group","Type_label")))]

write.csv(alpha, file.path(results_dir, "alpha_diversity_per_sample.csv"), row.names = FALSE)
cat("\n=== Alpha diversity per sample (rarefied to", min_depth, ") ===\n")
print(alpha)


type_summary <- aggregate(cbind(Observed, Chao1, ACE, Shannon, Simpson, InvSimpson, Pielou_evenness) ~ Type_label,
                           data = alpha, FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
cat("\n=== Type-level summary (mean/sd/n; descriptive, n=2-3 per Type) ===\n")
print(type_summary)

saveRDS(otu_rare, file.path(results_dir, "otu_table_rarefied.rds"))


kw_results <- lapply(c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","Pielou_evenness"), function(v) {
  kt <- kruskal.test(alpha[[v]] ~ alpha$Type_label)
  data.frame(metric = v, KW_chisq = unname(kt$statistic), df = unname(kt$parameter), p_value = kt$p.value)
})
kw_df <- do.call(rbind, kw_results)
write.csv(kw_df, file.path(results_dir, "alpha_diversity_kruskal_wallis_exploratory.csv"), row.names = FALSE)
cat("\n=== Kruskal-Wallis by Type (EXPLORATORY ONLY, n=2/3/2, not confirmatory) ===\n")
print(kw_df)

cat("\nDONE\n")
