suppressMessages(library(dada2))

input_root <- "/input/youre/data/"
output_root <- "/output/your/data"
results_dir <- file.path(output_root, "results")
ref_fasta <- file.path(input_root, "ref", "silva_nr99_v144_toSpecies_trainset.fa.gz")
tax_levels <- c("Kingdom", "Superphylum", "Phylum", "Class", "Order", "Family", "Genus", "Species")

seqtab.nochim <- readRDS(file.path(results_dir, "seqtab_nochim.rds"))

lens <- nchar(getSequences(seqtab.nochim))
cat("Length distribution before filter:\n"); print(table(lens))


keep <- lens >= 250 & lens <= 256
cat("Sequences kept:", sum(keep), "of", length(keep), "\n")
cat("Reads kept:", sum(seqtab.nochim[, keep]), "of", sum(seqtab.nochim),
    sprintf(" (%.2f%%)\n", 100 * sum(seqtab.nochim[, keep]) / sum(seqtab.nochim)))

seqtab.final <- seqtab.nochim[, keep, drop = FALSE]
saveRDS(seqtab.final, file.path(results_dir, "seqtab_final.rds"))

reads_per_sample_final <- rowSums(seqtab.final)
write.csv(data.frame(sample = rownames(seqtab.final), reads_final = reads_per_sample_final),
          file.path(results_dir, "reads_per_sample_final.csv"), row.names = FALSE)

cat("== assignTaxonomy (SILVA r144, minBoot = 80), run in 2 query batches ==\n")
set.seed(100)
n <- ncol(seqtab.final)
half1 <- seqtab.final[, 1:ceiling(n / 2), drop = FALSE]
half2 <- seqtab.final[, (ceiling(n / 2) + 1):n, drop = FALSE]

t0 <- Sys.time()
res1 <- assignTaxonomy(half1, ref_fasta, minBoot = 80, multithread = TRUE,
                        verbose = TRUE, outputBootstraps = TRUE, taxLevels = tax_levels)
res2 <- assignTaxonomy(half2, ref_fasta, minBoot = 80, multithread = TRUE,
                        verbose = TRUE, outputBootstraps = TRUE, taxLevels = tax_levels)
cat("assignTaxonomy total time:", as.numeric(difftime(Sys.time(), t0, units = "mins")), "min\n")

taxa <- rbind(res1$tax, res2$tax)
boot <- rbind(res1$boot, res2$boot)
stopifnot(identical(rownames(taxa), colnames(seqtab.final)))

saveRDS(taxa, file.path(results_dir, "taxa_silva144.rds"))
saveRDS(list(tax = taxa, boot = boot), file.path(results_dir, "taxa_silva144_with_bootstrap.rds"))
write.csv(data.frame(ASV_ID = rownames(taxa), taxa, sequence = rownames(taxa)),
          file.path(results_dir, "taxa_silva144.csv"), row.names = FALSE)

cat("Taxonomy table dim:", dim(taxa), "\n")
cat("Domain-level counts (n ASVs):\n")
print(table(taxa[, "Kingdom"], useNA = "ifany"))

cat("Genus-level NA count:", sum(is.na(taxa[, "Genus"])), "of", nrow(taxa),
    sprintf(" (%.1f%%)\n", 100 * sum(is.na(taxa[, "Genus"])) / nrow(taxa)))
cat("NOTE: see results/CRITICAL_FINDING_ferrovum_taxonomy_resolution.md - the\n")
cat("most abundant ASV (19.4%% of reads) is a 100%% sequence match to Ferrovum\n")
cat("but does not reach 80%% bootstrap support beyond Class level.\n")

cat("DONE\n")
