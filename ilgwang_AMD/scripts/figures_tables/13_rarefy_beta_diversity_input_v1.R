suppressMessages(library(vegan))
args <- commandArgs(trailingOnly=TRUE)
source_root <- if (length(args) >= 1L) args[[1]] else "/input/youre/data/"
output_root <- if (length(args) >= 2L) args[[2]] else "/output/your/data"
seq_path <- file.path(source_root, "output_data", "seqtab_final_long.csv")
out_path <- file.path(output_root, "source_data", "SFig2_rarefied_ASV_table_v1.tsv")
log_path <- file.path(output_root, "logs", "13_rarefy_beta_diversity_input_v1.log")
check_path <- file.path(output_root, "logs", "13_rarefy_beta_diversity_input_v1_checks.tsv")
if (file.exists(out_path) || file.exists(log_path) || file.exists(check_path)) stop("Refusing to overwrite existing outputs")
dir.create(dirname(out_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(log_path), recursive=TRUE, showWarnings=FALSE)
seqtab <- read.csv(seq_path, check.names=FALSE, stringsAsFactors=FALSE)
samples <- c("1st_BFP", "1st_BFW1", "1st_BFW2", "2nd_BFP", "2nd_BFW", "2nd_wt", "3rd_wt")
stopifnot(nrow(seqtab) == 663L, all(c("ASV_ID", "sequence", samples) %in% names(seqtab)), length(unique(seqtab$sequence)) == 663L)
raw <- as.matrix(seqtab[, samples])
storage.mode(raw) <- "integer"
rownames(raw) <- seqtab$ASV_ID
mat <- t(raw)
stopifnot(all(rowSums(mat) >= 88397L), min(rowSums(mat)) == 88397L, sum(mat) == 716915L)
set.seed(100)
rare <- rrarefy(mat, sample=88397L)
stopifnot(all(rowSums(rare) == 88397L))
out <- data.frame(ASV_ID=seqtab$ASV_ID, sequence=seqtab$sequence, t(rare), check.names=FALSE)
write.table(out, file=out_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
checks <- data.frame(
  check=c("ASV_rows", "samples", "raw_total_reads", "minimum_raw_depth", "rarefaction_depth", "all_rarefied_depths", "rarefied_total_reads", "seed"),
  observed=c(nrow(seqtab), nrow(mat), sum(mat), min(rowSums(mat)), 88397L, all(rowSums(rare) == 88397L), sum(rare), 100L),
  expected=c(663L, 7L, 716915L, 88397L, 88397L, TRUE, 618779L, 100L),
  status="PASS",
  stringsAsFactors=FALSE
)
write.table(checks, file=check_path, sep="\t", quote=FALSE, row.names=FALSE)
lines <- c("PASS: Rarefied beta-diversity input generated.", "Input: seqtab_final_long.csv", "Method: vegan::rrarefy, set.seed(100), target depth=88,397 per sample", sprintf("Raw total=%d; rarefied total=%d", sum(mat), sum(rare)), "Output is used only for descriptive S2 Bray-Curtis heatmap.")
writeLines(lines, log_path)
cat(paste(lines, collapse="\n"), "\n")
