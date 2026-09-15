suppressMessages({
  library(dada2)
})

set.seed(100)

input_root <- "/input/youre/data/"
output_root <- "/output/your/data"
trimmed_dir <- file.path(input_root, "trimmed")
filt_dir <- file.path(output_root, "filtered")
results_dir <- file.path(output_root, "results")
dir.create(filt_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)

sample_sheet <- read.delim(file.path(input_root, "metadata", "sample_sheet.tsv"), stringsAsFactors = FALSE)
samples <- sample_sheet$analysis_name

fnFs <- file.path(trimmed_dir, paste0(samples, "_1.fastq.gz"))
fnRs <- file.path(trimmed_dir, paste0(samples, "_2.fastq.gz"))
names(fnFs) <- samples
names(fnRs) <- samples
stopifnot(all(file.exists(fnFs)), all(file.exists(fnRs)))

filtFs <- file.path(filt_dir, paste0(samples, "_1_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(samples, "_2_filt.fastq.gz"))
names(filtFs) <- samples
names(filtRs) <- samples

cat("== filterAndTrim ==\n")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      truncLen = c(220, 180),
                      maxN = 0, maxEE = c(2, 2), truncQ = 2,
                      rm.phix = TRUE, compress = TRUE, multithread = TRUE)
print(out)
write.csv(out, file.path(results_dir, "track_01_filterAndTrim.csv"))

cat("== learnErrors (fwd) ==\n")
errF <- learnErrors(filtFs, multithread = TRUE)
cat("== learnErrors (rev) ==\n")
errR <- learnErrors(filtRs, multithread = TRUE)
saveRDS(errF, file.path(results_dir, "errF.rds"))
saveRDS(errR, file.path(results_dir, "errR.rds"))

cat("== dada (fwd) ==\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
cat("== dada (rev) ==\n")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

cat("== mergePairs ==\n")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
saveRDS(mergers, file.path(results_dir, "mergers.rds"))

cat("== makeSequenceTable ==\n")
seqtab <- makeSequenceTable(mergers)
cat("Dim seqtab:", dim(seqtab), "\n")
cat("Sequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

cat("== removeBimeraDenovo (consensus) ==\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
cat("Dim seqtab.nochim:", dim(seqtab.nochim), "\n")
cat("Fraction reads kept after chimera removal:", sum(seqtab.nochim) / sum(seqtab), "\n")

saveRDS(seqtab, file.path(results_dir, "seqtab.rds"))
saveRDS(seqtab.nochim, file.path(results_dir, "seqtab_nochim.rds"))


getN <- function(x) sum(getUniques(x))
track <- data.frame(
  sample = samples,
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)
write.csv(track, file.path(results_dir, "read_tracking_core.csv"), row.names = FALSE)
print(track)

cat("DONE\n")
