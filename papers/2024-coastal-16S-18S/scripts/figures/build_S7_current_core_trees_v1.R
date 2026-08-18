#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("usage: worker <package_root> <marker>")
root <- normalizePath(args[[1]], winslash="/", mustWork=TRUE); marker <- args[[2]]
suppressPackageStartupMessages({library(Biostrings); library(DECIPHER); library(ape); library(phangorn)})
set.seed(1)
sd <- file.path(root, "source_data")
inp <- file.path(sd, sprintf("Supplementary_Figure_S7_%s_representative_core_unaligned_v1.fasta", marker))
alnfile <- file.path(sd, sprintf("Supplementary_Figure_S7_%s_representative_core_aligned_v1.fasta", marker))
treefile <- file.path(sd, sprintf("Supplementary_Figure_S7_%s_representative_core_ML_tree_v1.newick", marker))
segfile <- file.path(sd, sprintf("Supplementary_Figure_S7_%s_representative_core_tree_segments_v1.tsv", marker))
seqs <- readDNAStringSet(inp); seqs <- seqs[order(names(seqs))]
if (length(seqs) < 3L) stop("At least 3 tips are required")
aln <- AlignSeqs(seqs, processors=1, verbose=FALSE); writeXStringSet(aln, alnfile, width=80)
dat <- phyDat(as.matrix(aln), type="DNA")
d <- dist.ml(dat); nj0 <- NJ(d)
fit0 <- pml(nj0, dat)
fit <- optim.pml(fit0, model="GTR", optInv=FALSE, optGamma=FALSE, optNni=TRUE, control=pml.control(trace=0))
tr <- ladderize(fit$tree, right=FALSE); write.tree(tr, treefile)
# Deterministic rectangular phylogram coordinates for exact R/Python parity.
nt <- length(tr$tip.label); nn <- tr$Nnode; y <- numeric(nt+nn)
y[seq_len(nt)] <- seq_len(nt)
repeat {
  changed <- FALSE
  for (node in (nt+1):(nt+nn)) if (y[node] == 0) {
    kids <- tr$edge[tr$edge[,1] == node, 2]
    if (length(kids) && all(y[kids] > 0)) { y[node] <- mean(range(y[kids])); changed <- TRUE }
  }
  if (all(y[(nt+1):(nt+nn)] > 0) || !changed) break
}
rootnode <- setdiff(tr$edge[,1], tr$edge[,2])[[1]]; x <- rep(NA_real_, nt+nn); x[rootnode] <- 0
repeat {
  changed <- FALSE
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i,1]; c <- tr$edge[i,2]
    if (is.finite(x[p]) && !is.finite(x[c])) { x[c] <- x[p] + tr$edge.length[i]; changed <- TRUE }
  }
  if (all(is.finite(x)) || !changed) break
}
rows <- lapply(seq_len(nrow(tr$edge)), function(i) {
 p <- tr$edge[i,1]; c <- tr$edge[i,2]
 data.frame(Marker=marker, EdgeOrder=i, ParentNode=p, ChildNode=c, ParentX=x[p], ChildX=x[c], ParentY=y[p], ChildY=y[c], BranchLength=tr$edge.length[i], ChildIsTip=c<=nt, ChildTipLabel=if(c<=nt) tr$tip.label[c] else "", stringsAsFactors=FALSE)
})
write.table(do.call(rbind, rows), segfile, sep="\t", quote=FALSE, row.names=FALSE)
