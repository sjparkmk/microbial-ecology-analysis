# Existing official vegan implementation; no DAA, bootstraps or source writes.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
tab <- normalizePath(args[1],winslash='/'); inputs <- normalizePath(args[2],winslash='/')
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(permute))
md <- read.delim(file.path(inputs,'metadata.tsv'),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
md$dataset <- factor(md$dataset); md$sample_type <- factor(md$sample_type)
md$log10_reads <- log10(md$raw_eukaryote_reads)
set.seed(20261001)
control <- how(nperm=9999,blocks=md$dataset)
perms <- shuffleSet(nrow(md),nset=9999,control=control)
stopifnot(nrow(perms)==9999)
for (i in seq_len(nrow(perms))) stopifnot(all(md$dataset==md$dataset[perms[i,]]))
saveRDS(perms,file.path(dirname(tab),'permutations_within_dataset_seed20261001.rds'))
results <- list(); dispersion <- list(); depth_terms <- list(); idx <- 0
for (rank in c('Genus','Family')) {
  for (representation in c('raw_relative','mean_14000','mean_10000')) {
    D <- as.matrix(read.delim(file.path(tab,paste0('Depth_distance_',rank,'_',representation,'.tsv')),row.names=1,check.names=FALSE))
    stopifnot(identical(rownames(D),rownames(md)),identical(colnames(D),rownames(md)),all(is.finite(D)),max(abs(D-t(D)))<1e-10)
    d <- as.dist(D)
    fit <- adonis2(d ~ dataset + sample_type,data=md,permutations=perms,by='margin')
    idx <- idx+1
    tmp <- data.frame(rank=rank,representation=representation,term=rownames(fit),fit,check.names=FALSE)
    results[[idx]] <- tmp
    bd <- betadisper(d,md$sample_type,type='median',bias.adjust=TRUE,add='lingoes')
    pd <- permutest(bd,permutations=perms)
    dispersion[[idx]] <- data.frame(rank=rank,representation=representation,n=140,F=pd$tab[1,'F'],P=pd$tab[1,'Pr(>F)'],restriction='within_dataset',geometry='median_bias_adjusted_Lingoes')
    if (representation=='raw_relative') {
      fdepth <- adonis2(d ~ dataset + sample_type + log10_reads,data=md,permutations=perms,by='margin')
      depth_terms[[rank]] <- data.frame(rank=rank,term=rownames(fdepth),fdepth,check.names=FALSE)
    }
    cat(rank,representation,'done\n'); flush.console()
  }
}
res <- do.call(rbind,results); res$q_sample_type_family <- NA_real_
is_habitat <- res$term=='sample_type'
res$q_sample_type_family[is_habitat] <- p.adjust(res[is_habitat,'Pr(>F)'],method='BH')
write.table(res,file.path(tab,'Depth_PERMANOVA_all_terms.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
write.table(do.call(rbind,dispersion),file.path(tab,'Depth_PERMDISP.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,depth_terms),file.path(tab,'Depth_covariate_diagnostic.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dirname(tab),'environment_R_v1.txt'))
cat('VEGAN STATISTICS COMPLETE\n')
