# A1 only: pre-existing distances, no abundance/model-prediction reconstruction.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L, args[1] %in% c('--check','--run'))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(permute))
src <- file.path(Sys.getenv('KSS_INPUT_ROOT','/input/ur/data'),'comparative')
out <- file.path(Sys.getenv('KSS_OUTPUT_ROOT','/output/ur/data'),'a1_sensitivity')
md <- read.delim(file.path(src,'inputs/metadata.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(md)==140L,!anyDuplicated(md$sample),all(nzchar(md$sample)),
          sum(md$dataset=='ours')==15L,sum(md$dataset=='d2024')==75L,sum(md$dataset=='eukbank')==50L)
md$dataset <- factor(md$dataset,levels=c('ours','d2024','eukbank'))
md$sample_type <- factor(md$sample_type,levels=c('land_sediment','land_soil','land_freshwater','marine_sediment','marine_water'))
md$fine_block <- factor(ifelse(md$dataset=='ours','KSS_independent_locations',
    ifelse(md$dataset=='d2024',paste('Coastal',md$site,md$season,sep=':'),paste('EukBank',md$project,sep=':'))))
stopifnot(!anyNA(md$sample_type),!anyNA(md$dataset),!anyNA(md$fine_block))
xb <- model.matrix(~fine_block,data=md)
xf <- model.matrix(~fine_block+sample_type,data=md)
rb <- qr(xb)$rank; rf <- qr(xf)$rank
stopifnot(rf-rb==2L, nrow(md)-rf>0)
tab <- as.data.frame.matrix(table(md$fine_block,md$sample_type))
tab$block <- rownames(tab);tab$n <- rowSums(tab[,levels(md$sample_type),drop=FALSE]);tab$n_types <- rowSums(tab[,levels(md$sample_type),drop=FALSE]>0)
dists <- list()
for (rank in c('Genus','Family')) {
  a <- as.matrix(read.delim(file.path(src,'community',paste0(rank,'_bray_distances.tsv')),row.names=1,check.names=FALSE))
  stopifnot(identical(rownames(a),md$sample),identical(colnames(a),md$sample),
    all(is.finite(a)),all(a>=0),max(abs(a-t(a)))<1e-10,max(abs(diag(a)))<1e-10)
  dists[[rank]] <- as.dist(a)
}
cat('PRECHECK PASS\n');cat('N=',nrow(md),' blocks=',nlevels(md$fine_block),' informative_blocks=',sum(tab$n_types>1),
    ' block_rank=',rb,' full_rank=',rf,' within_block_df=',rf-rb,' residual_df=',nrow(md)-rf,'\n')
print(tab[tab$n_types>1,])
if (args[1]=='--check') quit(status=0)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
put <- function(x,name) write.table(x,file.path(out,name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
put(md,'aligned_metadata_with_fine_blocks.tsv');put(tab,'block_sample_type_counts.tsv')
put(data.frame(n=140,blocks=nlevels(md$fine_block),informative_blocks=sum(tab$n_types>1),block_rank=rb,
    full_rank=rf,within_block_df=rf-rb,residual_df=nrow(md)-rf,seed=20261007,permutations=9999),'design_rank_check.tsv')
set.seed(20261007)
perms <- shuffleSet(nrow(md),nset=9999,control=how(nperm=9999,blocks=md$fine_block))
stopifnot(nrow(perms)==9999L,all(apply(perms,1,function(i)all(md$fine_block[i]==md$fine_block))))
stopifnot(any(apply(perms,1,function(i)any(md$sample_type[i]!=md$sample_type))))
saveRDS(perms,file.path(out,'permutations_fine_blocks_seed20261007.rds'))
primary <- list();secondary <- list();disp <- list()
for (rank in names(dists)) {
  distance <- dists[[rank]]
  fit <- dbrda(distance ~ sample_type + Condition(fine_block),data=md)
  stopifnot(fit$CCA$rank==2L)
  test <- anova(fit,permutations=perms,model='reduced',parallel=1)
  primary[[rank]] <- data.frame(rank=rank,model='distance ~ sample_type + Condition(fine_block)',
    df=test$Df[1],residual_df=test$Df[2],sum_squares=fit$CCA$tot.chi,total_sum_squares=fit$tot.chi,
    conditioned_sum_squares=fit$pCCA$tot.chi,residual_sum_squares=fit$CA$tot.chi,
    R2_total=fit$CCA$tot.chi/fit$tot.chi,
    R2_after_conditioning=fit$CCA$tot.chi/(fit$CCA$tot.chi+fit$CA$tot.chi),
    pseudo_F=test$F[1],P=test$`Pr(>F)`[1],permutations=nrow(perms),seed=20261007)
  orig <- adonis2(distance ~ dataset + sample_type,data=md,permutations=perms,by='margin',parallel=1)
  rr <- orig['sample_type',]
  secondary[[rank]] <- data.frame(rank=rank,model='distance ~ dataset + sample_type; margin',df=rr$Df,
    sum_squares=rr$SumOfSqs,R2_total=rr$R2,pseudo_F=rr$F,P=rr$`Pr(>F)`,permutations=nrow(perms))
  bd <- betadisper(distance,md$sample_type,type='median',bias.adjust=TRUE,add='lingoes')
  btest <- permutest(bd,permutations=perms,parallel=1)
  disp[[rank]] <- data.frame(rank=rank,geometry='median; bias.adjust=TRUE; Lingoes',df=btest$tab$Df[1],
    F=btest$tab$F[1],P=btest$tab$`Pr(>F)`[1],permutations=nrow(perms))
  saveRDS(list(partial_fit=fit,partial_test=test,original_formula_sensitivity=orig,
    dispersion=bd,dispersion_test=btest),file.path(out,paste0(rank,'_tests.rds')))
}
primary <- do.call(rbind,primary);primary$P_Holm_two_ranks <- p.adjust(primary$P,'holm')
secondary <- do.call(rbind,secondary);secondary$P_Holm_two_ranks <- p.adjust(secondary$P,'holm')
disp <- do.call(rbind,disp);disp$P_Holm_two_ranks <- p.adjust(disp$P,'holm')
put(primary,'A1_partial_within_block_sensitivity.tsv');put(secondary,'A1_original_formula_fine_permutations.tsv');put(disp,'A1_matched_restricted_dispersion.tsv')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(primary);print(secondary);print(disp)
cat('A1 COMPLETE; all outcomes retained\n')
