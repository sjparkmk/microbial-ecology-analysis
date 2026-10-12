# Official numerical community and compositional analyses; no graphics.
suppressPackageStartupMessages(library(vegan))
VERSION <- "1.0.0"
args <- commandArgs(trailingOnly=TRUE)
if ("--help" %in% args) {cat("--input DIR --out DIR --mode community|daa|aldex|check --version\n");quit(status=0)}
if ("--version" %in% args) {cat(VERSION,"\n");quit(status=0)}
arg <- function(key, default=NULL) {i<-match(key,args);if(is.na(i)) default else args[i+1]}
input <- arg("--input"); out <- arg("--out"); mode <- arg("--mode","check")
stopifnot(!is.null(input),!is.null(out))
SEED <- 20260923L
readtab <- function(path) read.delim(path, row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
md <- readtab(file.path(input,"metadata.tsv"))
md$sample_type<-factor(md$sample_type,levels=c("land_sediment","land_soil","land_freshwater","marine_sediment","marine_water"))
md$dataset<-factor(md$dataset,levels=c("ours","d2024","eukbank"))
md$cluster<-factor(md$cluster);md$site<-factor(md$site);md$season<-factor(md$season)
counts <- lapply(c("Genus","Family"),function(r) as.matrix(readtab(file.path(input,paste0("counts_",r,".tsv")))))
names(counts)<-c("Genus","Family")
for (x in counts) stopifnot(identical(colnames(x),rownames(md)),all(x>=0),all(x==floor(x)),all(colSums(x)>0))
if(mode=="check") {cat("PASS input alignment, nonnegative integer counts, positive totals\n");quit(status=0)}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,name) {path<-file.path(out,name);if(file.exists(path))stop("Refusing existing output: ",path);write.table(x,path,sep="\t",row.names=FALSE,quote=FALSE,na="NA")}
capture.output(sessionInfo(),file=file.path(out,paste0("sessionInfo_",mode,".txt")))
cat("Mode",mode,"R",R.version.string,"seed",SEED,"\n")

if(mode=="community") {
  results<-list();dispersions<-list();status<-list()
  for(rank in names(counts)) {
    x<-counts[[rank]]; rel<-sweep(x,2,colSums(x),"/"); d<-vegdist(t(rel),method="bray")
    mat<-as.matrix(d)
    write_tsv(data.frame(sample=rownames(mat),mat,check.names=FALSE),paste0(rank,"_bray_distances.tsv"))
    pc<-wcmdscale(d,k=2,eig=TRUE,add="lingoes")
    pts<-data.frame(sample=rownames(pc$points),PCoA1=pc$points[,1],PCoA2=pc$points[,2])
    write_tsv(pts,paste0(rank,"_pcoa.tsv"))
    write_tsv(data.frame(axis=seq_along(pc$eig),eigenvalue=pc$eig,
       fraction_of_positive=pc$eig/sum(pc$eig[pc$eig>0]),lingoes_constant=pc$ac),paste0(rank,"_pcoa_eigenvalues.tsv"))
    for(dataset in c("combined","ours","d2024","eukbank")) {
      keep<-if(dataset=="combined")rep(TRUE,nrow(md)) else md$dataset==dataset
      meta<-droplevels(md[keep,,drop=FALSE]);dist_sub<-as.dist(mat[keep,keep]);
      contrasts<-if(dataset=="combined")c("unrestricted","within_dataset") else if(dataset=="d2024")c("site_season_adjusted") else if(dataset=="eukbank")c("unrestricted","within_project") else c("unrestricted")
      for(design in contrasts) {
        perm<-9999
        if(design=="within_dataset")perm<-permute::how(nperm=9999,blocks=meta$dataset)
        if(design=="site_season_adjusted")perm<-permute::how(nperm=9999,blocks=interaction(meta$site,meta$season,drop=TRUE))
        if(design=="within_project")perm<-permute::how(nperm=9999,blocks=meta$cluster)
        form<-if(dataset=="combined")dist_sub~dataset+sample_type else if(dataset=="d2024")dist_sub~site+season+sample_type else dist_sub~sample_type
        set.seed(SEED)
        ans<-tryCatch(adonis2(form,data=meta,permutations=perm,by="margin"),error=function(e)e)
        if(inherits(ans,"error")) {status[[length(status)+1]]<-data.frame(rank,dataset,design,status=conditionMessage(ans));next}
        a<-as.data.frame(ans);a$term<-rownames(a);a$rank<-rank;a$dataset<-dataset;a$design<-design;a$n<-nrow(meta)
        a$permutation_scope<-design
        a$interpretation<-if(design=="unrestricted" && dataset %in% c("combined","eukbank"))"exploratory; source clustering not controlled" else "restricted association, not causal"
        # Restricted design cannot test a variable that never changes within its blocks.
        if(design=="within_dataset")a[a$term=="dataset","Pr(>F)"]<-NA_real_
        if(design=="site_season_adjusted")a[a$term%in%c("site","season"),"Pr(>F)"]<-NA_real_
        if(design=="within_project" && all(tapply(meta$sample_type,meta$cluster,function(z)length(unique(z)))<=1)) {
          a[a$term=="sample_type","Pr(>F)"]<-NA_real_
          a$interpretation<-"not identifiable conditional on originating project: no within-project habitat variation"
        }
        results[[length(results)+1]]<-a
      }
      for(factor_name in if(dataset=="combined")c("sample_type","dataset") else "sample_type") {
        set.seed(SEED)
        z<-tryCatch({bd<-betadisper(dist_sub,meta[[factor_name]],type="median",bias.adjust=TRUE,add="lingoes");pp<-permutest(bd,permutations=9999);data.frame(rank,dataset,factor=factor_name,n=nrow(meta),F=pp$tab[1,"F"],p=pp$tab[1,"Pr(>F)"],note="PERMDISP diagnostic; unrestricted permutations")},error=function(e)data.frame(rank,dataset,factor=factor_name,n=nrow(meta),F=NA,p=NA,note=conditionMessage(e)))
        dispersions[[length(dispersions)+1]]<-z
      }
    }
  }
  zz<-do.call(rbind,results);zz$q_BH_all_reported_terms<-p.adjust(zz[["Pr(>F)"]],"BH")
  write_tsv(zz,"permanova.tsv");write_tsv(do.call(rbind,dispersions),"permdisp.tsv")
  if(length(status))write_tsv(do.call(rbind,status),"community_failures.tsv")
  cat("COMMUNITY COMPLETE\n")
}

if(mode=="daa") {
  suppressPackageStartupMessages(library(ANCOMBC))
  summary<-list()
  run_fit<-function(x,meta,label,formula,random,prv,group="sample_type") {
    destination<-file.path(out,label);dir.create(destination,recursive=TRUE,showWarnings=FALSE)
    if(file.exists(file.path(destination,"completed.txt"))) {cat("Already complete",label,"\n");return(invisible(NULL))}
    meta<-droplevels(meta[colnames(x),,drop=FALSE]);x<-x[rowSums(x)>0,,drop=FALSE]
    storage.mode(x)<-"integer"
    cat("START",label,"n",ncol(x),"features",nrow(x),"formula",formula,"random",random,"\n");flush.console()
    fixed<-model.matrix(as.formula(paste("~",formula)),data=meta)
    if(qr(fixed)$rank<ncol(fixed))stop("Rank-deficient fixed effects: ",label)
    warnings<-character()
    set.seed(SEED)
    fit<-tryCatch(withCallingHandlers(ancombc2(data=x,meta_data=meta,taxa_are_rows=TRUE,
      fix_formula=formula,rand_formula=random,p_adj_method="BH",pseudo_sens=TRUE,
      prv_cut=prv,lib_cut=0,s0_perc=.05,group=group,struc_zero=TRUE,neg_lb=TRUE,
      alpha=.05,n_cl=1,verbose=FALSE,global=nlevels(meta[[group]])>2,
      pairwise=nlevels(meta[[group]])>2,dunnet=FALSE,trend=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
    config<-data.frame(label,n=ncol(x),input_features=nrow(x),prv_cut=prv,formula,random=if(is.null(random))"none" else random,package=as.character(packageVersion("ANCOMBC")),seed=SEED)
    write.table(config,file.path(destination,"model_contract.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    writeLines(unique(warnings),file.path(destination,"warnings.txt"))
    if(inherits(fit,"error")) {
      writeLines(conditionMessage(fit),file.path(destination,"failure.txt"));cat("FAILED",label,conditionMessage(fit),"\n")
      summary[[length(summary)+1]]<<-data.frame(label,status="failed",tested=NA,robust=NA,detail=conditionMessage(fit));return(invisible(NULL))
    }
    saveRDS(fit,file.path(destination,"fit.rds"))
    for(item in c("res","res_global","res_pair","zero_ind"))if(!is.null(fit[[item]])) {
      tab<-as.data.frame(fit[[item]]);if(!"taxon" %in% names(tab))tab<-data.frame(taxon=rownames(tab),tab,check.names=FALSE)
      write.table(tab,file.path(destination,paste0(item,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    }
    res<-fit$res
    robust_cols<-grep("^diff_robust_",names(res),value=TRUE)
    if(length(robust_cols)) nr<-sum(as.matrix(res[,robust_cols,drop=FALSE]),na.rm=TRUE) else {
      terms<-sub("^q_","",grep("^q_",names(res),value=TRUE));terms<-setdiff(terms,"(Intercept)")
      nr<-sum(vapply(terms,function(t) {qc<-paste0("q_",t);sc<-paste0("passed_ss_",t);if(!sc%in%names(res))return(0);sum(res[[qc]]<.05 & res[[sc]],na.rm=TRUE)},numeric(1)))
    }
    summary[[length(summary)+1]]<<-data.frame(label,status="computed",tested=nrow(res),robust=nr,detail="finite-effect results; structural zeros separate; exploratory candidates only")
    writeLines("numeric fit exported; interpretation and QA remain",file.path(destination,"completed.txt"))
    cat("DONE",label,"tested",nrow(res),"robust terms",nr,"\n");flush.console()
  }
  for(rank in names(counts))for(ds in c("ours","d2024","eukbank","merged")) {
    sel<-if(ds=="merged")rep(TRUE,nrow(md)) else md$dataset==ds
    x<-counts[[rank]][,sel,drop=FALSE];meta<-droplevels(md[sel,,drop=FALSE])
    prv<-c(ours=.15,d2024=.25,eukbank=.10,merged=.08)[ds]
    formula<-if(ds=="merged")"sample_type + dataset" else if(ds=="d2024")"sample_type + season" else "sample_type"
    run_fit(x,meta,paste(ds,rank,"planned",sep="_"),formula,NULL,prv)
    if(ds %in% c("d2024","eukbank","merged"))run_fit(x,meta,paste(ds,rank,"clustered",sep="_"),formula,"(1 | cluster)",prv)
  }
  x<-as.matrix(readtab(file.path(input,"ours_ASV_counts.tsv")));om<-md[colnames(x),,drop=FALSE];om$exposure<-factor(om$exposure,levels=c("exposed","submerged"))
  run_fit(x,om,"ours_ASV_exposure","exposure",NULL,.20,"exposure")
  if(length(summary))write_tsv(do.call(rbind,summary),"daa_model_summary.tsv")
  cat("DAA COMPLETE\n")
}

if(mode=="aldex") {
  suppressPackageStartupMessages(library(ALDEx2));summary<-list()
  for(rank in names(counts))for(ds in c("ours","d2024","eukbank")) {
    m<-droplevels(md[md$dataset==ds,,drop=FALSE]);pairs<-combn(levels(m$sample_type),2,simplify=FALSE)
    for(pair in pairs) {
      sm<-m[m$sample_type%in%pair,,drop=FALSE];sm<-sm[order(sm$sample_type),,drop=FALSE]
      x<-counts[[rank]][,rownames(sm),drop=FALSE]
      cut<-c(ours=.15,d2024=.25,eukbank=.10)[ds]
      keep<-rowMeans(x>0)>=cut & rowSums(x)>=10
      x<-x[keep,,drop=FALSE];label<-paste(ds,rank,pair[2],"vs",pair[1],sep="_")
      cat("ALDEX",label,"n",ncol(x),"features",nrow(x),"\n");flush.console()
      set.seed(SEED)
      a<-tryCatch(aldex(x,as.character(sm$sample_type),mc.samples=128,test="t",effect=TRUE,denom="all",verbose=FALSE),error=function(e)e)
      if(inherits(a,"error")){summary[[length(summary)+1]]<-data.frame(label,status="failed",features=nrow(x),significant=NA,note=conditionMessage(a));next}
      a<-data.frame(feature=rownames(a),a,check.names=FALSE)
      write_tsv(a,paste0(label,".tsv"))
      summary[[length(summary)+1]]<-data.frame(label,status="computed",features=nrow(x),significant=sum(a$we.eBH<.05,na.rm=TRUE),note="two-group compositional sensitivity; no site/project random effect; not validation")
    }
  }
  x<-as.matrix(readtab(file.path(input,"ours_ASV_counts.tsv")));m<-md[colnames(x),,drop=FALSE];keep<-rowMeans(x>0)>=.20
  set.seed(SEED);a<-tryCatch(aldex(x[keep,,drop=FALSE],m$exposure,mc.samples=128,test="t",effect=TRUE,denom="all",verbose=FALSE),error=function(e)e)
  if(!inherits(a,"error"))write_tsv(data.frame(feature=rownames(a),a,check.names=FALSE),"ours_ASV_exposure.tsv")
  write_tsv(do.call(rbind,summary),"aldex_summary.tsv")
  cat("ALDEX COMPLETE\n")
}
