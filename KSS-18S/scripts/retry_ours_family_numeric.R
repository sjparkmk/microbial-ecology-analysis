# One estimator-failure retry, never a significance-selected rerun.
# Use the frozen original implementation with only the failed feature excluded.
args <- commandArgs(trailingOnly=TRUE)
arg <- function(k) args[match(k,args)+1]
destination <- arg("--out")
stopifnot(grepl("daa_ours_family_retry_v1$",destination),!dir.exists(destination))
original <- file.path(Sys.getenv("KSS_SCRIPT_ROOT"),"comparison_statistics.R")
code <- paste(readLines(original,warn=FALSE),collapse="\n")
start <- "  for(rank in names(counts))for(ds in c(\"ours\",\"d2024\",\"eukbank\",\"merged\")) {"
finish <- "  if(length(summary))write_tsv(do.call(rbind,summary),\"daa_model_summary.tsv\")"
loc1 <- regexpr(start,code,fixed=TRUE)[1]
loc2 <- regexpr(finish,code,fixed=TRUE)[1]
stopifnot(loc1>0,loc2>loc1)
replacement <- '
  excluded <- "Eukaryota|Obazoa|Opisthokonta|Fungi|Basidiomycota|Agaricomycotina|Agaricomycetes"
  selected <- md$dataset=="ours"
  x <- counts[["Family"]][,selected,drop=FALSE]
  stopifnot(excluded %in% rownames(x))
  audit <- data.frame(sample=colnames(x),excluded_feature=excluded,count=x[excluded,],
    reason="ANCOMBC2 .bias_em coefficient variance zero in original fit; raw counts are NOT constant",
    selection_basis="exact taxon named by estimator failure, not P/q/effect direction")
  write_tsv(audit,"non_estimable_feature_exclusion.tsv")
  x <- x[rownames(x)!=excluded,,drop=FALSE]
  run_fit(x,droplevels(md[selected,,drop=FALSE]),"ours_Family_planned","sample_type",NULL,.15)
'
code <- paste0(substr(code,1,loc1-1),replacement,substr(code,loc2,nchar(code)))
eval(parse(text=code),envir=.GlobalEnv)
