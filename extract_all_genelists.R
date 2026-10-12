setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(glmnet); library(parallel)})
CLIN <- c("(Intercept)","Tissue_Sex","Tissue_Age","Tissue_race","Tissue_smoke_status","Tissue_Tumor_size","Tissue_Module20")
parse_key <- function(bn, pop){
  # <set>_<cohort...>_<datatype>_<featureset...>_<subset>_<model/fs>_v<vc>.rds
  s <- sub("\\.rds$","",bn); sub("_cvSummary","",s)
}
one <- function(f){
  bn <- basename(f); kind <- if (grepl("_rf_v|_xgb_v", bn)) "tree" else "glm"
  out <- tryCatch({
    if (kind=="glm"){
      g <- readRDS(f); b <- as.matrix(coef(g, s="lambda.min")); b <- b[rownames(b)!="(Intercept)",,drop=FALSE]
      genes <- setdiff(rownames(b), CLIN); gv <- b[genes,1]
      used <- names(gv)[gv!=0]; val <- gv[used]
      list(offered=length(genes), used=used, val=round(as.numeric(val),5), metric="glm_coef")
    } else {
      x <- readRDS(f); feat <- as.character(x$features); genes <- setdiff(feat, CLIN)
      m <- x$model
      if (inherits(m,"xgb.Booster")){
        imp <- xgboost::xgb.importance(model=m)
        fi <- imp$Feature; if (length(fi) && all(grepl("^f[0-9]+$", fi))) { idx <- as.integer(sub("f","",fi))+1L; gn <- feat[idx] } else gn <- fi
        keep <- !gn %in% CLIN; gn <- gn[keep]; gain <- imp$Gain[keep]
        list(offered=length(genes), used=gn, val=round(gain,5), metric="xgb_gain")
      } else {
        im <- m$importance; col <- ncol(im); ord <- order(-im[,col]); gn <- rownames(im)[ord]
        keep <- !gn %in% CLIN; gn <- gn[keep]; vv <- im[ord,col][keep]
        list(offered=length(genes), used=gn, val=round(as.numeric(vv),5), metric="rf_import")
      }
    }
  }, error=function(e) list(offered=NA, used=character(0), val=numeric(0), metric=paste0("ERR:",conditionMessage(e))))
  out$key <- sub("\\.rds$","",bn); out
}
process_dir <- function(dir, pop){
  fs <- list.files(dir, pattern="_v[0-9]+\\.rds$", full.names=TRUE)
  cat(pop, dir, ":", length(fs), "files\n")
  res <- mclapply(fs, one, mc.cores=8)
  res
}
dirs <- list(c("res","full"), c("res_rfxgb","full"), c("res_smoking","smoke"), c("res_rfxgb_smoking","smoke"))
all <- list()
for (d in dirs) all <- c(all, lapply(process_dir(d[1], d[2]), function(r){ r$pop <- d[2]; r }))
# write long CSV (used genes only, with value)
con <- gzfile("model_gene_lists_202609.csv.gz","w")
writeLines("pop,key,metric,n_offered,rank,gene,value", con)
# compact JSON of model-used genes (cap 60 by |value|)
jcon <- file("model_gene_lists_compact.json","w"); writeLines("{", jcon); first<-TRUE
emptyfs <- c()
for (r in all){
  if (length(r$used)){
    o <- order(-abs(r$val)); g <- r$used[o]; v <- r$val[o]
    for (i in seq_along(g)) writeLines(sprintf("%s,%s,%s,%s,%d,%s,%s", r$pop, r$key, r$metric, as.character(r$offered), i, g[i], v[i]), con)
    k <- min(length(g),60)
    entry <- paste0('"', r$key, '":[', paste(sprintf('["%s",%s]', g[1:k], v[1:k]), collapse=","), ']')
    writeLines(paste0(if(first)"" else ",", entry), jcon); first<-FALSE
  } else {
    if (!is.na(r$offered) && r$offered==0) emptyfs <- c(emptyfs, r$key)
  }
}
close(con); writeLines("}", jcon); close(jcon)
cat("TOTAL configs:", length(all), "  with 0 offered genes:", length(emptyfs), "\n")
if (length(emptyfs)) { cat("sample 0-gene keys:\n"); print(head(emptyfs, 8)) }
errs <- sum(grepl("^ERR", sapply(all, function(r) r$metric)))
cat("extraction errors:", errs, "\n")
cat("EXTRACT_DONE\n")
