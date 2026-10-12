setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(glmnet); library(parallel)})
TH <- c(Inf, seq(1,0,by=-0.01), -Inf)   # probability-cutoff ladder (103), same as pooled ROC
dcache <- new.env(parent=emptyenv())
get_data <- function(rds){ if(exists(rds,envir=dcache,inherits=FALSE)) return(get(rds,envir=dcache))
  v <- if(file.exists(rds)) readRDS(rds) else NULL; assign(rds,v,envir=dcache); v }
roc_from <- function(prob, y){
  ok <- is.finite(prob)&is.finite(y); prob<-prob[ok]; y<-y[ok]
  np<-sum(y==1); nn<-sum(y==0); if(np==0||nn==0) return(NULL)
  tpr <- sapply(TH, function(t) sum(prob>=t & y==1)/np)
  fpr <- sapply(TH, function(t) sum(prob>=t & y==0)/nn)
  auc <- abs(sum(diff(fpr)*(head(tpr,-1)+tail(tpr,-1))/2))
  list(fpr=paste(sprintf("%.4f",fpr),collapse="|"), tpr=paste(sprintf("%.4f",tpr),collapse="|"), auc=round(auc,4))
}
CLIN <- c("(Intercept)","Tissue_Sex","Tissue_Age","Tissue_race","Tissue_smoke_status","Tissue_Tumor_size","Tissue_Module20")
insample <- function(rds_model, data_rds, is_tree){
  if(!file.exists(rds_model)) return(NULL)
  d <- get_data(data_rds); if(is.null(d)) return(NULL)
  y <- as.numeric(as.character(d[["recurrence"]]))
  m <- readRDS(rds_model)
  prob <- tryCatch({
    if(is_tree){
      feats <- as.character(m$features); feats <- feats[feats %in% names(d)]
      mm <- m$model
      if(inherits(mm,"xgb.Booster")) predict(mm, as.matrix(d[, feats, drop=FALSE]))
      else predict(mm, d[, feats, drop=FALSE], type="prob")[,"1"]
    } else {
      b <- as.matrix(coef(m, s="lambda.min")); fn <- setdiff(rownames(b), "(Intercept)")
      fn <- fn[fn %in% names(d)]
      as.numeric(predict(m, newx=as.matrix(d[, fn, drop=FALSE]), s="lambda.min", type="response"))
    }
  }, error=function(e) NULL)
  if(is.null(prob)||length(prob)!=length(y)) return(NULL)
  roc_from(prob, y)
}
proc <- function(agg, resdir, md_base, pop, is_glm){
  a <- read.csv(agg, stringsAsFactors=FALSE)
  a <- if(is_glm) a[a$featureset %in% c("fpkm","fpkm_clinics"),] else a[a$model %in% c("rf","xgb") & a$featureset %in% c("fpkm","fpkm_clinics"),]
  cat(pop, resdir, nrow(a), "configs\n")
  res <- mclapply(seq_len(nrow(a)), function(i){ r <- a[i,]
    key <- if(is_glm) sprintf("%s_%s_%s_%s_%s_glm_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$fsmethod,r$varcut) else sprintf("%s_%s_%s_%s_%s_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$model,r$varcut)
    rm_file <- file.path(resdir, paste0(key,".rds"))
    data_rds <- file.path(md_base, r$set, r$cohort, r$datatype, paste0(r$featureset,"_",r$subset,".RDS"))
    v <- insample(rm_file, data_rds, !is_glm); if(is.null(v)) return(NULL)
    data.frame(pop=pop, key=key, auc_insample=v$auc, fpr=v$fpr, tpr=v$tpr, stringsAsFactors=FALSE)
  }, mc.cores=8)
  do.call(rbind, res[!sapply(res,is.null)])
}
out <- rbind(
  proc("results_fsmethod_comparison.csv","res","data/model_data","full",TRUE),
  proc("rf_xgb_auc_202607.csv","res_rfxgb","data/model_data","full",FALSE),
  proc("results_fsmethod_comparison_smoking.csv","res_smoking","data/model_data_smoking","smoke",TRUE),
  proc("rf_xgb_auc_smoking.csv","res_rfxgb_smoking","data/model_data_smoking","smoke",FALSE))
write.csv(out, "roc_insample_202609.csv", row.names=FALSE)
cat("wrote roc_insample_202609.csv rows", nrow(out), "\n"); cat("INSAMPLE_DONE\n")
