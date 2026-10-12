setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(data.table); library(parallel)})
TH <- c(Inf, seq(1,0,by=-0.01), -Inf)   # CALIBRATED probability-cutoff ladder (103)
labcache <- new.env(parent=emptyenv())
getlab <- function(rds){ if(exists(rds,envir=labcache,inherits=FALSE)) return(get(rds,envir=labcache))
  v <- if(file.exists(rds)){ m<-readRDS(rds); y<-as.numeric(as.character(m[,"recurrence"])); names(y)<-rownames(m); y } else NULL
  assign(rds,v,envir=labcache); v }
calib_roc <- function(pred_file, rds){
  if(!file.exists(pred_file)) return(NULL)
  y <- getlab(rds); if(is.null(y)) return(NULL)
  d <- tryCatch(fread(pred_file, showProgress=FALSE), error=function(e) NULL); if(is.null(d)||!nrow(d)) return(NULL)
  # per-subject out-of-fold mean predicted probability (honest: each pred from a model that didn't see that subject)
  ps <- d[, .(prob=mean(prob, na.rm=TRUE)), by=subj]
  ps$y <- y[as.character(ps$subj)]; ps <- ps[is.finite(ps$y) & is.finite(ps$prob)]
  np <- sum(ps$y==1); nn <- sum(ps$y==0); if(np<2||nn<2||nrow(ps)<5) return(NULL)
  # isotonic calibration fit on the OOF scores (monotone: preserves ranking/AUC, fixes the probability axis)
  o <- order(ps$prob); ir <- isoreg(ps$prob[o], ps$y[o]); calib <- numeric(nrow(ps)); calib[o] <- ir$yf
  # ROC on calibrated OOF probabilities at the calibrated-prob ladder
  tpr <- sapply(TH, function(t) sum(calib>=t & ps$y==1)/np)
  fpr <- sapply(TH, function(t) sum(calib>=t & ps$y==0)/nn)
  auc <- abs(sum(diff(fpr)*(head(tpr,-1)+tail(tpr,-1))/2))
  list(fpr=paste(sprintf("%.4f",fpr),collapse="|"), tpr=paste(sprintf("%.4f",tpr),collapse="|"), auc=round(auc,4), n=nrow(ps))
}
proc <- function(agg, resdir, md_base, pop, is_glm){
  a <- read.csv(agg, stringsAsFactors=FALSE)
  a <- if(is_glm) a[a$featureset %in% c("fpkm","fpkm_clinics"),] else a[a$model %in% c("rf","xgb") & a$featureset %in% c("fpkm","fpkm_clinics"),]
  cat(pop, resdir, nrow(a), "configs\n")
  res <- mclapply(seq_len(nrow(a)), function(i){ r <- a[i,]
    key <- if(is_glm) sprintf("%s_%s_%s_%s_%s_glm_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$fsmethod,r$varcut) else sprintf("%s_%s_%s_%s_%s_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$model,r$varcut)
    pf <- file.path(resdir, paste0(key,"_pred_prob.csv"))
    rds <- file.path(md_base, r$set, r$cohort, r$datatype, paste0(r$featureset,"_",r$subset,".RDS"))
    v <- calib_roc(pf, rds); if(is.null(v)) return(NULL)
    data.frame(pop=pop, key=key, auc_calib=v$auc, n=v$n, fpr=v$fpr, tpr=v$tpr, stringsAsFactors=FALSE)
  }, mc.cores=8)
  do.call(rbind, res[!sapply(res,is.null)])
}
out <- rbind(
  proc("results_fsmethod_comparison.csv","res","data/model_data","full",TRUE),
  proc("rf_xgb_auc_202607.csv","res_rfxgb","data/model_data","full",FALSE),
  proc("results_fsmethod_comparison_smoking.csv","res_smoking","data/model_data_smoking","smoke",TRUE),
  proc("rf_xgb_auc_smoking.csv","res_rfxgb_smoking","data/model_data_smoking","smoke",FALSE))
write.csv(out, "roc_calib_202609.csv", row.names=FALSE)
cat("wrote roc_calib_202609.csv rows", nrow(out), "\n"); cat("CALIB_DONE\n")
