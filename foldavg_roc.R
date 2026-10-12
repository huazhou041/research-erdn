setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(data.table); library(parallel)})
G <- seq(0,1,by=0.01)   # FPR grid (x-axis); specificity = 1 - FPR
labcache <- new.env(parent=emptyenv())
getlab <- function(rds){ if(exists(rds,envir=labcache,inherits=FALSE)) return(get(rds,envir=labcache))
  v <- if(file.exists(rds)){ m<-readRDS(rds); y<-as.numeric(as.character(m[,"recurrence"])); names(y)<-rownames(m); y } else NULL
  assign(rds,v,envir=labcache); v }
# vertical averaging: for each fold ROC, TPR interpolated at fixed FPR grid, then mean across folds
vavg <- function(pred_file, rds){
  if(!file.exists(pred_file)) return(NULL)
  y <- getlab(rds); if(is.null(y)) return(NULL)
  d <- tryCatch(fread(pred_file, showProgress=FALSE), error=function(e) NULL); if(is.null(d)||!nrow(d)) return(NULL)
  d$y <- y[as.character(d$subj)]; d <- d[is.finite(d$y) & is.finite(d$prob)]; if(!nrow(d)) return(NULL)
  grp <- split(seq_len(nrow(d)), paste(d$foldid, d$bootid))
  acc <- numeric(length(G)); cnt <- 0L
  for(ix in grp){
    yy <- d$y[ix]; pp <- d$prob[ix]; np <- sum(yy==1); nn <- sum(yy==0); if(np==0||nn==0) next
    o <- order(pp, decreasing=TRUE); ys <- yy[o]
    tp <- cumsum(ys==1); fp <- cumsum(ys==0)
    fpr <- c(0, fp/nn); tpr <- c(0, tp/np)   # ROC vertices for this fold
    # TPR at each FPR grid point (fold uses its own thresholds); linear interp, monotone
    ty <- approx(x=fpr, y=tpr, xout=G, ties="max", rule=2)$y
    acc <- acc + ty; cnt <- cnt + 1L
  }
  if(cnt==0) return(NULL)
  mt <- acc/cnt
  auc <- sum(diff(G)*(head(mt,-1)+tail(mt,-1))/2)
  list(fpr=paste(sprintf("%.4f",G),collapse="|"), tpr=paste(sprintf("%.4f",mt),collapse="|"), auc=round(auc,4), ngrp=cnt)
}
proc <- function(agg, resdir, md_base, pop, is_glm){
  a <- read.csv(agg, stringsAsFactors=FALSE)
  a <- if(is_glm) a[a$featureset %in% c("fpkm","fpkm_clinics"),] else a[a$model %in% c("rf","xgb") & a$featureset %in% c("fpkm","fpkm_clinics"),]
  cat(pop, resdir, nrow(a), "configs\n")
  res <- mclapply(seq_len(nrow(a)), function(i){ r <- a[i,]
    key <- if(is_glm) sprintf("%s_%s_%s_%s_%s_glm_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$fsmethod,r$varcut) else sprintf("%s_%s_%s_%s_%s_%s_v%s", r$set,r$cohort,r$datatype,r$featureset,r$subset,r$model,r$varcut)
    pf <- file.path(resdir, paste0(key,"_pred_prob.csv"))
    rds <- file.path(md_base, r$set, r$cohort, r$datatype, paste0(r$featureset,"_",r$subset,".RDS"))
    fa <- vavg(pf, rds); if(is.null(fa)) return(NULL)
    data.frame(pop=pop, key=key, auc_vavg=fa$auc, ngrp=fa$ngrp, fpr=fa$fpr, tpr=fa$tpr, stringsAsFactors=FALSE)
  }, mc.cores=8)
  do.call(rbind, res[!sapply(res,is.null)])
}
out <- rbind(
  proc("results_fsmethod_comparison.csv","res","data/model_data","full",TRUE),
  proc("rf_xgb_auc_202607.csv","res_rfxgb","data/model_data","full",FALSE),
  proc("results_fsmethod_comparison_smoking.csv","res_smoking","data/model_data_smoking","smoke",TRUE),
  proc("rf_xgb_auc_smoking.csv","res_rfxgb_smoking","data/model_data_smoking","smoke",FALSE))
write.csv(out, "roc_foldavg_202609.csv", row.names=FALSE)
cat("wrote roc_foldavg_202609.csv rows", nrow(out), "\n"); cat("FOLDAVG_DONE\n")
