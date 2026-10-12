# Per-config operating-point / deployment reference table.
# For every model config, records:
#   - per-fold-averaged AUC (the trusted headline statistic)
#   - POOLED Youden cutoff (a deployable probability threshold + its sens/spec) -- honest only
#     for well-powered, calibrated configs; can be < 0.5 area for small-event cells.
#   - VERTICAL-averaged optimal operating point (specificity/sensitivity; no single probability)
#   - event counts + a low-event reliability flag.
# No new modelling: reads existing ROC curves, fold-averaged curves, and model_data labels.
setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
pc <- function(s) as.numeric(strsplit(as.character(s),"|",fixed=TRUE)[[1]])
# thresholds ladder used by the pooled ROC: c(Inf, 1.00..0.00, -Inf); index k=2..102 -> prob (102-k)/100... use 1-based
pooled_youden <- function(fpr,tpr){
  # finite thresholds live at indices 2..102 (prob 1.00 down to 0.00)
  idx <- 2:102; J <- tpr[idx]-fpr[idx]; b <- idx[which.max(J)]
  prob <- (102-b)/100  # index2->1.00, index102->0.00
  list(cutoff=prob, sens=tpr[b], spec=1-fpr[b], J=max(J))
}
vavg_youden <- function(tpr){ # fpr grid = (0..100)/100
  i <- 0:100; J <- tpr[i+1]-(i/100); b <- i[which.max(J)]
  list(spec=1-b/100, sens=tpr[b+1], J=max(J))
}
# event counts per (pop,set,cohort,subset) from the fpkm label RDS
ev_cache <- new.env()
events <- function(pop,set,cohort,subset){
  k <- paste(pop,set,cohort,subset)
  if(exists(k,envir=ev_cache)) return(get(k,envir=ev_cache))
  md <- if(pop=="smoke") "data/model_data_smoking" else "data/model_data"
  f <- file.path(md,set,cohort,"normalized",paste0("fpkm_",subset,".RDS"))
  v <- if(file.exists(f)){ d<-readRDS(f); c(sum(d$recurrence==1),sum(d$recurrence==0)) } else c(NA,NA)
  assign(k,v,envir=ev_cache); v
}
fa <- read.csv("roc_foldavg_202609.csv", stringsAsFactors=FALSE)
fa_key <- paste(fa$pop, fa$key, sep="|"); names(fa) <- names(fa)
rows <- list()
do_roc <- function(csv, pop, is_glm){
  r <- read.csv(csv, stringsAsFactors=FALSE)
  for(i in seq_len(nrow(r))){
    x <- r[i,]; if(is.na(x$fpr) || x$fpr=="") next
    model <- if(is_glm) paste0("glm_",x$model) else x$model
    key <- sprintf("%s_%s_%s_%s_%s_%s_v%s", x$set,x$cohort,x$datatype,x$featureset,x$subset,model,x$varcut)
    pk <- paste(pop,key,sep="|")
    py <- pooled_youden(pc(x$fpr), pc(x$tpr))
    fi <- which(fa_key==pk); v <- if(length(fi)){ vy<-vavg_youden(pc(fa$tpr[fi[1]])); list(auc=fa$auc_vavg[fi[1]], spec=vy$spec, sens=vy$sens) } else list(auc=NA,spec=NA,sens=NA)
    ec <- events(pop,x$set,x$cohort,x$subset)
    rows[[length(rows)+1]] <<- data.frame(pop=pop, set=x$set, cohort=x$cohort, datatype=x$datatype,
      featureset=x$featureset, subset=x$subset, model=model, varcut=x$varcut,
      n_recur=ec[1], n_nonrecur=ec[2], auc_perfold=round(x$mean,4),
      pooled_youden_cutoff=round(py$cutoff,2), pooled_sens=round(py$sens,3), pooled_spec=round(py$spec,3),
      vavg_auc=v$auc, vavg_opt_specificity=round(v$spec,2), vavg_opt_sensitivity=round(v$sens,3),
      low_event_flag=ifelse(is.na(ec[1]), NA, ec[1]<10 | ec[2]<10), stringsAsFactors=FALSE)
  }
}
do_roc("results_roc_curves_fsmethods.csv","full",TRUE)
do_roc("rf_xgb_roc_202607.csv","full",FALSE)
do_roc("results_roc_curves_fsmethods_smoking.csv","smoke",TRUE)
do_roc("rf_xgb_roc_smoking.csv","smoke",FALSE)
out <- do.call(rbind, rows)
write.csv(out, "model_operating_points_202609.csv", row.names=FALSE)
cat("wrote model_operating_points_202609.csv rows", nrow(out), "\n")
cat("low-event configs flagged:", sum(out$low_event_flag, na.rm=TRUE), "of", nrow(out), "\n")
cat("THRESH_DONE\n")
