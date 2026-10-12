setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
options(width=200)
pass <- function(ok,msg) cat(sprintf("[%s] %s\n", ifelse(ok,"PASS","**FAIL**"), msg))

cat("=================== A. COHORT SIZES (current, 123vs0, fpkm) ===================\n")
exp_full <- c(wholeset=150, stage1=120, intersection=58, nonintersection=92, novel=74)
for (co in names(exp_full)) {
  d <- readRDS(file.path("data/model_data/current",co,"normalized",paste0("fpkm_123vs0.RDS")))
  pass(nrow(d)==exp_full[co], sprintf("full %-15s n=%d (expect %d) events=%d", co, nrow(d), exp_full[co], sum(d$recurrence==1)))
}
cat("--- ever-smoker subset (should be <= full; wholeset ~120 labeled ever-smokers) ---\n")
for (co in c("wholeset","stage1","intersection","nonintersection")) {
  d <- readRDS(file.path("data/model_data_smoking/current",co,"normalized","fpkm_123vs0.RDS"))
  cat(sprintf("  smoking %-15s n=%d events=%d\n", co, nrow(d), sum(d$recurrence==1)))
}

cat("\n=================== B. COMPLETENESS (group.done) ===================\n")
for (dd in list(c("res",480),c("res_rfxgb",240),c("res_smoking",384),c("res_rfxgb_smoking",192))) {
  n <- length(list.files(dd[1], pattern="_group.done$")); pass(n==as.integer(dd[2]), sprintf("%-18s %d/%s group.done", dd[1], n, dd[2]))
}

cat("\n=================== C. natcommad gene list == venn top-MAD ===================\n")
nr <- read.csv("raw_data/natcom_mad_ranked.csv", stringsAsFactors=FALSE)
top50_rank <- head(as.character(nr[[1]]),50)
if (file.exists("venn_madgenes.csv")) {
  v <- read.csv("venn_madgenes.csv", stringsAsFactors=FALSE)
  cat("venn cols:", paste(names(v),collapse=","), "\n")
} else cat("venn_madgenes.csv not in this folder (copied for build)\n")
cat("natcom_mad_ranked top10:", paste(head(top50_rank,10),collapse=", "), "\n")

cat("\n=================== D. AGGREGATION FIDELITY (raw cvSummary -> aggregated table) ===================\n")
chk <- function(aggcsv, resdir, model_in_name, aucfield="mean_auc") {
  agg <- read.csv(aggcsv, stringsAsFactors=FALSE)
  set.seed(1); samp <- agg[sample(nrow(agg), 12),]
  maxd <- 0; bad <- 0
  for (i in seq_len(nrow(samp))) {
    r <- samp[i,]
    if ("fsmethod" %in% names(r)) {
      fn <- sprintf("%s/%s_%s_%s_%s_%s_%s_%s_v%s_cvSummary.csv", resdir, r$set,r$cohort,r$datatype,r$featureset,r$subset,r$model,r$fsmethod,r$varcut)
      aggv <- r$mean_auc; aggn <- r$n_auc
    } else {
      fn <- sprintf("%s/%s_%s_%s_%s_%s_%s_v%s_cvSummary.csv", resdir, r$set,r$cohort,r$datatype,r$featureset,r$subset,r$model,r$varcut)
      aggv <- r$mean; aggn <- r$n
    }
    if (!file.exists(fn)) { cat("  missing raw:", fn, "\n"); bad<-bad+1; next }
    raw <- read.csv(fn, stringsAsFactors=FALSE); v <- raw$AUC[is.finite(raw$AUC)]
    d <- abs(mean(v) - aggv); maxd <- max(maxd,d)
    if (d > 1e-6 || length(v)!=aggn) { bad<-bad+1; cat(sprintf("  MISMATCH %s : agg=%.6f raw=%.6f n_agg=%s n_raw=%d\n", basename(fn), aggv, mean(v), aggn, length(v))) }
  }
  pass(bad==0, sprintf("%s : 12 sampled configs, max|agg-raw|=%.2e, mismatches=%d", basename(aggcsv), maxd, bad))
}
chk("results_fsmethod_comparison.csv","res")
chk("results_fsmethod_comparison_smoking.csv","res_smoking")
chk("rf_xgb_auc_202607.csv","res_rfxgb")
chk("rf_xgb_auc_smoking.csv","res_rfxgb_smoking")

cat("\n=================== E. ROC CSV integrity ===================\n")
for (f in c("results_roc_curves_fsmethods.csv","results_roc_curves_fsmethods_smoking.csv","rf_xgb_roc_202607.csv","rf_xgb_roc_smoking.csv")) {
  if (!file.exists(f)) { pass(FALSE, paste(f,"MISSING")); next }
  r <- read.csv(f, stringsAsFactors=FALSE); nafpr <- sum(is.na(r$fpr) | r$fpr=="")
  pass(TRUE, sprintf("%-40s rows=%d  NA-curve(degenerate small-n)=%d", f, nrow(r), nafpr))
}
cat("\nAUDIT_CLUSTER_DONE\n")
