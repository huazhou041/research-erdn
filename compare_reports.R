SRC <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_20260709"
NEW <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
cmp <- function(fname, keycols) {
  a <- read.csv(file.path(SRC, fname), stringsAsFactors=FALSE)
  b <- read.csv(file.path(NEW, fname), stringsAsFactors=FALSE)
  cat("\n=====", fname, "=====\n")
  cat(sprintf("prev rows=%d  new rows=%d\n", nrow(a), nrow(b)))
  ak <- do.call(paste, c(a[keycols], sep="|")); bk <- do.call(paste, c(b[keycols], sep="|"))
  common <- intersect(ak, bk)
  cat(sprintf("common configs=%d | only-in-prev=%d | only-in-new=%d\n",
              length(common), length(setdiff(ak,bk)), length(setdiff(bk,ak))))
  ai <- a[match(common, ak),]; bi <- b[match(common, bk),]
  d <- bi$mean_auc - ai$mean_auc
  d <- d[is.finite(d)]
  cat(sprintf("mean_auc delta (new-prev): mean=%.4f  meanabs=%.4f  max abs=%.4f  cor=%.4f\n",
              mean(d), mean(abs(d)), max(abs(d)), suppressWarnings(cor(ai$mean_auc, bi$mean_auc, use="complete.obs"))))
  cat(sprintf("identical(<0.001)=%.1f%%  within0.02=%.1f%%\n",
              100*mean(abs(d)<0.001), 100*mean(abs(d)<0.02)))
  # what's only in new (by fsmethod/model)
  onlynew <- b[!bk %in% ak,]
  if (nrow(onlynew)) { cat("only-in-NEW breakdown:\n"); print(table(onlynew$fsmethod, onlynew$model)) }
  onlyprev <- a[!ak %in% bk,]
  if (nrow(onlyprev)) { cat("only-in-PREV breakdown (featureset x varcut sample):\n"); print(head(table(onlyprev$featureset, onlyprev$varcut))) }
}
K <- c("set","cohort","datatype","featureset","subset","model","fsmethod","varcut")
cmp("results_fsmethod_comparison.csv", K)
cmp("results_fsmethod_comparison_smoking.csv", K)
# rfxgb: check columns first
rh <- names(read.csv(file.path(NEW,"rf_xgb_auc_202607.csv"), nrows=1))
cat("\nrfxgb columns:", paste(rh, collapse=","), "\n")
Krf <- intersect(c("set","cohort","datatype","featureset","subset","model","varcut"), rh)
cmp("rf_xgb_auc_202607.csv", Krf)
cmp("rf_xgb_auc_smoking.csv", Krf)
