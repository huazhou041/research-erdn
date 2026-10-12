SRC <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_20260709"
NEW <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
cmp <- function(fname) {
  a <- read.csv(file.path(SRC, fname), stringsAsFactors=FALSE)
  b <- read.csv(file.path(NEW, fname), stringsAsFactors=FALSE)
  K <- c("set","cohort","datatype","featureset","subset","model","varcut")
  ak <- do.call(paste, c(a[K], sep="|")); bk <- do.call(paste, c(b[K], sep="|"))
  common <- intersect(ak,bk)
  ai <- a[match(common,ak),]; bi <- b[match(common,bk),]
  d <- bi$mean - ai$mean; d <- d[is.finite(d)]
  cat("\n=====", fname, "=====\n")
  cat(sprintf("prev=%d new=%d common=%d only-prev=%d only-new=%d\n",
              nrow(a),nrow(b),length(common),length(setdiff(ak,bk)),length(setdiff(bk,ak))))
  cat(sprintf("delta mean: mean=%.4f meanabs=%.4f maxabs=%.4f cor=%.4f identical<0.001=%.1f%%\n",
              mean(d),mean(abs(d)),max(abs(d)),cor(ai$mean,bi$mean,use="complete.obs"),100*mean(abs(d)<0.001)))
  op <- a[!ak %in% bk,]
  if (nrow(op)) { cat("--- only-in-PREV by cohort ---\n"); print(table(op$cohort));
                  cat("--- only-in-PREV by model ---\n"); print(table(op$model)) }
}
cmp("rf_xgb_auc_202607.csv")
cmp("rf_xgb_auc_smoking.csv")
