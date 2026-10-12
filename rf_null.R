suppressPackageStartupMessages({library(randomForest); library(pROC); library(caret)})
set.seed(1)
pooled_rf_auc <- function(X, y, ntop=200, k=10) {
  n <- nrow(X); flds <- createFolds(factor(y), k=k, list=TRUE, returnTrain=FALSE)
  prob <- rep(NA_real_, n)
  for (i in seq_along(flds)) {
    te <- flds[[i]]; tr <- setdiff(seq_len(n), te)
    mads <- apply(X[tr,,drop=FALSE], 2, mad, na.rm=TRUE)
    sel <- order(mads, decreasing=TRUE)[seq_len(min(ntop, ncol(X)))]
    rf <- randomForest(x=X[tr,sel,drop=FALSE], y=factor(y[tr]), ntree=500)
    prob[te] <- predict(rf, X[te,sel,drop=FALSE], type="prob")[,"1"]
  }
  as.numeric(pROC::auc(y, prob, quiet=TRUE, direction="<"))
}
for (s in c("prev","current")) {
  m <- readRDS(sprintf("data/model_data/%s/nonintersection/normalized/fpkm_123vs0.RDS", s))
  m <- m[!is.na(m$recurrence),]; y <- as.numeric(as.character(m$recurrence))
  X <- as.matrix(m[, setdiff(colnames(m),"recurrence")]); X <- X[, apply(X,2,function(g) sd(g,na.rm=TRUE)>0)]
  gm <- apply(X,2,mad,na.rm=TRUE); X <- X[, order(gm,decreasing=TRUE)[1:3000]]  # superset of fold top-200
  obs <- mean(replicate(3, pooled_rf_auc(X,y)))
  nperm <- 30; nullv <- replicate(nperm, pooled_rf_auc(X, sample(y)))
  cat(sprintf("%s: OBSERVED=%.3f  NULL mean=%.3f 95pct=%.3f max=%.3f  p=%.3f\n",
      s, obs, mean(nullv), quantile(nullv,0.95), max(nullv), (1+sum(nullv>=obs))/(nperm+1)))
}
