# Permutation null for the nonintersection / 3vs0 headline cell (0.875).
#
# Replicates the pipeline's estimator: fold-local MAD top-N selection, cv.glmnet,
# per-fold AUC with fixed direction="<", averaged over the FINITE folds only
# (folds with 0 positives yield NA and are dropped -- exactly as tune_drug.R does).
#
# Observed and null are computed with the SAME code, so the comparison is valid
# even though this is an approximation of tune_drug.R's exact alpha grid.
suppressPackageStartupMessages({library(glmnet); library(pROC); library(caret)})

setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
m <- readRDS("data/model_data/prev/nonintersection/fpkm/fpkm_clinics_3vs0.RDS")
m <- m[!is.na(m$recurrence), ]
y <- as.numeric(as.character(m$recurrence))
X <- as.matrix(m[, setdiff(colnames(m), "recurrence"), drop = FALSE])
X <- X[, apply(X, 2, function(g) sd(g, na.rm = TRUE) > 0), drop = FALSE]
clin <- grep("^Tissue_", colnames(X))
cat("n=", nrow(X), " pos=", sum(y == 1), " neg=", sum(y == 0), " p=", ncol(X), "\n")

NTOP <- 100
fold_aucs <- function(y, nboot = 10, k = 10) {
  out <- c()
  for (b in seq_len(nboot)) {
    flds <- createFolds(factor(y), k = k, list = TRUE, returnTrain = FALSE)
    for (i in seq_along(flds)) {
      te <- flds[[i]]; tr <- setdiff(seq_len(nrow(X)), te)
      yte <- y[te]
      # fold with no positives (or no negatives) -> AUC undefined -> NA, dropped
      if (length(unique(yte)) < 2) { out <- c(out, NA_real_); next }
      gi <- setdiff(seq_len(ncol(X)), clin)
      mads <- apply(X[tr, gi, drop = FALSE], 2, mad, na.rm = TRUE)
      sel <- sort(c(clin, gi[order(mads, decreasing = TRUE)[seq_len(min(NTOP, length(gi)))]]))
      Xtr <- X[tr, sel, drop = FALSE]; Xte <- X[te, sel, drop = FALSE]
      ctr <- colMeans(Xtr); s <- apply(Xtr, 2, sd); s[s == 0 | is.na(s)] <- 1
      Xtr <- sweep(sweep(Xtr, 2, ctr, "-"), 2, s, "/")
      Xte <- sweep(sweep(Xte, 2, ctr, "-"), 2, s, "/")
      p <- tryCatch({
        cv <- cv.glmnet(Xtr, y[tr], family = "binomial", alpha = 0.5, nfolds = 5)
        as.numeric(predict(cv, Xte, s = "lambda.min", type = "response"))
      }, error = function(e) rep(0.5, length(te)))
      a <- tryCatch(as.numeric(pROC::auc(yte, p, quiet = TRUE, direction = "<")),
                    error = function(e) NA_real_)
      out <- c(out, a)
    }
  }
  out
}

set.seed(1)
obs_f <- fold_aucs(y)
obs <- mean(obs_f, na.rm = TRUE)
cat(sprintf("\nOBSERVED: mean-of-finite-folds AUC = %.4f  (finite %d / %d)\n",
            obs, sum(is.finite(obs_f)), length(obs_f)))

NPERM <- 40
nullv <- numeric(NPERM)
for (i in seq_len(NPERM)) {
  yp <- sample(y)
  nullv[i] <- mean(fold_aucs(yp), na.rm = TRUE)
  cat(sprintf("  perm %2d: %.4f\n", i, nullv[i])); flush.console()
}
p <- (1 + sum(nullv >= obs)) / (NPERM + 1)
cat(sprintf("\n=== 3vs0 nonintersection/prev/fpkm/fpkm_clinics/v100/glm_mad ===\n"))
cat(sprintf("OBSERVED=%.4f  NULLmean=%.4f  NULLsd=%.4f  95pct=%.4f  max=%.4f  p=%.4f\n",
            obs, mean(nullv), sd(nullv), quantile(nullv, 0.95), max(nullv), p))
cat("PERM3_DONE\n")
