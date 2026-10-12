suppressPackageStartupMessages({
  library(randomForest); library(pROC); library(caret); library(glmnet); library(xgboost)
})

# Pooled out-of-fold AUC, fixed direction="<", fold-local unsupervised MAD gene
# selection. Clinical covariates (the `Tissue_*` columns) are ALWAYS carried into
# the model -- they must never compete with ~22k genes for a top-200 MAD slot,
# which would silently drop them and make `fpkm_clinics` identical to `fpkm`.
CLIN_RE <- "^Tissue_"

select_local <- function(Xtr, ntop, clin_idx) {
  gene_idx <- setdiff(seq_len(ncol(Xtr)), clin_idx)
  if (!length(gene_idx)) return(clin_idx)
  if (length(gene_idx) <= ntop) return(sort(c(clin_idx, gene_idx)))
  mads <- apply(Xtr[, gene_idx, drop = FALSE], 2, mad, na.rm = TRUE)
  sort(c(clin_idx, gene_idx[order(mads, decreasing = TRUE)[seq_len(ntop)]]))
}

pooled_auc <- function(X, y, model, clin_idx, ntop = 200, k = 10, seed = 1) {
  set.seed(seed)
  flds <- createFolds(factor(y), k = k, list = TRUE, returnTrain = FALSE)
  prob <- rep(NA_real_, nrow(X))
  for (i in seq_along(flds)) {
    te <- flds[[i]]; tr <- setdiff(seq_len(nrow(X)), te)
    sel <- select_local(X[tr, , drop = FALSE], ntop, clin_idx)
    Xtr <- X[tr, sel, drop = FALSE]; Xte <- X[te, sel, drop = FALSE]; ytr <- y[tr]
    if (model == "glm") {
      ctr <- colMeans(Xtr); scl <- apply(Xtr, 2, sd); scl[scl == 0 | is.na(scl)] <- 1
      Xtr2 <- sweep(sweep(Xtr, 2, ctr, "-"), 2, scl, "/")
      Xte2 <- sweep(sweep(Xte, 2, ctr, "-"), 2, scl, "/")
      if (ncol(Xtr2) < 2) {
        fit <- glm(y ~ x, data.frame(y = factor(ytr), x = Xtr2[, 1]), family = binomial)
        prob[te] <- predict(fit, data.frame(x = Xte2[, 1]), type = "response")
      } else {
        cv <- cv.glmnet(Xtr2, ytr, family = "binomial", alpha = 0.5, nfolds = 5)
        prob[te] <- as.numeric(predict(cv, Xte2, s = "lambda.min", type = "response"))
      }
    } else if (model == "rf") {
      rf <- randomForest(x = Xtr, y = factor(ytr), ntree = 400)
      prob[te] <- predict(rf, Xte, type = "prob")[, "1"]
    } else if (model == "xgb") {
      bst <- xgb.train(params = list(max_depth = 3, eta = 0.1, subsample = 0.8,
                       colsample_bytree = 0.8, objective = "binary:logistic"),
                       data = xgb.DMatrix(as.matrix(Xtr), label = ytr),
                       nrounds = 100, verbose = 0)
      prob[te] <- predict(bst, xgb.DMatrix(as.matrix(Xte)))
    }
  }
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

co <- "nonintersection"; dt <- "normalized"
feats <- c("fpkm", "fpkm_clinics", "module20", "module20_clinics", "clinics")
tasks <- c("123vs0", "12vs0", "3vs0")
models <- c("glm", "rf", "xgb")

rows <- list()
for (s in c("prev", "current")) {
  for (ft in feats) for (tk in tasks) {
    f <- sprintf("data/model_data/%s/%s/%s/%s_%s.RDS", s, co, dt, ft, tk)
    if (!file.exists(f)) next
    m <- readRDS(f); m <- m[!is.na(m$recurrence), ]
    y <- as.numeric(as.character(m$recurrence))
    X <- as.matrix(m[, setdiff(colnames(m), "recurrence"), drop = FALSE])
    X <- X[, apply(X, 2, function(g) sd(g, na.rm = TRUE) > 0), drop = FALSE]
    if (length(unique(y)) < 2 || min(table(y)) < 3) next
    clin_idx <- grep(CLIN_RE, colnames(X))
    a <- vapply(models, function(mdl)
      tryCatch(mean(vapply(1:3, function(sd) pooled_auc(X, y, mdl, clin_idx, seed = sd), numeric(1))),
               error = function(e) NA_real_), numeric(1))
    rows[[paste(s, ft, tk)]] <- data.frame(
      quant = s, cohort = co, datatype = dt, feature = ft, task = tk,
      n = nrow(X), events = sum(y), p = ncol(X), n_clin = length(clin_idx),
      GLM = round(a[["glm"]], 4), RF = round(a[["rf"]], 4), XGB = round(a[["xgb"]], 4),
      stringsAsFactors = FALSE)
    cat(sprintf("%s %s %s: n=%d ev=%d p=%d clin=%d | GLM=%.3f RF=%.3f XGB=%.3f\n",
        s, ft, tk, nrow(X), sum(y), ncol(X), length(clin_idx),
        a[["glm"]], a[["rf"]], a[["xgb"]]))
    flush.console()
  }
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, "auc_table_nonint_92.csv", row.names = FALSE)
cat("TABLE_DONE rows=", nrow(tab), "\n")
