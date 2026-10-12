#!/usr/bin/env Rscript
# MLP / deep-learning arm (h2o.deeplearning), written standalone because the DL branch
# inside tune_drug.R is dead code: line ~322 selects models via `mod$cvm[mod$lambda ==
# mod$lambda.min]`, which are cv.glmnet-only fields, so it errors on h2o S4 objects.
#
# Correctness choices that differ from the legacy script:
#  * hyperparameters are chosen by h2o INNER cross-validation on the training folds
#    (nfolds), never by scoring the held-out fold -- the legacy GLM path used
#    cv.glmnet's internal CV, and the DL path must have an equivalent or it leaks.
#  * AUC uses direction="<" (fixed orientation), same as the patched GLM script.
#  * pooled out-of-fold AUC is written alongside per-fold AUC.
suppressPackageStartupMessages({
  library(argparser); library(jsonlite); library(caret); library(pROC); library(readr); library(h2o)
})
p = arg_parser("DL (MLP) training pipeline")
p = add_argument(p, "drug", help = "outcome column")
p = add_argument(p, "matrix", help = "feature RDS")
p = add_argument(p, "hyperparams", help = "grid json")
p = add_argument(p, "prefix", help = "output prefix")
p = add_argument(p, "--varcut", "-v", default = 100, help = "top-N genes (fold-local MAD)")
p = add_argument(p, "--kfold", "-k", default = 10L, help = "outer CV folds")
p = add_argument(p, "--nboot", default = 10L, help = "repeats of k-fold CV")
p = add_argument(p, "--max_models", default = 8L, help = "hyperparameter settings sampled from grid")
p = add_argument(p, "--inner_nfolds", default = 3L, help = "h2o inner CV for model selection")
p = add_argument(p, "--max_mem", default = "16G", help = "H2O max heap")
p = add_argument(p, "--cores", default = 4L, help = "H2O threads")
argv = parse_args(p)

drug <- argv$drug; varcut <- as.integer(argv$varcut)
kfold <- as.integer(argv$kfold); kboot <- as.integer(argv$nboot)
Nmodels <- as.integer(argv$max_models); inner <- as.integer(argv$inner_nfolds)
prefix <- argv$prefix

mat <- readRDS(argv$matrix)
mat <- mat[!is.na(mat[[drug]]), , drop = FALSE]
mat[[drug]] <- as.factor(as.character(mat[[drug]]))
feats0 <- setdiff(colnames(mat), drug)
message(sprintf("%s | n=%d pos=%d features=%d varcut=%d",
                basename(argv$matrix), nrow(mat), sum(mat[[drug]] == "1"), length(feats0), varcut))

# expand grid -> list of param lists, sample Nmodels
grid <- fromJSON(argv$hyperparams, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
combos <- list()
for (blk in grid) {
  keys <- names(blk)
  idx <- expand.grid(lapply(blk, function(v) seq_along(v)), KEEP.OUT.ATTRS = FALSE)
  for (r in seq_len(nrow(idx))) {
    prm <- lapply(keys, function(k) blk[[k]][[ idx[r, k] ]])
    names(prm) <- keys
    combos[[length(combos) + 1]] <- prm
  }
}
set.seed(123)
if (length(combos) > Nmodels) combos <- combos[sample(length(combos), Nmodels)]
message(sprintf("  hyperparameter settings evaluated: %d", length(combos)))

safe_auc <- function(y, prob) {
  y <- as.numeric(as.character(y))
  if (length(unique(y[is.finite(y)])) < 2 || sd(prob, na.rm = TRUE) == 0) return(NA_real_)
  as.numeric(pROC::auc(y, as.numeric(prob), quiet = TRUE, direction = "<"))
}
fs_top <- function(m, top) {
  g <- setdiff(colnames(m), drug)
  if (length(g) <= top) return(g)
  ix <- apply(m[, g, drop = FALSE], 2, mad, na.rm = TRUE)
  ix[!is.finite(ix)] <- 0
  g[order(ix, decreasing = TRUE)[seq_len(top)]]
}

h2o.init(nthreads = as.integer(argv$cores), max_mem_size = argv$max_mem)
h2o.no_progress()

# Targeted cleanup. A blanket h2o.removeAll() races against h2o's ASYNC cross-validation
# sub-jobs and deletes models mid-flight ("Job is missing" / H2OKeyNotFoundArgumentException),
# which killed 10 of 24 jobs on the first attempt. Remove only what we own, only when done.
cleanup <- function(models, ...) {
  for (m in models) try(h2o.rm(m), silent = TRUE)
  for (fr in list(...)) try(h2o.rm(fr), silent = TRUE)
}

set.seed(123)
per_fold <- list(); oof <- list()
for (b in seq_len(kboot)) {
  flds <- createFolds(mat[[drug]], k = kfold, list = TRUE, returnTrain = FALSE)
  for (i in seq_along(flds)) {
    te <- flds[[i]]; tr <- setdiff(seq_len(nrow(mat)), te)
    sel <- fs_top(mat[tr, , drop = FALSE], varcut)
    tr_h <- as.h2o(mat[tr, c(drug, sel), drop = FALSE])
    te_h <- as.h2o(mat[te, c(drug, sel), drop = FALSE])
    # inner-CV model selection: never touches te
    best <- NULL; best_auc <- -Inf; trained <- list()
    for (prm in combos) {
      m <- tryCatch(do.call(h2o.deeplearning, c(
             list(y = drug, x = sel, training_frame = tr_h,
                  nfolds = inner, seed = 1, reproducible = FALSE,
                  keep_cross_validation_models = FALSE,
                  keep_cross_validation_predictions = FALSE), prm)),
           error = function(e) NULL)
      if (is.null(m)) next
      trained[[length(trained) + 1]] <- m
      a <- tryCatch(h2o.auc(m, xval = TRUE), error = function(e) NA_real_)
      if (is.finite(a) && a > best_auc) { best_auc <- a; best <- m }
    }
    if (is.null(best)) { cleanup(trained, tr_h, te_h); next }
    pr <- tryCatch(as.data.frame(h2o.predict(best, te_h)), error = function(e) NULL)
    if (is.null(pr) || !("p1" %in% names(pr))) { cleanup(trained, tr_h, te_h); next }
    prob <- pr$p1
    ytrue <- mat[te, drug]
    per_fold[[length(per_fold) + 1]] <- data.frame(
      foldid = i, bootid = b, inner_auc = best_auc,
      AUC = safe_auc(ytrue, prob), n_test = length(te))
    oof[[length(oof) + 1]] <- data.frame(
      subj = rownames(mat)[te], foldid = i, bootid = b, prob = prob,
      y = as.numeric(as.character(ytrue)))
    cleanup(trained, tr_h, te_h)
  }
  message(sprintf("  repeat %d/%d done", b, kboot))
}
per_fold <- do.call(rbind, per_fold); oof <- do.call(rbind, oof)
if (is.null(per_fold)) { message("no usable folds"); h2o.shutdown(prompt = FALSE); quit(save="no") }

# pooled out-of-fold AUC: one AUC per repeat over ALL subjects (the honest estimator)
pooled <- do.call(rbind, lapply(split(oof, oof$bootid), function(d)
  data.frame(bootid = d$bootid[1], pooled_AUC = safe_auc(d$y, d$prob), n = nrow(d))))
summ <- data.frame(
  matrix = basename(argv$matrix), varcut = varcut,
  n = nrow(mat), pos = sum(mat[[drug]] == "1"), n_features = length(feats0),
  meanfold_AUC = mean(per_fold$AUC, na.rm = TRUE),
  pooled_AUC_mean = mean(pooled$pooled_AUC, na.rm = TRUE),
  pooled_AUC_sd = sd(pooled$pooled_AUC, na.rm = TRUE),
  repeats = nrow(pooled))
write_csv(per_fold, paste0(prefix, "_cvSummary.csv"))
write_csv(oof,      paste0(prefix, "_pred_prob.csv"))
write_csv(summ,     paste0(prefix, "_summary.csv"))
message(sprintf("  -> meanFold AUC=%.4f   POOLED AUC=%.4f +/- %.4f",
                summ$meanfold_AUC, summ$pooled_AUC_mean, summ$pooled_AUC_sd))
h2o.shutdown(prompt = FALSE)
