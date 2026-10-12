#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(caret)
  library(jsonlite)
  library(parallel)
  library(pROC)
  library(readr)
  library(xgboost)
})

p <- arg_parser("Run XGBoost training pipeline for recurrence RDS matrices")
p <- add_argument(p, "drug", help = "Outcome column, e.g. recurrence")
p <- add_argument(p, "matrix", help = "RDS model matrix")
p <- add_argument(p, "hyperparams", help = "JSON hyperparameter grid")
p <- add_argument(p, "prefix", help = "Output prefix")
p <- add_argument(p, "--max_models", default = 25L, help = "Max random grid models per fold")
p <- add_argument(p, "--varcut", default = 100L, help = "Number of top-MAD gene features to keep")
p <- add_argument(p, "--nboot", default = 50L, help = "Number of repeated CV runs")
p <- add_argument(p, "--kfold", default = 5L, help = "Number of CV folds")
p <- add_argument(p, "--nthread", default = 4L, help = "XGBoost threads")
p <- add_argument(p, "--inner_kfold", default = 10L, help = "Inner CV folds for nested hyperparameter selection (10 to mirror cv.glmnet)")
p <- add_argument(p, "--test", help = "Optional test RDS")
argv <- parse_args(p)

drug <- argv$drug
mat_file <- argv$matrix
hyperparam_file <- argv$hyperparams
prefix <- argv$prefix
Nmodels <- as.integer(argv$max_models)
varcut <- as.integer(argv$varcut)
kboot <- as.integer(argv$nboot)
kfold <- as.integer(argv$kfold)
nthread <- as.integer(argv$nthread)
inner_kfold <- as.integer(argv$inner_kfold)
test_file <- argv$test

cat("*****************\n")
cat("XGB :", mat_file, "\n")
cat("*****************\n")

rbind_data_frame <- function(...) {
  rbind.data.frame(..., make.row.names = FALSE, stringsAsFactors = FALSE)
}

expand_grid <- function(hyper_params) {
  if (is.data.frame(hyper_params)) {
    grids <- lapply(seq_len(nrow(hyper_params)), function(i) {
      row <- hyper_params[i, , drop = FALSE]
      values <- lapply(row, function(x) {
        if (is.list(x)) x[[1]] else x
      })
      do.call(expand.grid, c(values, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
    })
  } else {
    grids <- lapply(hyper_params, function(x) {
      do.call(expand.grid, c(x, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
    })
  }
  grid <- Reduce(rbind_data_frame, grids)
  rownames(grid) <- NULL
  grid
}

select_features <- function(mat, genes, tissues, top) {
  if (length(genes) == 0) {
    return(tissues)
  }
  top <- min(top, length(genes))
  ix <- apply(mat[, genes, drop = FALSE], 2, mad, na.rm = TRUE)
  selected_genes <- names(ix)[order(ix, decreasing = TRUE)[seq_len(top)]]
  c(tissues, selected_genes)
}

as_xgb_matrix <- function(df, x_cols) {
  x <- data.matrix(df[, x_cols, drop = FALSE])
  storage.mode(x) <- "double"
  # Drop feature names: xgboost 3.x predict() serialises them via jsonlite on every
  # call, which corrupts memory ("recursive gc invocation" -> segfault) under the
  # volume of the nested CV loop. Columns are matched by position (train & predict
  # use the same x_cols in the same order), so names are unnecessary.
  dimnames(x) <- NULL
  x
}

safe_auc <- function(y, prob) {
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  # direction="<": predictor (prob) assumed LOWER for recurrence=0 than recurrence=1;
  # see tune_drug.R for why "auto" (the default) is unsafe.
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

# Train from a pre-built DMatrix so one DMatrix can be reused across the whole
# grid instead of rebuilding one per model -- far fewer external-pointer objects
# for R's GC to churn, which was causing recursive-gc segfaults in the 10-fold
# nested loop on large matrices.
train_xgb_on <- function(dtrain, params) {
  xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      base_score = 0.5,   # fixed init; avoids crash when an inner fold is single-class
      eta = as.numeric(params$eta),
      max_depth = as.integer(params$max_depth),
      min_child_weight = as.numeric(params$min_child_weight),
      subsample = as.numeric(params$subsample),
      colsample_bytree = as.numeric(params$colsample_bytree),
      lambda = as.numeric(params$lambda),
      alpha = as.numeric(params$alpha),
      nthread = nthread
    ),
    data = dtrain,
    nrounds = as.integer(params$nrounds),
    verbose = 0
  )
}

train_xgb <- function(train_df, x_cols, params) {
  dtrain <- xgb.DMatrix(data = as_xgb_matrix(train_df, x_cols),
                        label = as.numeric(train_df[, drug]))
  train_xgb_on(dtrain, params)
}

hyper_params <- jsonlite::fromJSON(hyperparam_file, simplifyVector = TRUE)
hyper_grid <- expand_grid(hyper_params)
Nmodels <- min(Nmodels, nrow(hyper_grid))
set.seed(123)
hyper_grid <- hyper_grid[sample(nrow(hyper_grid), Nmodels), , drop = FALSE]
hyper_grid$model_id <- seq_len(nrow(hyper_grid))
hyper_grid <- hyper_grid[, c("model_id", setdiff(names(hyper_grid), "model_id")), drop = FALSE]

mat <- readRDS(mat_file)
mat <- mat[!is.na(mat[, drug]), , drop = FALSE]
mat[, drug] <- as.numeric(as.character(mat[, drug]))

tissues <- grep("^Tissue_", colnames(mat), value = TRUE)
genes <- setdiff(colnames(mat), c(drug, tissues))
origmat <- mat

logloss <- function(y, p) {
  eps <- 1e-15; p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# Nested hyperparameter selection to mirror GLM's cv.glmnet as closely as possible:
#   * INNER CV on the training data only (inner_kfold folds; 10 by default, like
#     cv.glmnet), so the config is never tuned on the fold it is scored on.
#   * config chosen by minimum mean binomial DEVIANCE (2*logloss), i.e. cv.glmnet's
#     `cvm` criterion -- always defined (unlike AUC on single-class inner folds).
# The outer-fold performance is still recorded/reported as AUC.
# (RF, by contrast, uses a single fixed config; GLM tunes alpha/lambda by cvm.)
select_best_params <- function(train_df, x_cols, inner_k) {
  y <- as.numeric(train_df[, drug])
  ik <- max(2L, min(inner_k, sum(y == 1), sum(y == 0)))
  inner <- caret::createFolds(y, k = ik, list = TRUE, returnTrain = FALSE)
  dev_sum <- numeric(nrow(hyper_grid)); dev_cnt <- numeric(nrow(hyper_grid))
  for (ii in seq_along(inner)) {
    itr <- train_df[-inner[[ii]], , drop = FALSE]
    iva <- train_df[ inner[[ii]], , drop = FALSE]
    iy  <- as.numeric(iva[, drug])
    dtr <- xgb.DMatrix(data = as_xgb_matrix(itr, x_cols), label = as.numeric(itr[, drug]))
    xva <- as_xgb_matrix(iva, x_cols)
    for (j in seq_len(nrow(hyper_grid))) {
      m  <- train_xgb_on(dtr, as.list(hyper_grid[j, , drop = FALSE]))
      pr <- predict(m, xva)
      dev_sum[j] <- dev_sum[j] + 2 * logloss(iy, pr)   # binomial deviance
      dev_cnt[j] <- dev_cnt[j] + 1
      rm(m)
    }
    rm(dtr, xva); gc(FALSE)   # free per-inner-fold objects; time GC outside xgb.train
  }
  which.min(dev_sum / dev_cnt)
}

cvtotal <- vector("list", kboot)
predtotal <- vector("list", kboot)

set.seed(123)
for (tt in seq_len(kboot)) {
  mat <- origmat
  flds <- caret::createFolds(mat[, drug], k = kfold, list = TRUE, returnTrain = FALSE)
  fold_id <- integer(nrow(mat))
  for (i in seq_len(kfold)) {
    fold_id[flds[[i]]] <- i
  }

  cv_rows <- list()
  pred_rows <- list()
  for (i in seq_len(kfold)) {
    x_cols <- select_features(mat[fold_id != i, , drop = FALSE], genes, tissues, varcut)
    train_df <- mat[fold_id != i, c(drug, x_cols), drop = FALSE]
    valid_df <- mat[fold_id == i, c(drug, x_cols), drop = FALSE]

    # pick hyperparameter on training data only (inner CV) ...
    best_j <- select_best_params(train_df, x_cols, inner_kfold)
    best_params <- as.list(hyper_grid[best_j, , drop = FALSE])
    # ... then train it on the full training fold and score ONCE on the outer fold
    model <- train_xgb(train_df, x_cols, best_params)
    prob  <- predict(model, as_xgb_matrix(valid_df, x_cols))
    y     <- as.numeric(valid_df[, drug])
    cv_rows[[i]] <- data.frame(
      fold_id = i, bootid = tt,
      hyper_grid[best_j, , drop = FALSE],
      MSE = mean((prob - y) ^ 2), AUC = safe_auc(y, prob),
      stringsAsFactors = FALSE
    )
    pred_rows[[i]] <- data.frame(
      subj = rownames(valid_df), foldid = i, bootid = tt, prob = prob,
      stringsAsFactors = FALSE
    )
  }

  cvtotal[[tt]] <- do.call(rbind, cv_rows)
  predtotal[[tt]] <- do.call(rbind, pred_rows)
}

cvtotal <- do.call(rbind, cvtotal)          # one row per (fold, boot) = kfold*kboot
predtotal <- do.call(rbind, predtotal)
readr::write_csv(cvtotal, paste0(prefix, "_cvSummary.csv"))
readr::write_csv(predtotal, paste0(prefix, "_pred_prob.csv"))

# Final deployed model: select config by inner CV on all data, train on all data.
x_cols <- select_features(origmat, genes, tissues, varcut)
final_df <- origmat[, c(drug, x_cols), drop = FALSE]
best_j <- select_best_params(final_df, x_cols, inner_kfold)
best_params <- as.list(hyper_grid[best_j, , drop = FALSE])
final_model <- train_xgb(final_df, x_cols, best_params)
saveRDS(list(model = final_model, features = x_cols, params = best_params), paste0(prefix, ".rds"))

if (!is.na(test_file)) {
  test <- readRDS(test_file)
  patient_id <- rownames(test)
  prob <- predict(final_model, as_xgb_matrix(test, x_cols))
  readr::write_csv(data.frame(patient_id = patient_id, recurrence_prob = prob), paste0(prefix, "_clinicalPrediction.csv"))
}

cat("done\n", date(), "\n")
