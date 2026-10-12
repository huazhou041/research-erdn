#!/usr/bin/env Rscript
# Random Forest training pipeline for recurrence RDS matrices.
# Mirrors the harness of tune_drug_xgb.R (bootstrap x k-fold CV, top-MAD feature
# selection, pred_prob + cvSummary outputs, final whole-data model) but uses the
# proven randomForest implementation from val_20260605/tune_drug_extended.R
# (fixed ntree, mtry = floor(sqrt(#features))). No h2o/Java required.

suppressPackageStartupMessages({
  library(argparser)
  library(caret)
  library(pROC)
  library(randomForest)
  library(readr)
})

p <- arg_parser("Run Random Forest training pipeline for recurrence RDS matrices")
p <- add_argument(p, "drug", help = "Outcome column, e.g. recurrence")
p <- add_argument(p, "matrix", help = "RDS model matrix")
p <- add_argument(p, "hyperparams", help = "JSON grid (unused; kept for CLI parity)")
p <- add_argument(p, "prefix", help = "Output prefix")
p <- add_argument(p, "--varcut", default = 100L, help = "Number of top-MAD gene features to keep")
p <- add_argument(p, "--nboot", default = 50L, help = "Number of repeated CV runs")
p <- add_argument(p, "--kfold", default = 5L, help = "Number of CV folds")
p <- add_argument(p, "--ntree", default = 500L, help = "Number of trees")
p <- add_argument(p, "--test", help = "Optional test RDS")
argv <- parse_args(p)

drug <- argv$drug
mat_file <- argv$matrix
prefix <- argv$prefix
varcut <- as.integer(argv$varcut)
kboot <- as.integer(argv$nboot)
kfold <- as.integer(argv$kfold)
ntree <- as.integer(argv$ntree)
test_file <- argv$test

cat("*****************\n")
cat("RF :", mat_file, "\n")
cat("*****************\n")

select_features <- function(mat, genes, tissues, top) {
  if (length(genes) == 0) {
    return(tissues)
  }
  top <- min(top, length(genes))
  ix <- apply(mat[, genes, drop = FALSE], 2, mad, na.rm = TRUE)
  selected_genes <- names(ix)[order(ix, decreasing = TRUE)[seq_len(top)]]
  c(tissues, selected_genes)
}

safe_auc <- function(y, prob) {
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  # direction="<": predictor (prob) assumed LOWER for recurrence=0 than recurrence=1;
  # see tune_drug.R for why "auto" (the default) is unsafe.
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

train_rf <- function(train_df, x_cols) {
  # x/y interface (not formula) so non-syntactic gene names are handled safely
  y <- factor(train_df[, drug], levels = c(0, 1))
  x <- train_df[, x_cols, drop = FALSE]
  randomForest(
    x = x,
    y = y,
    ntree = ntree,
    mtry = max(1L, floor(sqrt(length(x_cols)))),
    importance = FALSE
  )
}

predict_rf <- function(model, df, x_cols) {
  # P(class == 1)
  predict(model, df[, x_cols, drop = FALSE], type = "prob")[, "1"]
}

mat <- readRDS(mat_file)
mat <- mat[!is.na(mat[, drug]), , drop = FALSE]
mat[, drug] <- as.numeric(as.character(mat[, drug]))

tissues <- grep("^Tissue_", colnames(mat), value = TRUE)
genes <- setdiff(colnames(mat), c(drug, tissues))
origmat <- mat

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

    model <- train_rf(train_df, x_cols)
    prob <- predict_rf(model, valid_df, x_cols)
    y <- as.numeric(valid_df[, drug])

    cv_rows[[i]] <- data.frame(
      fold_id = i,
      bootid = tt,
      model_id = 1L,
      ntree = ntree,
      mtry = max(1L, floor(sqrt(length(x_cols)))),
      MSE = mean((prob - y) ^ 2),
      AUC = safe_auc(y, prob),
      stringsAsFactors = FALSE
    )
    pred_rows[[i]] <- data.frame(
      subj = rownames(valid_df),
      foldid = i,
      bootid = tt,
      prob = prob,
      stringsAsFactors = FALSE
    )
  }

  cvtotal[[tt]] <- do.call(rbind, cv_rows)
  predtotal[[tt]] <- do.call(rbind, pred_rows)
}

cvtotal <- do.call(rbind, cvtotal)
predtotal <- do.call(rbind, predtotal)
readr::write_csv(cvtotal, paste0(prefix, "_cvSummary.csv"))
readr::write_csv(predtotal, paste0(prefix, "_pred_prob.csv"))

# Final model on all data, whole-data feature set
x_cols <- select_features(origmat, genes, tissues, varcut)
final_model <- train_rf(origmat[, c(drug, x_cols), drop = FALSE], x_cols)
saveRDS(list(model = final_model, features = x_cols,
             params = list(ntree = ntree, mtry = max(1L, floor(sqrt(length(x_cols)))))),
        paste0(prefix, ".rds"))

if (!is.na(test_file)) {
  test <- readRDS(test_file)
  patient_id <- rownames(test)
  prob <- predict_rf(final_model, test, x_cols)
  readr::write_csv(data.frame(patient_id = patient_id, recurrence_prob = prob),
                   paste0(prefix, "_clinicalPrediction.csv"))
}

cat("done\n", date(), "\n")
