#!/bin/bash
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript -e '
suppressPackageStartupMessages(library(data.table))
base <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
setwd(base)
d <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors = FALSE)
stopifnot(all(d$model == "glm"))
names(d)[names(d) == "fsmethod"] <- "fs_method"
d <- d[1:20, ]  # smoke test subset

thresholds <- c(Inf, seq(1, 0, by = -0.01), -Inf)
fmt_vec <- function(x) paste(sprintf("%.3f", x), collapse = "|")
label_cache <- new.env(parent = emptyenv())
get_labels <- function(set, cohort, datatype, featureset, subset) {
  key <- paste(set, cohort, datatype, featureset, subset, sep = "|")
  if (exists(key, envir = label_cache, inherits = FALSE)) return(get(key, envir = label_cache, inherits = FALSE))
  rds <- file.path("data/model_data", set, cohort, datatype, paste0(featureset, "_", subset, ".RDS"))
  mat <- readRDS(rds)
  y <- as.numeric(as.character(mat[, "recurrence"]))
  names(y) <- rownames(mat)
  assign(key, y, envir = label_cache)
  y
}
calc_curve <- function(y, prob) {
  ok <- !is.na(y) & !is.na(prob)
  y <- y[ok]; prob <- prob[ok]
  pos <- sum(y == 1); neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(list(fpr=NA_character_,tpr=NA_character_,n_pred=length(y),n_pos_pred=pos,n_neg_pred=neg))
  tpr <- numeric(length(thresholds)); fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    pred <- prob >= thresholds[i]
    tpr[i] <- sum(pred & y == 1) / pos
    fpr[i] <- sum(pred & y == 0) / neg
  }
  list(fpr=fmt_vec(fpr), tpr=fmt_vec(tpr), n_pred=length(y), n_pos_pred=pos, n_neg_pred=neg)
}
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  pred_file <- file.path("res", sprintf("%s_%s_%s_%s_%s_glm_%s_v%s_pred_prob.csv",
                                          r$set, r$cohort, r$datatype, r$featureset, r$subset, r$fs_method, r$varcut))
  if (!file.exists(pred_file)) { cat("MISSING:", pred_file, "\n"); next }
  pred <- fread(pred_file, showProgress = FALSE)
  y_map <- get_labels(r$set, r$cohort, r$datatype, r$featureset, r$subset)
  y <- unname(y_map[as.character(pred$subj)])
  cat(sprintf("row %d: %s | n=%d, n_NA_label=%d, pos=%d neg=%d, reported_mean_auc=%.4f\n",
              i, basename(pred_file), length(y), sum(is.na(y)), sum(y==1,na.rm=T), sum(y==0,na.rm=T), r$mean_auc))
  curve <- calc_curve(y, as.numeric(pred$prob))
  cat(sprintf("   pooled curve: n_pred=%d pos=%d neg=%d\n", curve$n_pred, curve$n_pos_pred, curve$n_neg_pred))
}
'
