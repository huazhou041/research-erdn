# Pooled ROC curves for the RF/XGB nonintersection configs, matching
# make_roc_curves_fsmethods.R exactly (same thresholds, same pipe-delimited
# fpr/tpr string format, same pooled-across-all-fold-x-boot predictions).
# Curve shape only -- the AUC columns carried through are the pipeline's own
# 100x10 fold-averaged, fixed-direction values, not re-derived here.
#
# Emits rf_xgb_roc_smoking.csv in rf_xgb_roc_202607.csv's schema
# (mean/sd/n, not mean_auc/sd_auc/n_auc).
suppressPackageStartupMessages(library(data.table))

base <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
setwd(base)

d <- read.csv("rf_xgb_auc_smoking.csv", stringsAsFactors = FALSE)
d <- d[d$model %in% c("rf", "xgb"), ]
cat("configs to draw:", nrow(d), "\n")

thresholds <- c(Inf, seq(1, 0, by = -0.01), -Inf)
fmt_vec <- function(x) paste(sprintf("%.3f", x), collapse = "|")

label_cache <- new.env(parent = emptyenv())
get_labels <- function(set, cohort, datatype, featureset, subset) {
  key <- paste(set, cohort, datatype, featureset, subset, sep = "|")
  if (exists(key, envir = label_cache, inherits = FALSE)) return(get(key, envir = label_cache, inherits = FALSE))
  rds <- file.path("data/model_data_smoking", set, cohort, datatype, paste0(featureset, "_", subset, ".RDS"))
  if (!file.exists(rds)) { assign(key, NULL, envir = label_cache); return(NULL) }
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
  if (pos == 0 || neg == 0)
    return(list(fpr = NA_character_, tpr = NA_character_, n_pred = length(y), n_pos_pred = pos, n_neg_pred = neg))
  tpr <- numeric(length(thresholds)); fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    pred <- prob >= thresholds[i]
    tpr[i] <- sum(pred & y == 1) / pos
    fpr[i] <- sum(pred & y == 0) / neg
  }
  list(fpr = fmt_vec(fpr), tpr = fmt_vec(tpr), n_pred = length(y), n_pos_pred = pos, n_neg_pred = neg)
}

rows <- vector("list", nrow(d)); missing <- character()
for (i in seq_len(nrow(d))) {
  r <- d[i, ]
  pred_file <- file.path("res_rfxgb_smoking", sprintf("%s_%s_%s_%s_%s_%s_v%s_pred_prob.csv",
                          r$set, r$cohort, r$datatype, r$featureset, r$subset, r$model, r$varcut))
  out <- as.list(r[, c("set","cohort","datatype","featureset","subset","model","varcut","mean","sd","n")])
  na_row <- list(fpr = NA_character_, tpr = NA_character_, n_pred = NA_integer_,
                 n_pos_pred = NA_integer_, n_neg_pred = NA_integer_)
  if (!file.exists(pred_file)) { missing <- c(missing, pred_file); rows[[i]] <- c(out, na_row); next }
  pred <- tryCatch(fread(pred_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(pred)) { missing <- c(missing, pred_file); rows[[i]] <- c(out, na_row); next }
  y_map <- get_labels(r$set, r$cohort, r$datatype, r$featureset, r$subset)
  if (is.null(y_map)) { missing <- c(missing, pred_file); rows[[i]] <- c(out, na_row); next }
  y <- unname(y_map[as.character(pred$subj)])
  rows[[i]] <- c(out, calc_curve(y, as.numeric(pred$prob)))
  if (i %% 200 == 0) cat("processed", i, "of", nrow(d), "(", length(missing), "missing )\n")
}

roc <- do.call(rbind.data.frame, rows)
write.csv(roc, "rf_xgb_roc_smoking.csv", row.names = FALSE)
cat("wrote rf_xgb_roc_smoking.csv rows", nrow(roc), " missing", length(missing), "\n")
cat("ROC_RFXGB_DONE\n")
