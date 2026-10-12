# Pooled ROC curves for the GLM *nonintersection* configs only.
#
# Identical method to make_roc_curves_fsmethods.R (same thresholds, same
# pipe-delimited fpr/tpr format, same pooled-across-all-fold-x-boot predictions,
# curve shape only). Scoped to nonintersection so the existing 5400 curve rows in
# results_roc_curves_fsmethods.csv are left byte-identical and merely appended to,
# rather than regenerated.
suppressPackageStartupMessages(library(data.table))

base <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
setwd(base)

d <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors = FALSE)
stopifnot(all(d$model == "glm"))
names(d)[names(d) == "fsmethod"] <- "fs_method"
d <- d[grepl("^nonintersection", d$cohort), ]
cat("nonintersection GLM configs:", nrow(d), "\n")

thresholds <- c(Inf, seq(1, 0, by = -0.01), -Inf)
fmt_vec <- function(x) paste(sprintf("%.3f", x), collapse = "|")

label_cache <- new.env(parent = emptyenv())
get_labels <- function(set, cohort, datatype, featureset, subset) {
  key <- paste(set, cohort, datatype, featureset, subset, sep = "|")
  if (exists(key, envir = label_cache, inherits = FALSE)) return(get(key, envir = label_cache, inherits = FALSE))
  rds <- file.path("data/model_data", set, cohort, datatype, paste0(featureset, "_", subset, ".RDS"))
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
  pred_file <- file.path("res", sprintf("%s_%s_%s_%s_%s_glm_%s_v%s_pred_prob.csv",
                          r$set, r$cohort, r$datatype, r$featureset, r$subset, r$fs_method, r$varcut))
  out <- as.list(r[, c("set","cohort","datatype","featureset","subset","fs_method","varcut","mean_auc","sd_auc","n_auc")])
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
names(roc)[names(roc) == "fs_method"] <- "model"
write.csv(roc, "results_roc_curves_nonint.csv", row.names = FALSE)
cat("wrote results_roc_curves_nonint.csv rows", nrow(roc), " missing", length(missing), "\n")
cat("ROC_GLM_NONINT_DONE\n")
