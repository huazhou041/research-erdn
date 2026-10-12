#!/usr/bin/env Rscript
# ============================================================================
# Score a NEW cohort with a chosen full-cohort ("final") model, WITH calibration.
#
# Pipeline (matches the dashboard "Calibrated CV (deployable)" basis):
#   1. load the config's full-data model (its saved .rds IS the fit on ALL data)
#   2. score the new cohort -> raw probability
#   3. CALIBRATE: fit Platt scaling (logistic regression) on the config's OUT-OF-FOLD
#      CV predictions (honest: each patient scored by a model that didn't see it),
#      apply that smooth monotone mapping to the new raw scores -> calibrated probability
#   4. classify at <threshold> (a CALIBRATED-probability cutoff, e.g. the Youden
#      value from the "Calibrated CV" slider in the dashboard)
#
# USAGE:
#   Rscript predict_new_cohort.R <config_key> <new_expr.csv> <calibrated_threshold> <out.csv>
#     <config_key>  key from model_operating_points_202609.csv (use a WELL-POWERED,
#                   non-low-event, gene-based config: fpkm / fpkm_clinics)
#     <new_expr.csv> genes x samples, first col = gene symbol; normalized the SAME
#                   way as the training matrix named in the key (normalized / fpkm)
#     <calibrated_threshold> calibrated-probability cutoff from the dashboard slider
#
# CAVEATS: calibration (Platt) is fit on this cohort's OOF predictions (approximates the
# deployed model's scale) -> still VALIDATE on an independent cohort. Do NOT use
# low-event configs. Probabilities out are calibrated but not guaranteed transportable.
# ============================================================================
suppressMessages({ library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: predict_new_cohort.R <config_key> <new_expr.csv> <calibrated_threshold> <out.csv>")
key <- args[1]; exprf <- args[2]; thr <- as.numeric(args[3]); outf <- args[4]
base <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"
is_tree <- grepl("_rf_v|_xgb_v", key)
resdirs <- c("res","res_rfxgb","res_smoking","res_rfxgb_smoking")

# --- locate model + its training-label RDS + its OOF predictions ---
rds <- file.path(base, resdirs, paste0(key, ".rds")); rds <- rds[file.exists(rds)][1]
if (is.na(rds)) stop("model .rds not found for key: ", key)
resdir <- basename(dirname(rds))
pred_oof <- file.path(base, resdir, paste0(key, "_pred_prob.csv"))
# training RDS: data[_smoking]/<set>/<cohort>/<datatype>/<featureset>_<subset>.RDS
toks <- strsplit(key, "_")[[1]]
mdl <- readRDS(rds)
CLIN <- c("Tissue_Sex","Tissue_Age","Tissue_race","Tissue_smoke_status","Tissue_Tumor_size","Tissue_Module20")

score <- function(model, X) {
  if (is_tree) {
    feats <- as.character(model$features); feats <- feats[feats %in% colnames(X)]
    m <- model$model
    if (inherits(m, "xgb.Booster")) predict(m, X[, feats, drop=FALSE])
    else predict(m, as.data.frame(X[, feats, drop=FALSE]), type = "prob")[, "1"]
  } else {
    b <- as.matrix(coef(model, s = "lambda.min")); fn <- setdiff(rownames(b), "(Intercept)")
    fn <- fn[fn %in% colnames(X)]
    as.numeric(predict(model, newx = X[, fn, drop=FALSE], s = "lambda.min", type = "response"))
  }
}

# --- new cohort -> raw prob ---
e <- read.csv(exprf, check.names = FALSE, stringsAsFactors = FALSE); names(e)[1] <- "gene"
emat <- as.matrix(e[, -1, drop = FALSE]); rownames(emat) <- e$gene
Xnew <- t(emat)
raw_new <- score(mdl, Xnew)

# --- fit Platt (logistic) calibration on OOF predictions, apply to new raw scores ---
# Platt = glm(label ~ raw_oof_prob, binomial): smooth, strictly monotone (preserves
# ranking/AUC), standard for small n. Calibrated prob = plogis(a + b*raw).
platt <- NULL
if (file.exists(pred_oof)) {
  oof <- read.csv(pred_oof, stringsAsFactors = FALSE)
  ag <- aggregate(prob ~ subj, oof, mean)                 # per-subject OOF mean
  # labels: match any model_data RDS whose rownames cover these subjects and has 'recurrence'
  lab <- NULL
  cand <- Sys.glob(file.path(base, "data", "model_data*", "*", "*", "*", "*.RDS"))
  for (f in cand) {
    d <- tryCatch(readRDS(f), error = function(x) NULL); if (is.null(d) || !"recurrence" %in% names(d)) next
    if (all(ag$subj %in% rownames(d))) { lab <- setNames(as.numeric(as.character(d$recurrence)), rownames(d)); break }
  }
  if (!is.null(lab)) {
    y <- lab[as.character(ag$subj)]; ok <- is.finite(y) & is.finite(ag$prob)
    fit <- tryCatch(suppressWarnings(glm(y[ok] ~ ag$prob[ok], family = binomial())), error = function(x) NULL)
    if (!is.null(fit) && isTRUE(fit$converged)) platt <- as.numeric(coef(fit))  # c(intercept, slope)
  }
}
if (is.null(platt)) {
  message("WARNING: could not fit Platt calibration (OOF preds or labels missing / non-convergent); outputting RAW probabilities. Treat <threshold> as a RAW cutoff.")
  cal_new <- raw_new
} else {
  cal_new <- plogis(platt[1] + platt[2] * raw_new)
}

pred <- ifelse(cal_new >= thr, 1L, 0L)
res <- data.frame(sample = rownames(Xnew), prob_raw = round(raw_new,4),
                  prob_calibrated = round(cal_new,4), predicted_class = pred)
write.csv(res, outf, row.names = FALSE)
cat(sprintf("scored %d samples with %s; platt_calibrated=%s; cutoff %.2f -> %d positive\n",
            nrow(res), key, !is.null(platt), thr, sum(pred==1)))
