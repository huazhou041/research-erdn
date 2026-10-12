#!/usr/bin/env Rscript
# ============================================================================
# Score a NEW cohort with a chosen full-cohort ("final") model from this project.
#
# Each config's saved .rds IS the model fit on ALL available data (the CV was only
# used to estimate performance). This script loads one config's final model, aligns
# a new expression matrix to that model's genes, and outputs a recurrence probability
# per sample, plus a class call at a chosen threshold.
#
# USAGE:
#   Rscript predict_new_cohort.R <config_key> <new_expr.csv> <threshold> <out.csv>
#
#   <config_key>  e.g. current_wholeset_normalized_fpkm_123vs0_glm_mad_v200
#                 (the key column in model_operating_points_202609.csv; pick a
#                  WELL-POWERED, non-low-event config -- see that table.)
#   <new_expr.csv> genes x samples matrix: first column = gene symbol, remaining
#                 columns = samples. MUST be normalized the SAME way as training
#                 (the training matrices are data/model_data[_smoking]/<set>/<cohort>/
#                 <datatype>/, i.e. DESeq2 'normalized' or 'fpkm' as in the key).
#   <threshold>   probability cutoff for the positive (recurrence) class. Get a
#                 starting value from model_operating_points_202609.csv
#                 (pooled_youden_cutoff) -- but RE-CALIBRATE + EXTERNALLY VALIDATE
#                 before any real use (see README caveats).
#
# IMPORTANT CAVEATS (see README):
#   * Probabilities are NOT guaranteed calibrated -> a fixed probability threshold
#     transfers only if you calibrate (Platt/isotonic) on a reference set.
#   * Thresholds chosen on the training/CV data are optimistic; validate externally.
#   * Do NOT deploy low-event configs (3vs0, tiny cohorts) -- discrimination itself
#     is not established there.
# ============================================================================
suppressMessages({ library(glmnet) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: predict_new_cohort.R <config_key> <new_expr.csv> <threshold> <out.csv>")
key <- args[1]; exprf <- args[2]; thr <- as.numeric(args[3]); outf <- args[4]
base <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609"

# locate the config's saved final-model .rds across the result dirs
cands <- file.path(base, c("res","res_rfxgb","res_smoking","res_rfxgb_smoking"), paste0(key, ".rds"))
rds <- cands[file.exists(cands)][1]
if (is.na(rds)) stop("model .rds not found for key: ", key, " (looked in res*/, res_rfxgb*/)")
is_tree <- grepl("_rf_v|_xgb_v", key)
mdl <- readRDS(rds)

# load + align new expression (genes x samples) to the model's feature genes
e <- read.csv(exprf, check.names = FALSE, stringsAsFactors = FALSE)
names(e)[1] <- "gene"
emat <- as.matrix(e[, -1, drop = FALSE]); rownames(emat) <- e$gene        # genes x samples
X <- t(emat)                                                              # samples x genes

CLIN <- c("Tissue_Sex","Tissue_Age","Tissue_race","Tissue_smoke_status","Tissue_Tumor_size","Tissue_Module20")
if (is_tree) {
  feats <- as.character(mdl$features)
} else {
  b <- as.matrix(coef(mdl, s = "lambda.min")); feats <- setdiff(rownames(b), "(Intercept)")
}
gene_feats <- setdiff(feats, CLIN)
missing <- setdiff(gene_feats, colnames(X))
if (length(missing))
  message(sprintf("WARNING: %d of %d model genes absent in new data; filling 0 (align normalization!).",
                  length(missing), length(gene_feats)))
Xf <- matrix(0, nrow = nrow(X), ncol = length(feats), dimnames = list(rownames(X), feats))
common <- intersect(gene_feats, colnames(X)); Xf[, common] <- X[, common]
if (any(CLIN %in% feats))
  message("NOTE: this config also used clinical covariates (", paste(intersect(CLIN,feats),collapse=","),
          "); supply them as rows in new_expr.csv or they are treated as 0.")

if (is_tree) {
  m <- mdl$model
  if (inherits(m, "xgb.Booster")) { prob <- predict(m, Xf) }
  else { prob <- predict(m, Xf, type = "prob")[, "1"] }
} else {
  prob <- as.numeric(predict(mdl, newx = Xf, s = "lambda.min", type = "response"))
}
pred <- ifelse(prob >= thr, 1L, 0L)
res <- data.frame(sample = rownames(X), prob_recurrence = round(prob, 4), predicted_class = pred)
write.csv(res, outf, row.names = FALSE)
cat(sprintf("scored %d samples with %s at threshold %.2f -> %s\n", nrow(res), key, thr, outf))
cat(sprintf("  positives called: %d / %d\n", sum(pred == 1), length(pred)))
