#!/usr/bin/env Rscript
# Aggregate res_rfxgb/*_cvSummary.csv (RF/XGB pipeline) into the same long schema
# used by rf_xgb_auc_202607.csv, WITH the nonintersection cohort included.
#
# Filename: <set>_<cohort>_<datatype>_<featureset>_<subset>_<model>_v<varcut>_cvSummary.csv
# e.g. current_nonintersection_fpkm_fpkm_clinics_12vs0_rf_v1000_cvSummary.csv
#
# `mean` is the average of the per-fold AUC column across all 10 folds x 100 boot
# reps (the 100*10 estimator used everywhere else in the dashboard). The upstream
# pipeline already computes each fold's AUC with fixed direction="<"
# (tune_drug_rf.R), so nothing is re-derived here.
suppressPackageStartupMessages(library(data.table))

resdir <- "res_rfxgb_smoking"
out    <- "rf_xgb_auc_smoking.csv"

# longest-first: nonintersection* must be tried before intersection*
cohorts_re <- "novel_coding|novel|nonintersection_coding|nonintersection|wholeset_coding|wholeset|stage1_coding|stage1|intersection_coding|intersection"
subsets_re <- "123vs0|12vs0|3vs0"
models_re  <- "glm|rf|xgb"

parse_name <- function(path) {
  bn <- sub("_cvSummary\\.csv$", "", basename(path))
  varcut <- sub(".*_v", "", bn)
  rest   <- sub("_v[0-9]+$", "", bn)                       # ..._<model>
  model  <- regmatches(rest, regexpr(paste0("(", models_re, ")$"), rest))
  if (!length(model)) return(NULL)
  rest   <- sub(paste0("_", model, "$"), "", rest)         # set_cohort_datatype_featureset_subset
  subset <- regmatches(rest, regexpr(paste0("(", subsets_re, ")$"), rest))
  if (!length(subset)) return(NULL)
  rest   <- sub(paste0("_", subset, "$"), "", rest)
  set    <- sub("_.*", "", rest)
  r2     <- sub("^[^_]+_", "", rest)                       # cohort_datatype_featureset
  cohort <- regmatches(r2, regexpr(paste0("^(", cohorts_re, ")"), r2))
  if (!length(cohort)) return(NULL)
  r3     <- sub(paste0("^", cohort, "_"), "", r2)          # datatype_featureset
  datatype   <- sub("_.*", "", r3)
  featureset <- sub("^[^_]+_", "", r3)
  data.frame(set, cohort, datatype, featureset, subset, model, varcut,
             stringsAsFactors = FALSE)
}

files <- list.files(resdir, pattern = "_cvSummary\\.csv$", full.names = TRUE)
cat("found", length(files), "cvSummary files in", resdir, "\n")

rows <- vector("list", length(files)); n_bad <- 0
for (i in seq_along(files)) {
  meta <- tryCatch(parse_name(files[i]), error = function(e) NULL)
  if (is.null(meta)) { n_bad <- n_bad + 1; next }
  d <- tryCatch(fread(files[i], showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || !"AUC" %in% names(d)) { n_bad <- n_bad + 1; next }
  v <- d$AUC[is.finite(d$AUC)]
  if (!length(v)) { n_bad <- n_bad + 1; next }
  rows[[i]] <- cbind(meta, mean = mean(v), sd = sd(v), n = length(v))
  if (i %% 500 == 0) cat("  ...", i, "/", length(files), "\n")
}
res <- do.call(rbind, rows[!sapply(rows, is.null)])
cat("parsed:", nrow(res), " unparseable:", n_bad, "\n")

write.csv(res, out, row.names = FALSE)
cat("wrote", out, "\n\n")
cat("=== rows by cohort x model ===\n")
print(table(res$cohort, res$model))
cat("\n=== nonintersection coverage by set x datatype x subset (rf+xgb only) ===\n")
ni <- res[grepl("^nonintersection", res$cohort) & res$model %in% c("rf","xgb"), ]
print(table(ni$set, ni$subset, ni$model))
