#!/usr/bin/env Rscript
# Aggregate the full GLM run (10 cohorts x 3 fs-methods x ...; ~9000 configs) (mad vs cor vs mad_nosexmt) from
# rnaseq_202609/res/*_cvSummary.csv into one long comparison table.
#
# Filename pattern: <set>_<cohort>_<datatype>_<featureset>_<subset>_glm_<fsmethod>_v<varcut>_cvSummary.csv
# e.g. prev_wholeset_normalized_fpkm_clinics_12vs0_glm_mad_nosexmt_v200_cvSummary.csv
#
# Mirrors verify_auc.R's parsing approach (anchored regex extraction since cohort/
# featureset/fsmethod names themselves contain underscores) but extended for the
# extra fs_method segment introduced in this run.

suppressPackageStartupMessages(library(data.table))

resdir <- "res"
out <- "results_fsmethod_comparison.csv"

cohorts_re  <- "novel_coding|novel|nonintersection_coding|nonintersection|wholeset_coding|wholeset|stage1_coding|stage1|intersection_coding|intersection"
subsets_re  <- "123vs0|12vs0|3vs0"
fsmethod_re <- "natcommad|mad_nosexmt|mad|cor"

files <- list.files(resdir, pattern = "_cvSummary\\.csv$", full.names = TRUE)
cat("found", length(files), "cvSummary files\n")

parse_name <- function(path) {
  bn <- sub("_cvSummary\\.csv$", "", basename(path))

  varcut <- sub(".*_v", "", bn)
  rest   <- sub("_v[0-9]+$", "", bn)                                   # ..._glm_<fsmethod>

  fsmethod <- regmatches(rest, regexpr(paste0("(", fsmethod_re, ")$"), rest))
  rest     <- sub(paste0("_", fsmethod, "$"), "", rest)                # ..._glm

  model <- sub(".*_", "", rest)                                       # "glm"
  rest  <- sub("_glm$", "", rest)                                     # set_cohort_datatype_featureset_subset

  subset <- regmatches(rest, regexpr(paste0("(", subsets_re, ")$"), rest))
  rest   <- sub(paste0("_", subset, "$"), "", rest)                   # set_cohort_datatype_featureset

  set <- sub("_.*", "", rest)
  r2  <- sub("^[^_]+_", "", rest)                                     # cohort_datatype_featureset

  cohort <- regmatches(r2, regexpr(paste0("^(", cohorts_re, ")"), r2))
  r3     <- sub(paste0("^", cohort, "_"), "", r2)                     # datatype_featureset

  datatype   <- sub("_.*", "", r3)
  featureset <- sub("^[^_]+_", "", r3)

  data.frame(set, cohort, datatype, featureset, subset, model, fsmethod, varcut,
             stringsAsFactors = FALSE)
}

rows <- vector("list", length(files))
n_bad <- 0
for (i in seq_along(files)) {
  meta <- tryCatch(parse_name(files[i]), error = function(e) NULL)
  if (is.null(meta) || any(sapply(meta, function(x) length(x) == 0 || x == ""))) {
    n_bad <- n_bad + 1
    next
  }
  d <- tryCatch(fread(files[i], showProgress = FALSE), error = function(e) NULL)
  if (is.null(d) || !"AUC" %in% names(d)) { n_bad <- n_bad + 1; next }
  v <- d$AUC[is.finite(d$AUC)]
  rows[[i]] <- cbind(meta, mean_auc = mean(v), sd_auc = sd(v), n_auc = length(v))
  if (i %% 500 == 0) cat("  ...processed", i, "/", length(files), "\n")
}
res <- do.call(rbind, rows[!sapply(rows, is.null)])
cat("parsed rows:", nrow(res), " unparseable/unreadable files:", n_bad, "\n")

write.csv(res, out, row.names = FALSE)
cat("wrote", out, "\n")

# ---------------------------------------------------------------------------
# Quick summary: mean AUC by fsmethod, averaged across the whole grid, and
# separately for the "real" gene-based feature sets (fpkm, fpkm_clinics,
# module20, module20_clinics) vs the gene-free "clinics" baseline.
# ---------------------------------------------------------------------------
cat("\n=== mean AUC by fs_method, ALL configs ===\n")
print(aggregate(mean_auc ~ fsmethod, res, mean))

cat("\n=== mean AUC by fs_method, gene-based featuresets only (excl. 'clinics') ===\n")
gene_based <- res[res$featureset != "clinics", ]
print(aggregate(mean_auc ~ fsmethod, gene_based, mean))

cat("\n=== mean AUC by fs_method x set (gene-based featuresets only) ===\n")
print(aggregate(mean_auc ~ fsmethod + set, gene_based, mean))

cat("\n=== mean AUC by fs_method x subset (gene-based featuresets only) ===\n")
print(aggregate(mean_auc ~ fsmethod + subset, gene_based, mean))

cat("\n=== mean AUC by fs_method x cohort (gene-based featuresets only) ===\n")
print(aggregate(mean_auc ~ fsmethod + cohort, gene_based, mean))

cat("\n=== best (max) config per fs_method (gene-based featuresets only) ===\n")
best <- do.call(rbind, lapply(split(gene_based, gene_based$fsmethod), function(d) d[which.max(d$mean_auc), ]))
print(best)
