# Gene correlation explorer (current RNA-seq reprocessing only -- scope matches
# the request, which named conditions/cohorts, not the prev/current axis).
# For each (cohort, subset), rank genes by |Pearson r| with the binary
# recurrence outcome (same statistic as the existing "cor" fs_method elsewhere
# in this pipeline -- cor(gene, as.numeric(recurrence)), pairwise-complete),
# take the top 100, and union those 18 lists into one master gene list. Then
# build a compact scatter dataset (expr/recurrence/subj, pipe-delimited --
# same convention as make_roc_curves_fsmethods.R's fpr/tpr strings) for that
# master gene list across every (cohort, subset), so switching gene/cohort/
# subset in the dashboard is instant client-side lookup, no reload.
rm(list = ls())

# --- CONFIG (override via env vars) -------------------------------------------
# GENECORR_RAW must contain the per-(cohort,subset) current-RNA-seq expression
# subsets this script reads. The original run used an ephemeral session scratchpad
# that no longer exists; set this to a persistent location holding regenerated
# inputs before re-running. Default: a provenance_inputs/ folder beside res_20260713.
proj <- Sys.getenv("ERDN_PROJ", ".")
base <- Sys.getenv("GENECORR_RAW", file.path(proj, "provenance_inputs", "genecorr_raw", "current"))
out_dir <- Sys.getenv("GENECORR_OUT", file.path(proj, "htmls"))
stopifnot(dir.exists(base))  # fail loudly instead of silently reading nothing
# ------------------------------------------------------------------------------

cohorts <- c("wholeset","wholeset_coding","stage1","stage1_coding","intersection","intersection_coding")
subsets <- c("12vs0","123vs0","3vs0")
topk <- 100

membership <- list()
loaded <- list()

for (coh in cohorts) for (sub in subsets) {
  f <- file.path(base, coh, "fpkm", paste0("fpkm_", sub, ".RDS"))
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); next }
  x <- readRDS(f)
  y <- as.numeric(x$recurrence)
  genes <- setdiff(colnames(x), "recurrence")
  mat <- x[, genes, drop = FALSE]
  r <- apply(mat, 2, function(g) suppressWarnings(cor(as.numeric(g), y, use = "pairwise.complete.obs")))
  r[is.na(r)] <- 0
  ord <- order(-abs(r))[seq_len(min(topk, length(r)))]
  key <- paste(coh, sub)
  membership[[key]] <- data.frame(cohort = coh, subset = sub, gene = genes[ord],
                                   r = round(r[ord], 4), rank = seq_along(ord), stringsAsFactors = FALSE)
  loaded[[key]] <- x
  cat(sprintf("%s: n=%d, top |r| range %.3f-%.3f\n", key, nrow(x), abs(r[ord[length(ord)]]), abs(r[ord[1]])))
}

membership_tab <- do.call(rbind, membership)
master_genes <- sort(unique(membership_tab$gene))
cat("\n=== master gene union:", length(master_genes), "genes, from", nrow(membership_tab), "membership rows ===\n")
write.csv(membership_tab, file.path(out_dir, "gene_corr_membership.csv"), row.names = FALSE)

## per-gene summary for the picker: how many combos, best |r|, mean |r|
agg_n   <- aggregate(r ~ gene, membership_tab, length)
agg_max <- aggregate(r ~ gene, membership_tab, function(v) max(abs(v)))
master_tab <- merge(agg_n, agg_max, by = "gene", suffixes = c("_n","_maxr"))
names(master_tab) <- c("gene","n_combos","max_abs_r")
master_tab <- master_tab[order(-master_tab$n_combos, -master_tab$max_abs_r), ]
master_tab$max_abs_r <- round(master_tab$max_abs_r, 3)
write.csv(master_tab, file.path(out_dir, "gene_corr_master.csv"), row.names = FALSE)

## ---- compact scatter data: master gene list x every (cohort, subset) ----
fmt_num <- function(x) paste(sprintf("%.2f", x), collapse = "|")
fmt_int <- function(x) paste(x, collapse = "|")
fmt_chr <- function(x) paste(x, collapse = "|")

scatter_rows <- list()
for (key in names(loaded)) {
  x <- loaded[[key]]
  parts <- strsplit(key, " ")[[1]]
  coh <- parts[1]; sub <- parts[2]
  present <- intersect(master_genes, colnames(x))
  y <- as.integer(x$recurrence)
  subj <- rownames(x)
  for (g in present) {
    e <- as.numeric(x[[g]])
    e0 <- e[y == 0]; e1 <- e[y == 1]
    ## Wilcoxon rank-sum (Mann-Whitney) p for the group difference -- pairs
    ## naturally with a boxplot (nonparametric, robust to FPKM's typical
    ## right-skew), unlike the Pearson r above which is a linear-association
    ## statistic, not a group-difference test.
    ## suppressWarnings: wilcox.test warns "cannot compute exact p-value with
    ## ties" for tied FPKM values (common, e.g. zeros) and silently falls back
    ## to the normal approximation -- expected and harmless here, not a bug.
    pval <- tryCatch(suppressWarnings(wilcox.test(e0, e1)$p.value), error = function(e) NA_real_)
    scatter_rows[[length(scatter_rows) + 1]] <- data.frame(
      cohort = coh, subset = sub, gene = g,
      expr = fmt_num(e), recur = fmt_int(y), subj = fmt_chr(subj),
      n = length(y), median0 = round(median(e0), 2), median1 = round(median(e1), 2),
      pval = signif(pval, 3), stringsAsFactors = FALSE)
  }
}
scatter_tab <- do.call(rbind, scatter_rows)
cat("scatter_tab rows (gene x cohort x subset combos):", nrow(scatter_tab), "\n")
write.csv(scatter_tab, file.path(out_dir, "gene_corr_scatter.csv"), row.names = FALSE)
cat("wrote gene_corr_membership.csv, gene_corr_master.csv, gene_corr_scatter.csv\n")

file_sizes <- file.info(file.path(out_dir, c("gene_corr_membership.csv","gene_corr_master.csv","gene_corr_scatter.csv")))$size
cat("file sizes (MB):", round(file_sizes / 1e6, 2), "\n")
