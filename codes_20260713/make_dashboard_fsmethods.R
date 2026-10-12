# Dashboard for the mad / cor / mad_nosexmt gene-prefilter comparison (GLM only),
# adapted from erdn_202607/codes/make_dashboard_method.R (the latest/most complete
# of that codebase's dashboard generators, by file timestamp and by using the
# richer dashboard_template_v2.html). Structural differences from the original:
#
#  - Only the "byfold" method is available here (mean of per-fold AUC over the
#    10-fold x 100-boot repeats) -- this run never built the pooled-ROC-curve or
#    per-patient-averaged aggregates the original's "pooled"/"perpatient" methods
#    need, so this script only produces that one (the pipeline's official metric
#    per the original script's own comments).
#  - 5-way "model" axis: glm_mad/glm_cor/glm_mad_nosexmt (this run's 3 gene-
#    prefilter methods, GLM only) PLUS rf_mad/xgb_mad, pulled in from the earlier
#    erdn_202607 pipeline (Section 9) which ran mad-only feature selection across
#    all three learners. glm_mad/glm_cor are "stable" (shown by default);
#    glm_mad_nosexmt/rf_mad/xgb_mad sit behind the toggle -- purely a display
#    filter (default-collapsed set of "additional variants"), not a trust flag.
#    RF/XGB results are NOT re-derived here -- they're the same fixed-direction
#    numbers already validated in Section 9, just relabeled ("rf"/"xgb" -> "rf_mad"
#    /"xgb_mad") so they don't collide with GLM's own "mad" naming once combined
#    into one long table.
#  - nature_com (the original manuscript-era model) is NOT part of this pipeline's
#    prev/current grid -- it's a separate 5-varcut series (Normal tissue only,
#    per instruction), cohort=matched_rm5years/12vs0/fpkm/normalized, computed by
#    make_nature_com_roc_series.R (fixed-direction per-fold AUC + pooled ROC
#    curve, same methodology as everything else here) from the same raw
#    predictions used throughout Section 2/3. Added as 5 extra LONG_DATA rows so
#    it's visible/filterable in the Explorer, Condition-summary, and Varcut tabs
#    (added to vcSet's dropdown in the template) and downloadable via CSV export,
#    WITHOUT forcing it into the Headline/Heatmap/Best-model tabs, which are
#    structurally a prev-vs-current comparison tool and assume exactly those two
#    `set` values throughout.
#  - ROC_DATA: pooled curves for the full prev/current GLM grid
#    (make_roc_curves_fsmethods.R, adapted from res/make_roc_curves_202607.R),
#    the RF/XGB curves from that same original res/make_roc_curves_202607.R run,
#    and nature_com's 5-varcut series -- all merged against `d$mean` so the
#    curve's displayed/ranked AUC always matches this dashboard's official
#    per-fold metric, never the pooled curve's own raw area -- exactly as the
#    ROC tab's own description (carried over unchanged from the original
#    template) explains.
#
# Run from the rnaseq_202609 root (where data/model_data/ lives) via srun/sbatch,
# never the login node.
suppressPackageStartupMessages({ library(jsonlite) })

coh_lv <- c("matched_rm5years",  # nature_com-only, not in the prev/current grid -- listed first per instruction, so its row is the top row of every Cohort x subset heatmap panel
            "wholeset","wholeset_coding","stage1","stage1_coding","intersection","intersection_coding",
            "nonintersection","nonintersection_coding")
## 5-value "model" axis: GLM's 3 fs_methods (this run) + RF/XGB from the earlier
## erdn_202607 pipeline (mad-only feature selection there -- that pipeline never
## tried cor/mad_nosexmt with RF/XGB). "glm_mad"/"glm_cor"/"glm_mad_nosexmt" so
## RF/XGB's own "mad" doesn't collide with GLM's.
mod_lv <- c("glm_mad","glm_cor","glm_mad_nosexmt","rf_mad","xgb_mad","glm_natcommad")
stable_lv <- c("glm_mad","glm_cor")
dm_lv  <- c("fpkm","normalized")
dn_lv  <- c("clinics","fpkm","fpkm_clinics","module20","module20_clinics")
cmp_lv <- c("123vs0","12vs0","3vs0")
gene_fs <- c("fpkm","fpkm_clinics")
key <- c("set","cohort","datatype","featureset","subset","model","varcut")

## ---- load prev/current GLM results, reshape to the schema the template expects ----
## NOTE: this CSV already has a "model" column (always "glm", GLM-only run) --
## drop it before promoting "fsmethod" to "model", or the rename creates a
## duplicate "model" column and d$model silently resolves to the wrong
## (always-"glm") one, filtering everything to zero rows downstream.
d <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors = FALSE)
stopifnot(all(d$model == "glm"))  # sanity: confirms which column we're safe to drop
d$model <- NULL
names(d)[names(d) == "mean_auc"] <- "mean"
names(d)[names(d) == "sd_auc"]   <- "sd"
names(d)[names(d) == "n_auc"]    <- "n"
names(d)[names(d) == "fsmethod"] <- "model"
d$model <- paste0("glm_", d$model)  # mad/cor/mad_nosexmt -> glm_mad/glm_cor/glm_mad_nosexmt
se <- d$sd / sqrt(10)
d$lo <- pmax(0, d$mean - 1.96 * se); d$hi <- pmin(1, d$mean + 1.96 * se)
d$tissue <- NA_character_

## ---- RF/XGB from the earlier erdn_202607 pipeline (mad feature selection) ----
rfxgb <- read.csv("rf_xgb_auc_202607.csv", stringsAsFactors = FALSE)
rfxgb <- rfxgb[rfxgb$model %in% c("rf","xgb"), ]
rfxgb$model <- paste0(rfxgb$model, "_mad")
se2 <- rfxgb$sd / sqrt(10)
rfxgb$lo <- pmax(0, rfxgb$mean - 1.96 * se2); rfxgb$hi <- pmin(1, rfxgb$mean + 1.96 * se2)
rfxgb$tissue <- NA_character_

## ---- nature_com reference rows: full 5-varcut series, Normal tissue only ----
## (built by make_nature_com_roc_series.R from the same raw res_withmodel
## predictions used in Section 2; v200 here reproduces 0.8283/CI 0.740-0.917
## exactly, confirming consistency with the earlier hand-computed number)
nature_com <- read.csv("nature_com_normal_varcut_series_auc.csv", stringsAsFactors = FALSE)
nature_com$model <- paste0("glm_", nature_com$model)  # "mad" -> "glm_mad", consistent with the rest
nature_com$tissue <- "Normal"

d <- rbind(d[, c(key,"mean","sd","n","lo","hi","tissue")],
           rfxgb[, c(key,"mean","sd","n","lo","hi","tissue")],
           nature_com[, c(key,"mean","sd","n","lo","hi","tissue")])

d <- d[!is.na(d$mean), ]
for (col in intersect(c("mean","sd","lo","hi"), names(d))) d[[col]] <- round(d[[col]], 4)

## ---- LONG_DATA ----
long_cols <- c("set","cohort","datatype","featureset","subset","model","varcut","mean","sd","n","lo","hi")
long_df <- d[, long_cols]
long_json <- toJSON(long_df, dataframe = "values", auto_unbox = TRUE, na = "null")

## ---- headline (best/median AUC per model x set), stable (mad+cor) vs all ----
## generic over set -- nature_com (model="mad" only) shows up as its own bar
## group automatically, alongside prev/current, wherever it has data
mk_headline <- function(models) {
  dd <- d[d$model %in% models, ]
  agg_best <- aggregate(mean ~ model + set, dd, max)
  agg_med  <- aggregate(mean ~ model + set, dd, median)
  merge(agg_best, agg_med, by = c("model","set"), suffixes = c("_best","_median"))
}
headline_stable_json <- toJSON(mk_headline(stable_lv), dataframe = "rows", auto_unbox = TRUE, na = "null")
headline_all_json    <- toJSON(mk_headline(mod_lv),    dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- best model per cohort x subset x featureset, stable vs all (prev/current only) ----
md_base <- "data/model_data/current"
rec_n <- function(x) { r <- as.numeric(as.character(x$recurrence)); c(tot=length(r), pos=sum(r==1), neg=sum(r==0)) }
noclin <- list(); clin <- list(); nf_clinics <- list()
for (coh in setdiff(coh_lv, "matched_rm5years")) for (sub in cmp_lv) {
  k <- paste(coh, sub, sep = "|")
  m20  <- readRDS(file.path(md_base, coh, "fpkm", paste0("module20_", sub, ".RDS")))
  m20c <- readRDS(file.path(md_base, coh, "fpkm", paste0("module20_clinics_", sub, ".RDS")))
  cl   <- readRDS(file.path(md_base, coh, "fpkm", paste0("clinics_", sub, ".RDS")))
  noclin[[k]] <- rec_n(m20); clin[[k]] <- rec_n(m20c); nf_clinics[[k]] <- ncol(cl) - 1
}
with_clin <- c("clinics","fpkm_clinics","module20_clinics")

best_one <- function(dd, set, coh, sub, fs) {
  x <- dd[dd$set==set & dd$cohort==coh & dd$subset==sub & dd$featureset==fs, ]
  if (!nrow(x)) return(data.frame(auc=NA, model=NA, datatype=NA, varcut=NA, ci=NA, stringsAsFactors=FALSE))
  i <- which.max(x$mean)
  ci_txt <- sprintf("%.2f (%.2f-%.2f)", x$mean[i], x$lo[i], x$hi[i])
  data.frame(auc = round(x$mean[i],3), model = x$model[i], datatype = x$datatype[i],
             varcut = x$varcut[i], ci = ci_txt, stringsAsFactors = FALSE)
}
## nature_com's known sample counts, per subset (matched_rm5years, Normal tissue)
## -- from make_nature_com_roc_series.R's own label-file counts, confirmed
## against Section 2/3 for 12vs0. Not readable from data/model_data (that's
## prev/current only), so hardcoded here rather than left blank.
nc_cnt_by_sub <- list(
  "12vs0"  = c(tot = 91,  pos = 23, neg = 68),
  "123vs0" = c(tot = 113, pos = 45, neg = 68),
  "3vs0"   = c(tot = 90,  pos = 22, neg = 68)
)

mk_best <- function(models) {
  dd <- d[d$model %in% models, ]
  rows <- list()
  for (coh in coh_lv) for (sub in cmp_lv) for (fs in dn_lv) {
    k <- paste(coh, sub, sep = "|")
    is_nc <- coh == "matched_rm5years"
    if (is_nc) {
      # nature_com only ever has fpkm-featureset data (all 3 subsets, per
      # make_nature_com_roc_series.R) -- skip the other 4 featuresets for this
      # cohort entirely rather than emitting empty rows.
      if (fs != "fpkm") next
      cnt <- nc_cnt_by_sub[[sub]]; nfeat_fixed <- NA
    } else {
      cnt <- if (fs %in% with_clin) clin[[k]] else noclin[[k]]
      nfeat_fixed <- switch(fs, clinics = nf_clinics[[k]], module20 = 170, module20_clinics = 175, NA)
    }
    disp <- function(vc) if (fs %in% gene_fs) as.character(vc) else as.character(nfeat_fixed)
    p  <- best_one(dd, "prev", coh, sub, fs)
    c  <- best_one(dd, "current", coh, sub, fs)
    nc <- best_one(dd, "nature_com", coh, sub, fs)
    rows[[length(rows)+1]] <- data.frame(
      cohort = coh, subset = sub, featureset = fs,
      n_total = unname(cnt["tot"]), n_recur = unname(cnt["pos"]), n_nonrecur = unname(cnt["neg"]),
      prev_auc = p$auc, prev_auc_ci = p$ci, prev_model = p$model, prev_datatype = p$datatype,
      prev_nfeat = disp(p$varcut), prev_varcut = p$varcut,
      current_auc = c$auc, current_auc_ci = c$ci, current_model = c$model, current_datatype = c$datatype,
      current_nfeat = disp(c$varcut), current_varcut = c$varcut,
      nature_com_auc = nc$auc, nature_com_auc_ci = nc$ci, nature_com_varcut = nc$varcut,
      delta_current_minus_prev = round(c$auc - p$auc, 3),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
best_stable <- mk_best(stable_lv); best_all <- mk_best(mod_lv)
best_stable_json <- toJSON(best_stable, dataframe = "rows", auto_unbox = TRUE, na = "null")
best_all_json    <- toJSON(best_all,    dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- condition summary (family table): set x cohort x featureset x datatype x model ----
## generic over `set`, so nature_com's 2 rows DO show up here automatically
mk_family <- function(models) {
  dd <- d[d$model %in% models, ]
  mx  <- aggregate(mean ~ set + cohort + featureset + datatype + model, dd, max)
  med <- aggregate(mean ~ set + cohort + featureset + datatype + model, dd, median)
  cnt <- aggregate(mean ~ set + cohort + featureset + datatype + model, dd, length)
  fam <- merge(merge(mx, med, by = c("set","cohort","featureset","datatype","model"), suffixes = c("_max","_median")),
               cnt, by = c("set","cohort","featureset","datatype","model"))
  names(fam)[names(fam)=="mean"] <- "n_configs"
  names(fam)[names(fam)=="mean_max"] <- "max_auc"
  names(fam)[names(fam)=="mean_median"] <- "median_auc"
  fam$max_auc <- round(fam$max_auc, 3); fam$median_auc <- round(fam$median_auc, 3)
  fam
}
family_stable <- mk_family(stable_lv); family_all <- mk_family(mod_lv)
family_stable_json <- toJSON(family_stable, dataframe = "rows", auto_unbox = TRUE, na = "null")
family_all_json    <- toJSON(family_all,    dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- cohorts tab (prev/current only; nature_com's single config described in the banner instead) ----
gc1 <- function(set, coh, dt) ncol(readRDS(file.path("data/model_data", set, coh, dt, "fpkm_123vs0.RDS"))) - 1
gAf_c <- gc1("current","wholeset","fpkm");        gAn_c <- gc1("current","wholeset","normalized")
gAf_p <- gc1("prev","wholeset","fpkm");           gAn_p <- gc1("prev","wholeset","normalized")
gC_c  <- gc1("current","wholeset_coding","fpkm"); gC_p  <- gc1("prev","wholeset_coding","fpkm")
sW <- noclin[["wholeset|123vs0"]]; sS <- noclin[["stage1|123vs0"]]; sI <- noclin[["intersection|123vs0"]]
sN <- noclin[["nonintersection|123vs0"]]
cohorts_tab <- data.frame(
  cohort = setdiff(coh_lv, "matched_rm5years"),
  n_total = c(sW["tot"],sW["tot"], sS["tot"],sS["tot"], sI["tot"],sI["tot"], sN["tot"],sN["tot"]),
  n_pos   = c(sW["pos"],sW["pos"], sS["pos"],sS["pos"], sI["pos"],sI["pos"], sN["pos"],sN["pos"]),
  n_neg   = c(sW["neg"],sW["neg"], sS["neg"],sS["neg"], sI["neg"],sI["neg"], sN["neg"],sN["neg"]),
  description = c(
    sprintf("All recurrence-labelled samples (no stage filter). Full gene set: current %d (fpkm) / %d (normalized) genes; prev %d / %d.", gAf_c, gAn_c, gAf_p, gAn_p),
    sprintf("Same samples as wholeset, restricted to protein-coding genes: current %d -> %d, prev %d -> %d genes.", gAf_c, gC_c, gAf_p, gC_p),
    "Only stage IA1/IA2/IB samples; full gene set (same gene counts as wholeset).",
    "Stage IA1/IA2/IB samples, protein-coding genes only (same gene counts as wholeset_coding).",
    "Stage-I samples that also overlap the prior 246-sample paper cohort; full gene set (same gene counts as wholeset).",
    "Intersection samples, protein-coding genes only (same gene counts as wholeset_coding).",
    "Wholeset samples that do NOT overlap the prior 246-sample paper cohort (wholeset minus intersection). Note this is not fully nature_com-independent: 18 nature_com subjects still remain in it. Full gene set (same gene counts as wholeset).",
    "Nonintersection samples, protein-coding genes only (same gene counts as wholeset_coding)."),
  stringsAsFactors = FALSE)
cohorts_json <- toJSON(cohorts_tab, dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- ROC_DATA: pooled curves for prev/current GLM grid + RF/XGB + nature_com series ----
## Same pattern as make_dashboard_method.R: drop the curve file's own auc/n and
## re-merge against d$mean, so curve labeling/ranking always matches this
## dashboard's official per-fold AUC, never the pooled curve's own raw area.
roc_pc <- read.csv("results_roc_curves_fsmethods.csv", stringsAsFactors = FALSE)
roc_pc$model <- paste0("glm_", roc_pc$model)  # mad/cor/mad_nosexmt -> glm_mad/glm_cor/glm_mad_nosexmt, matching d
roc_nc <- read.csv("nature_com_normal_varcut_series_roc.csv", stringsAsFactors = FALSE)
roc_nc$model <- paste0("glm_", roc_nc$model)  # "mad" -> "glm_mad", must match d's relabeling above
roc_rx <- read.csv("rf_xgb_roc_202607.csv", stringsAsFactors = FALSE)
roc_rx <- roc_rx[roc_rx$model %in% c("rf","xgb"), ]
roc_rx$model <- paste0(roc_rx$model, "_mad")
roc <- rbind(roc_pc[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")],
             roc_rx[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")],
             roc_nc[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")])
roc <- merge(roc, d[, c(key,"mean")], by = key, all.x = TRUE)
names(roc)[names(roc) == "mean"] <- "auc"
roc <- roc[!is.na(roc$auc) & !is.na(roc$fpr), ]
roc$auc <- round(roc$auc, 4)
roc <- roc[, c(key, "auc","n_pred","n_pos_pred","n_neg_pred","fpr","tpr")]
roc_json <- toJSON(roc, dataframe = "rows", auto_unbox = TRUE, na = "null", digits = NA)

## ---- gene correlation explorer data (current RNA-seq only, see make_gene_correlation_data.R) ----
gene_master_json     <- toJSON(read.csv("gene_corr_master.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
gene_membership_json <- toJSON(read.csv("gene_corr_membership.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
gene_scatter_json    <- toJSON(read.csv("gene_corr_scatter.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- sample correlation heatmap images (prev-highreads QC + nature_com vs prev/current) ----
## Static PNGs (converted from existing PDFs via ghostscript, base64-embedded)
## -- not per-config generated data, just a pre-built JSON manifest.
## See erdn_20260709/correlation/ (Section 7) and erdn_202502/result/correlation_highreads/.
corr_plots_json <- paste(readLines("corr_plots.json", warn = FALSE), collapse = "")

## ---- top-MAD gene overlap (nature_com vs prev vs current), 3-way Venn by varcut ----
## See make_venn_data.R / erdn_20260709/correlation/venn_madgenes.csv
venn_json <- toJSON(read.csv("venn_madgenes.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")

## ---- meta ----
meta <- list(
  generated = as.character(Sys.Date()),
  method = "byfold",
  method_label = "Per-fold AUC averaged over 10-fold x 100-boot repeats (same OFFICIAL metric as the erdn_202607 dashboards). 5-way model axis: glm_mad/glm_cor/glm_mad_nosexmt (this run, GLM only, 3 gene-prefilter methods) + rf_mad/xgb_mad (earlier erdn_202607 pipeline, mad-only feature selection).",
  source = "results_fsmethod_comparison.csv + results_roc_curves_fsmethods.csv (rnaseq_202609, GLM) + rf_xgb_auc_202607.csv + rf_xgb_roc_202607.csv (erdn_202607, RF/XGB) + nature_com_normal_varcut_series_{auc,roc}.csv (5-varcut Normal-tissue series)",
  n_configs = nrow(d),
  stable_configs = sum(d$model %in% stable_lv),
  provisional_xgb_configs = sum(!d$model %in% stable_lv),  # reusing the template's toggle-count slot for "everything behind the toggle" (glm_mad_nosexmt + rf_mad + xgb_mad)
  models = mod_lv, stable_models = stable_lv,
  cohorts = coh_lv, datatypes = dm_lv, featuresets = dn_lv, subsets = cmp_lv,
  varcuts = sort(unique(d$varcut))
)
meta_json <- toJSON(meta, auto_unbox = TRUE)

cat("long:", nrow(long_df), " best_stable:", nrow(best_stable), " best_all:", nrow(best_all),
    " family_stable:", nrow(family_stable), " family_all:", nrow(family_all),
    " cohorts:", nrow(cohorts_tab), " nature_com rows:", sum(d$set=="nature_com"),
    " roc curves:", nrow(roc), " roc nature_com:", sum(roc$set=="nature_com"),
    " gene_master:", nrow(read.csv("gene_corr_master.csv")),
    " gene_scatter combos:", nrow(read.csv("gene_corr_scatter.csv")), "\n")

## ---- Plotly.js, inlined for full offline use (no CDN dependency) ----
plotly_js <- paste(readLines("codes_20260713/plotly-2.32.0.min.js", warn = FALSE), collapse = "\n")

## ---- assemble ----
tpl <- paste(readLines("codes_20260713/dashboard_template_fsmethods.html", warn = FALSE), collapse = "\n")
tpl <- sub("__PLOTLY_JS__", plotly_js, tpl, fixed = TRUE)
tpl <- sub("__LONG_DATA__", long_json, tpl, fixed = TRUE)
tpl <- sub("__COHORTS_TAB__", cohorts_json, tpl, fixed = TRUE)
tpl <- sub("__HEADLINE_STABLE__", headline_stable_json, tpl, fixed = TRUE)
tpl <- sub("__HEADLINE_ALL__", headline_all_json, tpl, fixed = TRUE)
tpl <- sub("__BEST_STABLE__", best_stable_json, tpl, fixed = TRUE)
tpl <- sub("__BEST_ALL__", best_all_json, tpl, fixed = TRUE)
tpl <- sub("__FAMILY_STABLE__", family_stable_json, tpl, fixed = TRUE)
tpl <- sub("__FAMILY_ALL__", family_all_json, tpl, fixed = TRUE)
tpl <- sub("__ROC_DATA__", roc_json, tpl, fixed = TRUE)
tpl <- sub("__META__", meta_json, tpl, fixed = TRUE)
tpl <- sub("__GENE_MASTER__", gene_master_json, tpl, fixed = TRUE)
tpl <- sub("__GENE_MEMBERSHIP__", gene_membership_json, tpl, fixed = TRUE)
tpl <- sub("__GENE_SCATTER__", gene_scatter_json, tpl, fixed = TRUE)
tpl <- sub("__CORR_PLOTS__", corr_plots_json, tpl, fixed = TRUE)
tpl <- sub("__VENN_DATA__", venn_json, tpl, fixed = TRUE)

out_file <- "htmls_20260713/auc_dashboard_20260713_fsmethods.html"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(tpl, out_file)
cat("Wrote", out_file, "\n")

write.csv(d, "htmls_20260713/dashboard_long_data_20260713.csv", row.names = FALSE)
cat("Wrote htmls_20260713/dashboard_long_data_20260713.csv\n")
