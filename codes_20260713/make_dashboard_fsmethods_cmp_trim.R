# Combined FULL-vs-EVER-SMOKER dashboard: one HTML with a top-bar checkbox that
# swaps ALL data in every tab between the full cohort and the ever-smoker subset.
# Builds every population-specific blob twice (full + smoke) and injects both into
# dashboard_template_fsmethods_cmp.html (which carries the swap JS). Gene/corr/venn
# panels are population-independent provenance and loaded once.
#
# Run from the rnaseq_202609 root via srun/sbatch (never the login node).
suppressPackageStartupMessages({ library(jsonlite) })

## ---- TRIM VARIANT: normalized-only, 3vs0 removed -------------------------
## Shadow read.csv so EVERY data frame with a datatype/subset column is filtered
## identically at the source. datatype=="normalized" only; subset!="3vs0".
.read_csv_base <- utils::read.csv
read.csv <- function(file, ...) {
  df <- .read_csv_base(file, ...)
  if (is.data.frame(df)) {
    if ("datatype" %in% names(df)) df <- df[df$datatype == "normalized", , drop = FALSE]
    if ("subset"   %in% names(df)) df <- df[df$subset != "3vs0",       , drop = FALSE]
  }
  df
}
## -------------------------------------------------------------------------

coh_lv <- c("matched_rm5years",
            "wholeset","wholeset_coding","stage1","stage1_coding","intersection","intersection_coding",
            "nonintersection","nonintersection_coding")
mod_lv <- c("glm_mad","glm_cor","glm_mad_nosexmt","rf_mad","xgb_mad","glm_natcommad")
stable_lv <- c("glm_mad","glm_cor")
dm_lv  <- c("normalized")
dn_lv  <- c("clinics","fpkm","fpkm_clinics","module20","module20_clinics")
cmp_lv <- c("123vs0","12vs0")
gene_fs <- c("fpkm","fpkm_clinics")
key <- c("set","cohort","datatype","featureset","subset","model","varcut")
with_clin <- c("clinics","fpkm_clinics","module20_clinics")
nc_cnt_by_sub <- list(
  "12vs0"  = c(tot = 91,  pos = 23, neg = 68),
  "123vs0" = c(tot = 113, pos = 45, neg = 68),
  "3vs0"   = c(tot = 90,  pos = 22, neg = 68))

## ============================ per-population builder ========================
build_pop <- function(glm_csv, rfxgb_csv, roc_glm_csv, roc_rfxgb_csv, md_base, gc_root, src_label) {
  d <- read.csv(glm_csv, stringsAsFactors = FALSE)
  stopifnot(all(d$model == "glm")); d$model <- NULL
  names(d)[names(d) == "mean_auc"] <- "mean"
  names(d)[names(d) == "sd_auc"]   <- "sd"
  names(d)[names(d) == "n_auc"]    <- "n"
  names(d)[names(d) == "fsmethod"] <- "model"
  d$model <- paste0("glm_", d$model)
  se <- d$sd / sqrt(10)
  d$lo <- pmax(0, d$mean - 1.96 * se); d$hi <- pmin(1, d$mean + 1.96 * se)
  d$tissue <- NA_character_

  rfxgb <- read.csv(rfxgb_csv, stringsAsFactors = FALSE)
  rfxgb <- rfxgb[rfxgb$model %in% c("rf","xgb"), ]
  rfxgb$model <- paste0(rfxgb$model, "_mad")
  se2 <- rfxgb$sd / sqrt(10)
  rfxgb$lo <- pmax(0, rfxgb$mean - 1.96 * se2); rfxgb$hi <- pmin(1, rfxgb$mean + 1.96 * se2)
  rfxgb$tissue <- NA_character_

  nature_com <- read.csv("nature_com_normal_varcut_series_auc.csv", stringsAsFactors = FALSE)
  nature_com$model <- paste0("glm_", nature_com$model)
  nature_com$tissue <- "Normal"

  d <- rbind(d[, c(key,"mean","sd","n","lo","hi","tissue")],
             rfxgb[, c(key,"mean","sd","n","lo","hi","tissue")],
             nature_com[, c(key,"mean","sd","n","lo","hi","tissue")])
  d <- d[!is.na(d$mean), ]
  for (col in intersect(c("mean","sd","lo","hi"), names(d))) d[[col]] <- round(d[[col]], 4)

  long_cols <- c("set","cohort","datatype","featureset","subset","model","varcut","mean","sd","n","lo","hi")
  long_json <- toJSON(d[, long_cols], dataframe = "values", auto_unbox = TRUE, na = "null")

  ## counts from this population's model_data
  rec_n <- function(x) { r <- as.numeric(as.character(x$recurrence)); c(tot=length(r), pos=sum(r==1), neg=sum(r==0)) }
  noclin <- list(); clin <- list(); nf_clinics <- list()
  for (coh in setdiff(coh_lv, "matched_rm5years")) for (sub in cmp_lv) {
    k <- paste(coh, sub, sep = "|")
    m20  <- readRDS(file.path(md_base, coh, "fpkm", paste0("module20_", sub, ".RDS")))
    m20c <- readRDS(file.path(md_base, coh, "fpkm", paste0("module20_clinics_", sub, ".RDS")))
    cl   <- readRDS(file.path(md_base, coh, "fpkm", paste0("clinics_", sub, ".RDS")))
    noclin[[k]] <- rec_n(m20); clin[[k]] <- rec_n(m20c); nf_clinics[[k]] <- ncol(cl) - 1
  }

  mk_headline <- function(models) {
    dd <- d[d$model %in% models, ]
    agg_best <- aggregate(mean ~ model + set, dd, max)
    agg_med  <- aggregate(mean ~ model + set, dd, median)
    merge(agg_best, agg_med, by = c("model","set"), suffixes = c("_best","_median"))
  }
  best_one <- function(dd, set, coh, sub, fs) {
    x <- dd[dd$set==set & dd$cohort==coh & dd$subset==sub & dd$featureset==fs, ]
    if (!nrow(x)) return(data.frame(auc=NA, model=NA, datatype=NA, varcut=NA, ci=NA, stringsAsFactors=FALSE))
    i <- which.max(x$mean)
    ci_txt <- sprintf("%.2f (%.2f-%.2f)", x$mean[i], x$lo[i], x$hi[i])
    data.frame(auc = round(x$mean[i],3), model = x$model[i], datatype = x$datatype[i],
               varcut = x$varcut[i], ci = ci_txt, stringsAsFactors = FALSE)
  }
  mk_best <- function(models) {
    dd <- d[d$model %in% models, ]; rows <- list()
    for (coh in coh_lv) for (sub in cmp_lv) for (fs in dn_lv) {
      k <- paste(coh, sub, sep = "|"); is_nc <- coh == "matched_rm5years"
      if (is_nc) {
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
        delta_current_minus_prev = round(c$auc - p$auc, 3), stringsAsFactors = FALSE)
    }
    do.call(rbind, rows)
  }
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
    fam$max_auc <- round(fam$max_auc, 3); fam$median_auc <- round(fam$median_auc, 3); fam
  }

  headline_stable <- mk_headline(stable_lv); headline_all <- mk_headline(mod_lv)
  best_stable <- mk_best(stable_lv); best_all <- mk_best(mod_lv)
  family_stable <- mk_family(stable_lv); family_all <- mk_family(mod_lv)

  ## cohorts tab
  gc1 <- function(set, coh, dt) ncol(readRDS(file.path(gc_root, set, coh, dt, "fpkm_123vs0.RDS"))) - 1
  gAf_c <- gc1("current","wholeset","fpkm");        gAn_c <- gc1("current","wholeset","normalized")
  gAf_p <- gc1("prev","wholeset","fpkm");           gAn_p <- gc1("prev","wholeset","normalized")
  gC_c  <- gc1("current","wholeset_coding","fpkm"); gC_p  <- gc1("prev","wholeset_coding","fpkm")
  sW <- noclin[["wholeset|123vs0"]]; sS <- noclin[["stage1|123vs0"]]; sI <- noclin[["intersection|123vs0"]]
  sN <- noclin[["nonintersection|123vs0"]]
  pfx <- if (src_label == "smoke") "EVER-SMOKER subset. " else ""
  cohorts_tab <- data.frame(
    cohort = setdiff(coh_lv, "matched_rm5years"),
    n_total = c(sW["tot"],sW["tot"], sS["tot"],sS["tot"], sI["tot"],sI["tot"], sN["tot"],sN["tot"]),
    n_pos   = c(sW["pos"],sW["pos"], sS["pos"],sS["pos"], sI["pos"],sI["pos"], sN["pos"],sN["pos"]),
    n_neg   = c(sW["neg"],sW["neg"], sS["neg"],sS["neg"], sI["neg"],sI["neg"], sN["neg"],sN["neg"]),
    description = c(
      sprintf("%sAll recurrence-labelled samples (no stage filter). Full gene set: current %d (fpkm) / %d (normalized) genes; prev %d / %d.", pfx, gAf_c, gAn_c, gAf_p, gAn_p),
      sprintf("%sSame samples as wholeset, restricted to protein-coding genes: current %d -> %d, prev %d -> %d genes.", pfx, gAf_c, gC_c, gAf_p, gC_p),
      sprintf("%sOnly stage IA1/IA2/IB samples; full gene set (same gene counts as wholeset).", pfx),
      sprintf("%sStage IA1/IA2/IB samples, protein-coding genes only (same gene counts as wholeset_coding).", pfx),
      sprintf("%sStage-I samples that also overlap the prior 246-sample paper cohort; full gene set.", pfx),
      sprintf("%sIntersection samples, protein-coding genes only.", pfx),
      sprintf("%sWholeset samples that do NOT overlap the prior 246-sample paper cohort (wholeset minus intersection). Full gene set.", pfx),
      sprintf("%sNonintersection samples, protein-coding genes only.", pfx)),
    stringsAsFactors = FALSE)

  ## ROC
  roc_pc <- read.csv(roc_glm_csv, stringsAsFactors = FALSE); roc_pc$model <- paste0("glm_", roc_pc$model)
  roc_nc <- read.csv("nature_com_normal_varcut_series_roc.csv", stringsAsFactors = FALSE); roc_nc$model <- paste0("glm_", roc_nc$model)
  roc_rx <- read.csv(roc_rfxgb_csv, stringsAsFactors = FALSE); roc_rx <- roc_rx[roc_rx$model %in% c("rf","xgb"), ]; roc_rx$model <- paste0(roc_rx$model, "_mad")
  roc <- rbind(roc_pc[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")],
               roc_rx[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")],
               roc_nc[, c(key,"fpr","tpr","n_pred","n_pos_pred","n_neg_pred")])
  roc <- merge(roc, d[, c(key,"mean")], by = key, all.x = TRUE)
  names(roc)[names(roc) == "mean"] <- "auc"
  roc <- roc[!is.na(roc$auc) & !is.na(roc$fpr), ]; roc$auc <- round(roc$auc, 4)
  roc <- roc[, c(key, "auc","n_pred","n_pos_pred","n_neg_pred","fpr","tpr")]

  meta <- list(
    generated = as.character(Sys.Date()), method = "byfold",
    method_label = sprintf("%s Per-fold AUC averaged over 10-fold x 100-boot repeats, fixed direction. 5-way model axis (glm_mad/glm_cor/glm_mad_nosexmt/rf_mad/xgb_mad).",
                           if (src_label=="smoke") "EVER-SMOKER SUBSET (Current+Former; never-smokers excluded)." else "FULL COHORT (all recurrence-labelled samples)."),
    source = if (src_label=="smoke") "ever-smoker subset of rnaseq_202609 (results_*_smoking.csv + rf_xgb_*_smoking.csv)" else "full rnaseq_202609 (results_fsmethod_comparison.csv + rf_xgb_auc_202607.csv)",
    population = src_label,
    n_configs = nrow(d), stable_configs = sum(d$model %in% stable_lv),
    provisional_xgb_configs = sum(!d$model %in% stable_lv),
    models = mod_lv, stable_models = stable_lv,
    cohorts = coh_lv, datatypes = dm_lv, featuresets = dn_lv, subsets = cmp_lv,
    varcuts = sort(unique(d$varcut)))

  J <- function(x, df="rows") toJSON(x, dataframe = df, auto_unbox = TRUE, na = "null")
  list(
    long = long_json, d = d,
    cohorts = J(cohorts_tab), headline_stable = J(headline_stable), headline_all = J(headline_all),
    best_stable = J(best_stable), best_all = J(best_all),
    family_stable = J(family_stable), family_all = J(family_all),
    roc = toJSON(roc, dataframe="rows", auto_unbox=TRUE, na="null", digits=NA),
    meta = toJSON(meta, auto_unbox = TRUE))
}

## ---- Full-vs-ever-smoker summary comparison table (static; shown at top of Overview) ----
## Best AUC over the stable GLMs (mad/cor) and, separately, RF/XGB, across the gene
## feature sets (fpkm, fpkm_clinics) and all varcuts, current RNA-seq set. "Best over
## configs" is optimistic (winner's-curse), stated in the caption.
mk_cmp_table <- function(df, ds) {
  cohs <- c("wholeset","stage1","intersection","nonintersection")
  tasks <- c("123vs0","12vs0")
  bestv <- function(d, co, su, models) {
    x <- d[d$set=="current" & d$cohort==co & d$subset==su & d$model %in% models & d$featureset %in% c("fpkm","fpkm_clinics"), ]
    if (!nrow(x)) return(NA_real_); round(max(x$mean), 3)
  }
  f3 <- function(x) if (is.na(x)) "&ndash;" else sprintf("%.3f", x)
  fd <- function(x) {
    if (is.na(x)) return("&ndash;")
    col <- if (x > 0.02) "#1a7f37" else if (x < -0.02) "#c0392b" else "#666"
    sprintf('<span style="color:%s;font-weight:600">%+.3f</span>', col, x)
  }
  rows <- ""
  for (co in cohs) for (su in tasks) {
    gf <- bestv(df,co,su,c("glm_mad","glm_cor")); gs <- bestv(ds,co,su,c("glm_mad","glm_cor"))
    mf <- bestv(df,co,su,c("rf_mad","xgb_mad"));   ms <- bestv(ds,co,su,c("rf_mad","xgb_mad"))
    lowev <- if (su == "3vs0") ' <span style="color:#b8860b" title="&le;12 recurrence events; winner-curse artifact">&#9888;</span>' else ''
    rows <- paste0(rows, sprintf(
      '<tr><td>%s</td><td>%s%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
      co, su, lowev, f3(gf), f3(gs), fd(gs-gf), f3(mf), f3(ms), fd(ms-mf)))
  }
  paste0(
    '<div style="overflow-x:auto"><table style="border-collapse:collapse;font-size:13px;min-width:640px">',
    '<thead>',
    '<tr><th rowspan="2" style="text-align:left;padding:5px 9px;border-bottom:2px solid #ccc">Cohort</th>',
    '<th rowspan="2" style="text-align:left;padding:5px 9px;border-bottom:2px solid #ccc">Task</th>',
    '<th colspan="3" style="padding:5px 9px;border-bottom:1px solid #ccc;background:#eef4fb">GLM (mad/cor)</th>',
    '<th colspan="3" style="padding:5px 9px;border-bottom:1px solid #ccc;background:#fdf0e8">RF / XGB</th></tr>',
    '<tr><th style="padding:4px 9px;border-bottom:2px solid #ccc">Full</th><th style="padding:4px 9px;border-bottom:2px solid #ccc">Ever-smoker</th><th style="padding:4px 9px;border-bottom:2px solid #ccc">&Delta;</th>',
    '<th style="padding:4px 9px;border-bottom:2px solid #ccc">Full</th><th style="padding:4px 9px;border-bottom:2px solid #ccc">Ever-smoker</th><th style="padding:4px 9px;border-bottom:2px solid #ccc">&Delta;</th></tr>',
    '</thead><tbody>', rows, '</tbody></table></div>',
    '<p class="sub" style="margin-top:7px">Best AUC over the stable GLMs (GLM+MAD / GLM+COR) and, separately, RF/XGB, across the gene feature sets (fpkm, fpkm_clinics) and all varcuts, <b>current</b> RNA-seq. &Delta; = ever-smoker &minus; full (green &gt;+0.02, red &lt;&minus;0.02). These are max-over-configs values, so they are optimistic (winner&rsquo;s curse); the small <b>intersection</b> cohort (n=37&ndash;43 after subsetting) is especially noisy. <b>Takeaway:</b> restricting to ever-smokers does not reveal a stronger recurrence signal &mdash; the overall grid mean AUC drifts slightly down and the well-powered cells barely move.</p>')
}

## ============================ build both populations ========================
cat("building FULL...\n")
full <- build_pop("results_fsmethod_comparison.csv", "rf_xgb_auc_202607.csv",
                  "results_roc_curves_fsmethods.csv", "rf_xgb_roc_202607.csv",
                  "data/model_data/current", "data/model_data", "full")
cat("building SMOKE...\n")
smoke <- build_pop("results_fsmethod_comparison_smoking.csv", "rf_xgb_auc_smoking.csv",
                   "results_roc_curves_fsmethods_smoking.csv", "rf_xgb_roc_smoking.csv",
                   "data/model_data_smoking/current", "data/model_data_smoking", "smoke")

## ============================ shared provenance panels ======================
gene_master_json     <- toJSON(read.csv("gene_corr_master.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
gene_membership_json <- toJSON(read.csv("gene_corr_membership.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
gene_scatter_json    <- toJSON(read.csv("gene_corr_scatter.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
corr_plots_json <- paste(readLines("corr_plots.json", warn = FALSE), collapse = "")
venn_json <- toJSON(read.csv("venn_madgenes.csv", stringsAsFactors = FALSE), dataframe = "rows", auto_unbox = TRUE, na = "null")
plotly_js <- paste(readLines("codes_20260713/plotly-2.32.0.min.js", warn = FALSE), collapse = "\n")

## ============================ assemble ======================================
tpl <- paste(readLines("codes_20260713/dashboard_template_fsmethods_cmp.html", warn = FALSE), collapse = "\n")
tpl <- gsub("getElementById('vcDatatype').value='fpkm'", "getElementById('vcDatatype').value='normalized'", tpl, fixed = TRUE)
sub1 <- function(tpl, tag, val) sub(tag, val, tpl, fixed = TRUE)
tpl <- sub1(tpl, "__CMP_TABLE__", mk_cmp_table(full$d, smoke$d))
tpl <- sub1(tpl, "__PLOTLY_JS__", plotly_js)
for (popname in c("FULL","SMOKE")) {
  p <- if (popname=="FULL") full else smoke
  tpl <- sub1(tpl, paste0("__LONG_DATA_",popname,"__"), p$long)
  tpl <- sub1(tpl, paste0("__COHORTS_TAB_",popname,"__"), p$cohorts)
  tpl <- sub1(tpl, paste0("__HEADLINE_STABLE_",popname,"__"), p$headline_stable)
  tpl <- sub1(tpl, paste0("__HEADLINE_ALL_",popname,"__"), p$headline_all)
  tpl <- sub1(tpl, paste0("__BEST_STABLE_",popname,"__"), p$best_stable)
  tpl <- sub1(tpl, paste0("__BEST_ALL_",popname,"__"), p$best_all)
  tpl <- sub1(tpl, paste0("__FAMILY_STABLE_",popname,"__"), p$family_stable)
  tpl <- sub1(tpl, paste0("__FAMILY_ALL_",popname,"__"), p$family_all)
  tpl <- sub1(tpl, paste0("__ROC_DATA_",popname,"__"), p$roc)
  tpl <- sub1(tpl, paste0("__META_",popname,"__"), p$meta)
}
tpl <- sub1(tpl, "__GENE_MASTER__", gene_master_json)
tpl <- sub1(tpl, "__GENE_MEMBERSHIP__", gene_membership_json)
tpl <- sub1(tpl, "__GENE_SCATTER__", gene_scatter_json)
tpl <- sub1(tpl, "__CORR_PLOTS__", corr_plots_json)
tpl <- sub1(tpl, "__VENN_DATA__", venn_json)
genelists_json <- paste(readLines("model_gene_lists_compact.json", warn = FALSE), collapse = "")
tpl <- sub1(tpl, "__GENELISTS__", genelists_json)
foldavg_json <- paste(readLines("roc_foldavg_compact.json", warn = FALSE), collapse = "")
tpl <- sub1(tpl, "__FOLDAVG__", foldavg_json)
insample_json <- paste(readLines("roc_insample_compact.json", warn = FALSE), collapse = "")
tpl <- sub1(tpl, "__INSAMPLE__", insample_json)
calib_json <- paste(readLines("roc_calib_compact.json", warn = FALSE), collapse = "")
tpl <- sub1(tpl, "__CALIB__", calib_json)

## sanity: none of OUR known placeholder tokens left behind (a generic
## __[A-Z_]+__ scan false-matches tokens inside injected Plotly/data payloads)
known <- c(outer(c("LONG_DATA","COHORTS_TAB","HEADLINE_STABLE","HEADLINE_ALL","BEST_STABLE",
                   "BEST_ALL","FAMILY_STABLE","FAMILY_ALL","ROC_DATA","META"),
                 c("FULL","SMOKE"), function(a,b) sprintf("__%s_%s__", a, b)),
           "__PLOTLY_JS__","__GENE_MASTER__","__GENE_MEMBERSHIP__","__GENE_SCATTER__","__CORR_PLOTS__","__VENN_DATA__","__CMP_TABLE__")
leftover <- known[vapply(known, function(t) grepl(t, tpl, fixed = TRUE), logical(1))]
if (length(leftover)) { cat("ERROR: unfilled placeholders:", paste(leftover, collapse=", "), "\n"); quit(status=1) }

out_file <- "htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking_normOnly_no3vs0.html"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
writeLines(tpl, out_file)
cat("Wrote", out_file, "\n")
cat("FULL n_configs:", jsonlite::fromJSON(full$meta)$n_configs,
    " SMOKE n_configs:", jsonlite::fromJSON(smoke$meta)$n_configs, "\n")
