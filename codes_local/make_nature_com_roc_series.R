# Build a full varcut x subset series for nature_com Normal tissue
# (matched_rm5years/normalized/fpkm/glm), matching the same fixed-direction AUC
# methodology and pooled-ROC-curve construction used for prev/current, so
# nature_com becomes a genuinely comparable unit across the Heatmap/Varcut/ROC
# tabs, not just a single floating point.
#
# Each subset (12vs0/123vs0/3vs0) needs its OWN dedicated label file -- they are
# NOT interchangeable subject sets. Found at
# nested_cv/data_20220524/10fold/rdata/normalized/matched_rm5years/control/
# control_fpkm_<subset>.RDS on the cluster (the actual training-time label
# source, one level more "canonical" than the models_202205 copy used earlier
# for 12vs0 -- confirmed to give byte-identical results: 91 subjects, 68 neg/23
# pos, matching the original Section 2 numbers exactly).
#
# AUC: per-(foldid,bootid) fixed-direction (direction="<"), mean/sd/n over up to
#   1000 fold-level values -- identical to Section 2's methodology.
# ROC curve: pooled across ALL fold x boot predictions, threshold sweep
#   seq(1,0,by=-0.01) + Inf/-Inf, matching res/make_roc_curves_202607.R exactly
#   (same thresholds, same pipe-delimited fpr/tpr string format) so this plugs
#   straight into the same dashboard ROC_DATA schema/JS parser.
rm(list = ls())
suppressPackageStartupMessages(library(pROC))

base <- Sys.getenv("NATURECOM_DIR", "nature_com")
preds_dir <- file.path(base, "revise_202302/test/res_withmodel")
label_dir <- file.path(base, "rdata_labels")

thresholds <- c(Inf, seq(1, 0, by = -0.01), -Inf)
fmt_vec <- function(x) paste(sprintf("%.3f", x), collapse = "|")

calc_curve <- function(y, prob) {
  ok <- !is.na(y) & !is.na(prob)
  y <- y[ok]; prob <- prob[ok]
  pos <- sum(y == 1); neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(list(fpr = NA_character_, tpr = NA_character_, n_pred = length(y), n_pos_pred = pos, n_neg_pred = neg))
  tpr <- numeric(length(thresholds)); fpr <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    pred <- prob >= thresholds[i]
    tpr[i] <- sum(pred & y == 1) / pos
    fpr[i] <- sum(pred & y == 0) / neg
  }
  list(fpr = fmt_vec(fpr), tpr = fmt_vec(tpr), n_pred = length(y), n_pos_pred = pos, n_neg_pred = neg)
}

fixed_auc <- function(y, prob) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(y, prob, quiet = TRUE, direction = "<"))
}

varcuts <- c(50, 100, 200, 500, 1000)
subsets <- c("12vs0", "123vs0", "3vs0")
auc_rows <- list(); roc_rows <- list()

for (sub in subsets) {
  rds <- readRDS(file.path(label_dir, sprintf("control_fpkm_%s.RDS", sub)))
  y_map <- setNames(as.numeric(as.character(rds$recurrence)), rownames(rds))
  cat(sprintf("\n=== subset %s: %d labeled subjects, %d pos / %d neg ===\n",
              sub, length(y_map), sum(y_map == 1), sum(y_map == 0)))

  for (vc in varcuts) {
    pf <- file.path(preds_dir, sprintf("normalized_matched_rm5years_control_fpkm_%s_glm_v%d_pred_prob.csv", sub, vc))
    pred <- read.csv(pf, stringsAsFactors = FALSE)
    # y_map[pred$subj] would silently do POSITIONAL indexing (pred$subj is integer) instead
    # of name-based lookup, returning near-garbage -- must coerce to character first.
    pred$y <- y_map[as.character(pred$subj)]
    pred <- pred[!is.na(pred$y), ]

    key <- paste(sub, vc)
    grp <- split(pred, list(pred$foldid, pred$bootid), drop = TRUE)
    aucs <- vapply(grp, function(g) fixed_auc(g$y, g$prob), numeric(1))
    aucs <- aucs[is.finite(aucs)]
    auc_rows[[key]] <- data.frame(
      set = "nature_com", cohort = "matched_rm5years", datatype = "normalized",
      featureset = "fpkm", subset = sub, model = "mad", varcut = vc,
      mean = mean(aucs), sd = sd(aucs), n = length(aucs), tissue = "Normal",
      stringsAsFactors = FALSE)

    curve <- calc_curve(pred$y, pred$prob)
    roc_rows[[key]] <- data.frame(
      set = "nature_com", cohort = "matched_rm5years", datatype = "normalized",
      featureset = "fpkm", subset = sub, model = "mad", varcut = vc,
      fpr = curve$fpr, tpr = curve$tpr,
      n_pred = curve$n_pred, n_pos_pred = curve$n_pos_pred, n_neg_pred = curve$n_neg_pred,
      stringsAsFactors = FALSE)

    cat(sprintf("  v%d: mean AUC=%.4f (n=%d folds), pooled n=%d (%d pos/%d neg)\n",
                vc, mean(aucs), length(aucs), curve$n_pred, curve$n_pos_pred, curve$n_neg_pred))
  }
}

## ---- 3vs0/v200 override: published figure value, from a frozen single-bootid snapshot ----
## The 100-boot fold-averaged fixed-direction AUC computed above for 3vs0 is flat at
## chance (~0.47-0.49) at every varcut -- but the actual published manuscript figure
## (revise_202302/reviewer4/roc_control_3vs0.png, from review4_roc_3vs0_202308.R)
## reports AUC = 0.71 (0.61~0.81) for this exact config (control/Normal, 3vs0, v200).
## That script uses a different, narrower methodology: a SINGLE bootstrap replicate
## (bootid==1 only, not the 100-boot average), pooled across its 10 folds.
##
## Traced and confirmed genuine (see nature_com_summary_src.md Section 11.15):
## revise_202302/ctrlroc_3vs0.csv is a frozen snapshot of exactly that bootid==1
## pooled prediction set (90 subjects, 22 pos/68 neg, matching the canonical 3vs0
## label set used above). Recomputing AUC directly from it reproduces the published
## number almost exactly (0.6972 vs. published 0.71/0.70). It is NOT reproducible
## from the *current* pred_prob file for this config any more -- that file has since
## been regenerated (different seed/fold split entirely: 0/90 probabilities match by
## subject) and its bootid==1 draw happens to be a degenerate, collapsed-to-intercept
## GLM fit (54/90 subjects share one constant probability) -- hence the ~0.49 chance
## reading computed above from that file is an honest description of a DIFFERENT
## (degenerate) model fit, not a refutation of the published 0.70.
##
## Per instruction, this frozen-snapshot value replaces the standard loop's 3vs0/v200
## row below (only v200 -- the other 4 varcuts have no equivalent published snapshot
## and are left as the standard 100-boot computation). CI uses the published figure's
## own reported bounds (0.61-0.81) rather than a freshly recomputed DeLong CI (which
## gives a wider 0.55-0.85 on this same frozen n=90 snapshot) so the dashboard shows
## exactly what was published.
snap <- read.csv(file.path(base, "revise_202302", "ctrlroc_3vs0.csv"), stringsAsFactors = FALSE)
snap_auc <- as.numeric(pROC::auc(snap$prog, snap$prob, quiet = TRUE))
cat(sprintf("\n=== 3vs0/v200 override: published snapshot AUC=%.4f (n=%d, %d pos/%d neg) ===\n",
            snap_auc, nrow(snap), sum(snap$prog==1), sum(snap$prog==0)))
auc_rows[["3vs0 200"]] <- data.frame(
  set = "nature_com", cohort = "matched_rm5years", datatype = "normalized",
  featureset = "fpkm", subset = "3vs0", model = "mad", varcut = 200,
  mean = round(snap_auc, 4), sd = NA_real_, n = nrow(snap), tissue = "Normal",
  stringsAsFactors = FALSE)
snap_curve <- calc_curve(snap$prog, snap$prob)
roc_rows[["3vs0 200"]] <- data.frame(
  set = "nature_com", cohort = "matched_rm5years", datatype = "normalized",
  featureset = "fpkm", subset = "3vs0", model = "mad", varcut = 200,
  fpr = snap_curve$fpr, tpr = snap_curve$tpr,
  n_pred = snap_curve$n_pred, n_pos_pred = snap_curve$n_pos_pred, n_neg_pred = snap_curve$n_neg_pred,
  stringsAsFactors = FALSE)

auc_tab <- do.call(rbind, auc_rows)
se <- auc_tab$sd / sqrt(10)
auc_tab$lo <- pmax(0, auc_tab$mean - 1.96 * se); auc_tab$hi <- pmin(1, auc_tab$mean + 1.96 * se)
## 3vs0/v200 override row has sd=NA (single-snapshot point estimate, not a
## fold-level spread) -- use the published figure's own reported CI bounds
## directly instead of the NA that the se-based formula above would produce.
is_override <- auc_tab$subset == "3vs0" & auc_tab$varcut == 200
auc_tab$lo[is_override] <- 0.61; auc_tab$hi[is_override] <- 0.81
roc_tab <- do.call(rbind, roc_rows)

out_dir <- Sys.getenv("GENECORR_OUT", "htmls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(auc_tab, file.path(out_dir, "nature_com_normal_varcut_series_auc.csv"), row.names = FALSE)
write.csv(roc_tab, file.path(out_dir, "nature_com_normal_varcut_series_roc.csv"), row.names = FALSE)
cat("\nwrote nature_com_normal_varcut_series_auc.csv and _roc.csv (", nrow(auc_tab), "rows )\n")
print(auc_tab[order(auc_tab$subset, auc_tab$varcut), c("subset","varcut","mean","sd","n","lo","hi")])
