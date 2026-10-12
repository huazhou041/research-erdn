#!/usr/bin/env Rscript
# Cox (glmnet family="cox") counterpart of tune_drug.R.
#
# Differences from the binary pipeline that matter:
#  * outcome is (TTP, event) instead of binary recurrence  -> uses the follow-up time
#  * performance = Harrell's C (headline) + Uno's C (IPCW, censoring-robust)
#  * C is computed on POOLED out-of-fold risk scores within each repeat (every subject
#    predicted exactly once per repeat), NOT by averaging per-fold C. Averaging per-fold
#    metrics on ~15-sample folds is what inflated the binary 3vs0 AUCs.
#  * risk-score orientation is FIXED (reverse=TRUE: higher score = higher hazard).
#    Never let the metric pick the direction that flatters the model.
suppressPackageStartupMessages({
  library(argparser); library(glmnet); library(survival); library(readr); library(parallel)
})
p = arg_parser("Cox training pipeline (glmnet family=cox)")
p = add_argument(p, "matrix",  help = "feature RDS (rownames = Study_ID)")
p = add_argument(p, "surv",    help = "survival csv: Study_ID,TTP,event_all,event_recur,ProgType")
p = add_argument(p, "prefix",  help = "output prefix")
p = add_argument(p, "--endpoint", default = "recur", help = "recur (systemic+locoregional; 2nd primary censored) | all")
p = add_argument(p, "--varcut",   default = 100,  help = "top-N genes by fold-local ranking")
p = add_argument(p, "--kfold",    default = 10L,  help = "CV folds")
p = add_argument(p, "--nboot",    default = 100L, help = "repeats of k-fold CV")
p = add_argument(p, "--fs_method",default = "mad", help = "mad | cor (cor = |Cox score| via univariate concordance)")
p = add_argument(p, "--alpha",    default = 0.5,  help = "elastic-net alpha")
argv = parse_args(p)

mat_file <- argv$matrix; surv_file <- argv$surv; prefix <- argv$prefix
endpoint <- argv$endpoint; varcut <- as.integer(argv$varcut)
kfold <- as.integer(argv$kfold); kboot <- as.integer(argv$nboot)
fs_method <- argv$fs_method; alpha <- as.numeric(argv$alpha)
stopifnot(endpoint %in% c("recur","all"), fs_method %in% c("mad","cor"))

mat <- readRDS(mat_file)
sv  <- read.csv(surv_file, stringsAsFactors = FALSE)
rownames(sv) <- sv$Study_ID
common <- intersect(rownames(mat), rownames(sv))
if (length(common) < 20) stop("too few samples joined to survival table: ", length(common))
mat <- mat[common, , drop = FALSE]
sv  <- sv[common, , drop = FALSE]

time   <- as.numeric(sv$TTP)
status <- as.numeric(if (endpoint == "recur") sv$event_recur else sv$event_all)
ok <- is.finite(time) & time > 0 & is.finite(status)
mat <- mat[ok, , drop = FALSE]; time <- time[ok]; status <- status[ok]
feats <- setdiff(colnames(mat), "recurrence")
X <- as.matrix(mat[, feats, drop = FALSE])
storage.mode(X) <- "numeric"
# drop zero-variance / all-NA features up front
keep <- apply(X, 2, function(g) { s <- sd(g, na.rm = TRUE); is.finite(s) && s > 0 })
X <- X[, keep, drop = FALSE]
n <- nrow(X)
message(sprintf("%s | endpoint=%s n=%d events=%d features=%d varcut=%d",
                basename(mat_file), endpoint, n, sum(status), ncol(X), varcut))
if (sum(status) < 10) { message("fewer than 10 events -- skipping"); quit(save="no", status=0) }

# fold-local feature ranking (never sees the held-out fold)
rank_feats <- function(Xtr, ytime, ystat, top) {
  if (ncol(Xtr) <= top) return(colnames(Xtr))
  sc <- if (fs_method == "cor") {
    # supervised: univariate Cox score = |concordance - 0.5|, fold-local
    apply(Xtr, 2, function(g) {
      if (sd(g, na.rm = TRUE) == 0) return(0)
      cc <- tryCatch(survival::concordance(survival::Surv(ytime, ystat) ~ g, reverse = TRUE)$concordance,
                     error = function(e) 0.5)
      abs(cc - 0.5)
    })
  } else {
    apply(Xtr, 2, mad, na.rm = TRUE)
  }
  sc[!is.finite(sc)] <- 0
  colnames(Xtr)[order(sc, decreasing = TRUE)[seq_len(top)]]
}

# Harrell's + Uno's C for a RISK score (higher = higher hazard). Orientation fixed.
c_index <- function(ytime, ystat, risk, uno = FALSE) {
  if (length(unique(ystat)) < 2 || sd(risk) == 0) return(NA_real_)
  tryCatch({
    if (uno) survival::concordance(survival::Surv(ytime, ystat) ~ risk,
                                   reverse = TRUE, timewt = "n/G2")$concordance
    else     survival::concordance(survival::Surv(ytime, ystat) ~ risk,
                                   reverse = TRUE)$concordance
  }, error = function(e) NA_real_)
}

set.seed(123)
per_repeat <- vector("list", kboot)
oof_all    <- vector("list", kboot)
for (b in seq_len(kboot)) {
  # stratify folds on event status so no fold is event-free
  flds <- caret::createFolds(factor(status), k = kfold, list = TRUE, returnTrain = FALSE)
  risk <- rep(NA_real_, n)
  for (i in seq_along(flds)) {
    te <- flds[[i]]; tr <- setdiff(seq_len(n), te)
    if (sum(status[tr]) < 3) next
    sel <- rank_feats(X[tr, , drop = FALSE], time[tr], status[tr], varcut)
    Xtr <- X[tr, sel, drop = FALSE]; Xte <- X[te, sel, drop = FALSE]
    ytr <- survival::Surv(time[tr], status[tr])
    fit <- tryCatch(
      cv.glmnet(Xtr, ytr, family = "cox", alpha = alpha, nfolds = 5),
      error = function(e) NULL)
    if (is.null(fit)) next
    lp <- tryCatch(as.numeric(predict(fit, newx = Xte, s = "lambda.min", type = "link")),
                   error = function(e) rep(NA_real_, length(te)))
    risk[te] <- lp
  }
  good <- is.finite(risk)
  if (sum(good) < 20 || length(unique(status[good])) < 2) next
  per_repeat[[b]] <- data.frame(
    bootid = b,
    C_harrell = c_index(time[good], status[good], risk[good], uno = FALSE),
    C_uno     = c_index(time[good], status[good], risk[good], uno = TRUE),
    n_scored  = sum(good))
  oof_all[[b]] <- data.frame(subj = rownames(X)[good], bootid = b, risk = risk[good],
                             time = time[good], status = status[good])
}
per_repeat <- do.call(rbind, per_repeat)
oof_all    <- do.call(rbind, oof_all)
if (is.null(per_repeat) || !nrow(per_repeat)) { message("no usable repeats"); quit(save="no", status=0) }

summary_df <- data.frame(
  matrix = basename(mat_file), endpoint = endpoint, fs_method = fs_method, varcut = varcut,
  n = n, events = sum(status), n_features = ncol(X), repeats = nrow(per_repeat),
  C_harrell_mean = mean(per_repeat$C_harrell, na.rm = TRUE),
  C_harrell_sd   = sd(per_repeat$C_harrell,  na.rm = TRUE),
  C_uno_mean     = mean(per_repeat$C_uno,    na.rm = TRUE),
  C_uno_sd       = sd(per_repeat$C_uno,      na.rm = TRUE))
write_csv(per_repeat, paste0(prefix, "_cvSummary.csv"))
write_csv(oof_all,    paste0(prefix, "_oof_risk.csv"))
write_csv(summary_df, paste0(prefix, "_summary.csv"))
message(sprintf("  -> C(Harrell)=%.4f +/- %.4f   C(Uno)=%.4f   [%d repeats]",
                summary_df$C_harrell_mean, summary_df$C_harrell_sd,
                summary_df$C_uno_mean, nrow(per_repeat)))
