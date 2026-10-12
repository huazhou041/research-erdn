# Version history & decision log

Chronological record of what we tried, whether it worked, and the path finally adopted.
Status key: ✅ adopted/worked · ⚠️ built but flawed · ❌ tried and rejected · 🔧 fix.

---

## 1. Modeling pipeline (`rnaseq_202609`)
- ✅ Clean re-clone of the prior `rnaseq_20260709` work (code + data, no results); paths repointed; self-contained + regenerable.
- ✅ Reran every arm from scratch: GLM (`mad`/`cor`/`mad_nosexmt`), RF + XGB, over 10 full cohorts + 8 ever-smoker cohorts × prev/current × 3 tasks × 5 varcuts. 1,296 model groups (res 480, res_rfxgb 240, res_smoking 384, res_rfxgb_smoking 192).
- ✅ Added the **`natcommad`** arm — nature_com paper's top-MAD gene list transferred to our data (full + smoking). Verified the gene list is order-identical to the dashboard Venn ranking.
- ✅ Reproduction check vs previous run: GLM + smoking arms essentially identical (corr ≈ 1.0); RF/XGB-full corr 0.992 (expected stochasticity).
- 🔧 Fixed a prior bug: the old `rf_xgb_auc` table was contaminated with ~1,800 stray GLM rows; this rerun's aggregation is clean (RF/XGB only).
- 🔧 One genuinely degenerate cell (`novel/3vs0 fpkm_clinics` XGB, too few events) aborted its group under `set -e`; hardened the RF/XGB wrapper to log + skip instead of abort.

## 2. Monitoring / infrastructure
- 🔧 Self-healing SLURM monitor; completion defined by **group.done counts**, not an empty queue (an empty queue could false-signal DONE). Made it **self-renewing** after the first monitor hit its 48 h wall-time mid-run.
- 🔧 Flaky VPN / login nodes: pinned `IdentitiesOnly=yes -i ~/.ssh/keys/nyu_cluster` to beat "too many authentication failures"; handled a host-key rotation.

## 3. Dashboards & features
- ✅ Rebuilt the three dashboards (full / smoking / full-vs-smoking toggle) with `GLM+ncMAD` wired into the model axis + legend.
- ✅ **nature_com top-MAD results table** after the Venn, with an own-cohort-MAD comparison column and CSV download.
- ✅ **Per-model gene-list download** in the explorer — GLM coefficients / RF-XGB importances, extracted from all 29,159 models.

## 4. Operating-point / optimal-cutoff panel (the main methodological iteration)
- ⚠️ **v1 — Pooled probability-cutoff slider + Youden.** Worked, but the drawn curve sat **below the diagonal** for small-event configs while the headline AUC was high.
  - Diagnosis: the headline AUC is the **mean of per-fold AUCs**; the drawn curve is the **pooled** CV predictions under one global threshold. They diverge when cross-fold probability calibration is poor (small/imbalanced cohorts).
- ❌ **Threshold-averaging** (sens/spec at a fixed probability cutoff per fold, then averaged). Rejected — mathematically ≈ pooling (differs only by fold weighting), so it does **not** fix the calibration/below-diagonal problem. (Compute job cancelled.)
- ✅ **Vertical averaging** (each fold's TPR at fixed FPR, averaged). Honest trade-off, stays above the diagonal, area matches the per-fold-averaged AUC — **but** gives no single probability cutoff (operating point parameterised by specificity).
- ⚠️ **Full-data in-sample basis** (whole model scored on its own training data). Added to expose the **deployable probability scale** — but it's **optimistic** (memorisation; e.g. `novel/3vs0` in-sample AUC = 1.0, specificity pinned at 1.0). RF/XGB in-sample predictions wouldn't reconstruct cleanly → **GLM-only** in this basis.
- 🔑 **Calibration question resolved:** calibration is fit on the **out-of-fold CV predictions** (honest), NOT on the whole-set model's own predictions (that would be the in-sample memorisation trap); it is then **applied to the full-data model** at deployment.
- ✅ **FINAL — "Calibrated CV (deployable)" basis** (default): isotonic calibration on the OOF CV predictions → honest CV sensitivity/specificity on a **meaningful probability axis** → a deployable Youden cutoff. The other three bases are kept as diagnostics.
- 🔧 UI fixes: legend/x-axis overlap; readout shows **both** sensitivity and specificity explicitly.

## 5. Deployment
- ✅ Full-data ("final") models already exist per config (the saved `.rds` is the fit on all data).
- ✅ `model_operating_points_202609.csv`: per-config AUC, pooled Youden cutoff + sens/spec, vertical-averaged optimum, **event counts**, and a **`low_event_flag`** (recur or non-recur < 10).
- ✅ `predict_new_cohort.R`: scores a new cohort with the full model, **fits isotonic calibration on the config's OOF predictions**, applies it, and classifies at a chosen calibrated cutoff.
- ✅ README **Deployment recipe** box.

### ✅ Final adopted deployment flow
1. Pick a **well-powered, non-low-event, gene-based** config.
2. Read the **Calibrated-CV** Youden cutoff from the dashboard slider.
3. `predict_new_cohort.R <config_key> <new_expr.csv> <calibrated_cutoff> <out.csv>` — full model + OOF calibration.
4. **Validate on an independent cohort** before any real use.

## 6. Version control / GitHub
- ✅ Local **code-only** git repo; `.gitignore` excludes all `data/`, `raw_data/`, `res*/`, `*.rds`, `*.csv`, built HTML. Verified **0 PHI across history** (128 files, ~5.5 MB).
- 🔧 No `gh`/token available and GitHub port 22 blocked → used `ssh.github.com:443`.
- ❌→🔧 First key test authenticated as the wrong account (`huazhou2`) due to **ssh-agent contamination**; fixed with `IdentityAgent none` on a dedicated `github-huazhou041` host alias → authenticates as **huazhou041**.
- ✅ Published via `git bundle` → local clone → `git@github-huazhou041:huazhou041/research-erdn.git` (`master` → `main`).

## 7. Key scientific takeaways
- Recurrence is only **weakly predictable** from expression here: best well-powered AUC ≈ **0.6**.
- `natcommad` (nature_com genes) is **competitive-to-neutral**, no broad lift — slightly below fold-local MAD on average.
- Small-event cells (`3vs0`, `novel`) produce **winner's-curse artifacts**: inflated AUC, deceptively narrow CV CIs (resampling reproducibility, not true uncertainty), in-sample memorisation. Flagged and **not deployable**.
