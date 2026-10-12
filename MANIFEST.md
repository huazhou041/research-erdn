# rnaseq_202609 — corrected, self-contained pipeline (code + data + raw_data)

Clone of rnaseq_20260709: code + `data/` (no `res*/`/`results/`/`htmls/`/output CSVs),
plus `raw_data/` so model data regenerates on the cluster. Paths repointed to the folder;
log dirs pre-created.

## Cohort glossary (verified vs raw metadata; built RDS set-equal to defs)
- wholeset (150): all recurrence-labelled. stage1 (120): IA1/IA2/IB.
- intersection (58): stage-1 AND old-NatComm normal ("-N") RNA-seq (narrow).
- nonintersection (92): wholeset - intersection; still 18 NatComm inside (9 not stage-1;
  9 stage-1 NatComm with no old NatComm RNA-seq). NOT NatComm-free.
- novel (74): wholeset - full NatComm(76 = old-"-N" union "Nature"-flag). Clean/independent.
- ever-smoker subset: cohort restricted to smoke_status in {Current,Former} (150->120);
  3-level clinical field, not binary Tissue_smoke_status.
- tasks: 123vs0=all; 12vs0=ProgType{0,1,2}; 3vs0=ProgType{0,3}. recurrence=1 iff ProgType in {1,2,3}.

## nature_com alignment (verified in code)
- Fold-local unsupervised MAD: FS_top(mat[fold_id != i, ], top) (train folds only).
- Fixed direction = "<" AUC. 100x10 fold-averaged (GLM nboot 100; RF/XGB nboot 50; k=10).

## Corrections applied
1. Aggregators (GLM & RF/XGB) include novel_coding|novel; aggregate_rfxgb.R writes the
   complete rf_xgb_auc_202607.csv directly (no merge; merge script removed).
2. tune_drug.R errors on --learner DL/RF (use tune_drug_rf.R / _dl.R / _xgb.R).
3. mlp_pipeline.py global top-3000 MAD pre-filter removed -> strictly fold-local.
4. Provenance generators use env-overridable paths.
5. Submit defaults reproduce the FULL run: GLM = 10 cohorts x {mad,cor,mad_nosexmt} (9000);
   RF/XGB = 10 cohorts.
6. Ever-smoker data generated directly from raw metadata by
   codes_local/build_modeldata_smoking_202607.R (inline ever = Current+Former; no id files).
7. Provenance: raw_data/ holds all raw inputs; data/model_data_smoking regenerated on the
   cluster with the new builder, verified value-identical to prior (96/96, 480 RDS, 806M).
8. NO INTERMEDIATE exclude_sex_mt_genes.txt. fs_method=mad_nosexmt now derives the exclusion
   INLINE in tune_drug.R: sex-chromosome genes from the annotation chr column (chrX/chrY) +
   an embedded 86-gene mitochondrial list (annotation has no chrM). Reproduces the former
   exclude_sex_mt_genes.txt exactly (3055 genes, set-equal). GLM wrappers pass
   --exclude_sexmt_annot raw_data/gene_type.csv.
9. ALL FOUR model-data builders (prev_current, nonintersection, novel, smoking) now read
   raw_data/ (0 /Volumes refs) -> the full model data is cluster-regenerable. compare_smoking.R
   baseline repointed to the source folder's Aug-10 data (re-runnable). aggregate_fsmethods.R
   stale "216-job/5400-config" comment corrected to the 10-cohort/3-fs (~9000) grid.

Auxiliary QC/plotting scripts in codes_local/ (make_gene_correlation_data.R, make_venn_data.R,
build_corr_plots_json.py, make_nature_com_roc_series.R) now contain NO /Volumes paths: inputs that
exist here point at raw_data/ (gene_type.csv, prev|current/counts.norm.csv), ERDN_PROJ defaults to ".",
and inputs NOT shipped in this folder are fail-loud env vars -- GENECORR_RAW (per-cohort expr subsets),
CORR_PNG_DIR (correlation PNGs), NATURECOM_NORM_FILE / NATURECOM_123_RDS / NATURECOM_DIR (old nature_com
expression + prediction files). Set those env vars to run these OPTIONAL panels; they build gene-correlation /
venn / sample-correlation QC, not the model data or AUC pipeline, so they are not needed to reproduce results.

NOTE: this folder ships code + data only -- no training outputs yet (res/, res_rfxgb/, res_smoking/,
res_rfxgb_smoking/ are empty). Regenerate results with the submit/aggregate/dashboard steps below.

## Layout
- root: submit drivers, aggregators, ROC builders, tune wrappers, sbatch runners
  (regen_smoking.sbatch, compare_smoking.R, finalize_rfxgb.sbatch), log dirs
- slurm/deep-drug-response/: tune_drug*.R (+ inline sex/MT exclusion) + hyperparam grids
- codes_20260713/: dashboard builders + templates + plotly
- codes_local/: build_modeldata_{prev_current,nonintersection,novel,smoking}_202607.R
  (all read raw_data/), 3 provenance generators, make_nature_com_roc_series.R
- data/: model_data (10 cohorts) + model_data_smoking (8, regenerated) + metadata/counts
- raw_data/: metadata, gene_type.csv, module genes/scores, old_counts.normalized.csv,
  NatComm mapping, prev/ & current/ counts.{fpkm,norm,vst}.csv

## Regenerate (all on the cluster from raw_data/, run from the folder root)
- Model data: `Rscript codes_local/build_modeldata_prev_current_202607.R` (+ nonintersection,
  novel, smoking). Smoking convenience: `sbatch regen_smoking.sbatch`.
- GLM: `./run_submit_20260709_glm_fsmethods.sh` (smoking: `./submit_glm_smoking.sh`).
- RF/XGB: `./submit_rfxgb.sh` (smoking: `./submit_rfxgb_smoking.sh`).
- Aggregate: `Rscript aggregate_fsmethods.R`; `sbatch finalize_rfxgb.sbatch` (+ smoking).
- ROC + dashboards: codes_20260713/make_roc_curves_fsmethods.R, make_dashboard_fsmethods*.R.

## NEW arm: natcommad (nature_com top-MAD gene list, transferred to our data)
Trains GLM on prev/current data using a FIXED external feature set = the top-N
nature_com MAD genes (no fold-local ranking, no outcome use). Pieces:
- codes_local/build_natcom_madgenes.R -> raw_data/natcom_mad_ranked.csv. Computed
  IDENTICALLY to the dashboard Venn (make_venn_data.R mad_rank): per-gene MAD of
  log2(normalized+1e-4) over nature_com NORMAL-tissue samples, restricted to the fixed
  123 (121 found) nature_com training patients, sorted descending. Verified: top-50 ==
  Venn nature_com top-50, order-exact.
- tune_drug.R fs_method=natcommad + --rank_genes_file: FS_top uses head(natcom_ranked
  intersect data-genes, varcut) for every fold (so _coding cohorts naturally get the
  coding subset). Clinical columns kept as usual.
- submit_glm_natcommad.sh + run_one_natcommad.sh: GLM natcommad over prev/current x 10
  cohorts x {fpkm,fpkm_clinics} x {123vs0,12vs0,3vs0} x varcut{50,100,200,500,1000},
  written to res/ as ..._glm_natcommad_v<vc>_cvSummary.csv.
- aggregate_fsmethods.R fsmethod_re now includes natcommad -> results_fsmethod_comparison.csv
  gains glm_natcommad rows (the transferred-nature_com-genes AUC).
- raw_data/ now also holds naturecom_counts.normalized.csv (270-sample) +
  naturecom_control_fpkm_123vs0.RDS -> the ranking (and the Venn's nature_com side) are
  reproducible on the cluster.
- TO DISPLAY in the dashboard: add "glm_natcommad" to mod_lv/labelModel/colors in
  codes_20260713/make_dashboard_fsmethods.R (AUC is produced regardless).

## Run the natcommad arm
1. `Rscript codes_local/build_natcom_madgenes.R`  (-> raw_data/natcom_mad_ranked.csv)
2. `./submit_glm_natcommad.sh`  (prev,current; add args to narrow cohorts)
3. `Rscript aggregate_fsmethods.R`  (glm_natcommad rows appear in results_fsmethod_comparison.csv)
