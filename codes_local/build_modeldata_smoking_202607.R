rm(list = ls())

# ---------------------------------------------------------------------------
# Ever-smoker subset of every cohort, generated DIRECTLY from raw metadata +
# raw expression (no intermediate id files, no row-subsetting of prebuilt RDS).
#
# "Ever-smoker" = raw metadata smoke_status %in% {Current, Former} (never-smokers
# dropped). NOTE this is the 3-level clinical field, NOT the pipeline's binary
# Tissue_smoke_status (which is Current-vs-rest); ever/never is the biologically
# meaningful split.
#
# Output: data/model_data_smoking/{prev,current}/{cohort}/{fpkm,normalized}/
#         {clinics,module20,module20_clinics,fpkm,fpkm_clinics}_{123vs0,12vs0,3vs0}.RDS
# Cohorts: wholeset, stage1, intersection, nonintersection (+ _coding) -- the 8 used
#          by the smoking analysis. `novel` is included in the membership block and
#          can be added to `cohorts` below if a smoking `novel` set is wanted.
#
# make_model_sets / compute_mod20 / coding_filter / task splits are byte-identical
# to build_modeldata_prev_current/nonintersection/novel so the ever-smoker subset is
# directly comparable to the full-cohort data (same features, only rows differ).
# ---------------------------------------------------------------------------

base_dir <- normalizePath(".", mustWork = TRUE)
out_root <- file.path(base_dir, "data", "model_data_smoking")

# Raw inputs live under ./raw_data so this folder is self-contained on the cluster.
# raw_data/{prev,current}/counts.{fpkm,norm,vst}.csv + metadata/module/gene-type files.
# Override with ERDN_RAW to point elsewhere (e.g. the /Volumes deseq dirs on the Mac).
RAW <- Sys.getenv("ERDN_RAW", file.path(base_dir, "raw_data"))
metadata_file <- file.path(RAW, "metadata_model_202606.csv")
gene_type_file <- file.path(RAW, "gene_type.csv")
module_gene_file <- file.path(RAW, "246samples.coding10000.modules.dim10.k20.genes.csv")
old_module_score_file <- file.path(RAW, "246samples.coding10000.modules.dim10.k20.scores.csv")

deseq_sets <- list(
  prev = file.path(RAW, "prev"),
  current = file.path(RAW, "current")
)

read_expr <- function(path) {
  expr <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(expr)[1] <- "gene"; expr
}
drop_missing_rows <- function(x) x[stats::complete.cases(x), , drop = FALSE]
extract_numeric_id <- function(x) suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))

make_model_sets <- function(expr, clinics, mod20, outdir, restrict_samples = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  expr_t <- data.frame(t(expr[, -1, drop = FALSE]), check.names = FALSE)
  names(expr_t) <- expr$gene
  expr_t$Study_ID <- rownames(expr_t)

  common_samples <- intersect(clinics$Study_ID, expr_t$Study_ID)
  if (!is.null(restrict_samples)) common_samples <- intersect(restrict_samples, common_samples)
  temp <- data.frame(
    clinics[match(common_samples, clinics$Study_ID), ],
    Tissue_Module20 = mod20[match(common_samples, rownames(mod20)), 1],
    expr_t[match(common_samples, expr_t$Study_ID), setdiff(names(expr_t), "Study_ID"), drop = FALSE],
    check.names = FALSE
  )
  rownames(temp) <- common_samples

  saveRDS(drop_missing_rows(temp[, c("recurrence", "Tissue_Module20")]),
          file.path(outdir, "module20_123vs0.RDS"))
  saveRDS(drop_missing_rows(temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                     "Tissue_smoke_status", "Tissue_Tumor_size")]),
          file.path(outdir, "clinics_123vs0.RDS"))
  saveRDS(drop_missing_rows(temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                     "Tissue_smoke_status", "Tissue_Tumor_size", "Tissue_Module20")]),
          file.path(outdir, "module20_clinics_123vs0.RDS"))
  saveRDS(drop_missing_rows(temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                     "Tissue_smoke_status", "Tissue_Tumor_size", expr$gene)]),
          file.path(outdir, "fpkm_clinics_123vs0.RDS"))
  saveRDS(drop_missing_rows(temp[, c("recurrence", expr$gene)]),
          file.path(outdir, "fpkm_123vs0.RDS"))

  task_defs <- list("12vs0" = c(0, 1, 2), "3vs0" = c(0, 3))
  for (task in names(task_defs)) {
    task_temp <- temp[temp$ProgType %in% task_defs[[task]], , drop = FALSE]
    saveRDS(drop_missing_rows(task_temp[, c("recurrence", "Tissue_Module20")]),
            file.path(outdir, paste0("module20_", task, ".RDS")))
    saveRDS(drop_missing_rows(task_temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                            "Tissue_smoke_status", "Tissue_Tumor_size")]),
            file.path(outdir, paste0("clinics_", task, ".RDS")))
    saveRDS(drop_missing_rows(task_temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                            "Tissue_smoke_status", "Tissue_Tumor_size", "Tissue_Module20")]),
            file.path(outdir, paste0("module20_clinics_", task, ".RDS")))
    saveRDS(drop_missing_rows(task_temp[, c("recurrence", "Tissue_Sex", "Tissue_Age", "Tissue_race",
                                            "Tissue_smoke_status", "Tissue_Tumor_size", expr$gene)]),
            file.path(outdir, paste0("fpkm_clinics_", task, ".RDS")))
    saveRDS(drop_missing_rows(task_temp[, c("recurrence", expr$gene)]),
            file.path(outdir, paste0("fpkm_", task, ".RDS")))
  }
}

# --- metadata + cohort memberships (all defined inline from raw metadata) -----
metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
clinic_cols <- c("Study_ID", "recurrence", "TTP", "ProgType", "Tissue_Sex", "Tissue_Age",
                 "Tissue_race", "Tissue_smoke_status", "Tissue_Tumor_size")
clinics_whole <- metadata[!is.na(metadata$recurrence), clinic_cols]
rownames(clinics_whole) <- clinics_whole$Study_ID
clinics_stage1 <- metadata[metadata$stage_lungca %in% c("IA1", "IA2", "IB") & !is.na(metadata$recurrence), clinic_cols]
rownames(clinics_stage1) <- clinics_stage1$Study_ID

# EVER-SMOKER membership, inline from the raw 3-level smoke_status
ever_samples <- metadata$Study_ID[metadata$smoke_status %in% c("Current", "Former") & !is.na(metadata$recurrence)]

# intersection = stage-1 AND old-NatComm normal-tissue RNA-seq IDs
old_scores <- read.csv(old_module_score_file, stringsAsFactors = FALSE, check.names = FALSE)
old_ids <- extract_numeric_id(names(old_scores)[grep("-N", names(old_scores))])
curr_s1 <- extract_numeric_id(clinics_stage1$Study_ID)
intersection_samples <- clinics_stage1$Study_ID[match(intersect(old_ids, curr_s1), curr_s1)]
intersection_samples <- intersection_samples[!is.na(intersection_samples)]
nonintersection_samples <- setdiff(clinics_whole$Study_ID, intersection_samples)

module_genes <- read.csv(module_gene_file, check.names = FALSE, stringsAsFactors = FALSE)
module20_genes <- module_genes[module_genes[, 2] == 20, 1]
gene_type_map <- read.csv(gene_type_file, check.names = FALSE, stringsAsFactors = FALSE)

compute_mod20 <- function(deseq_dir) {
  vst <- read_expr(file.path(deseq_dir, "counts.vst.csv"))
  vst_mat <- as.matrix(vst[, -1, drop = FALSE]); rownames(vst_mat) <- vst$gene
  scaled_expr <- t(apply(vst_mat, 1, scale))
  mod20 <- data.frame(module20 = colMeans(scaled_expr[rownames(scaled_expr) %in% module20_genes, , drop = FALSE], na.rm = TRUE))
  rownames(mod20) <- colnames(vst_mat); mod20
}
coding_filter <- function(expr) {
  gt <- gene_type_map$gene_type[match(expr$gene, gene_type_map$gene_name)]
  expr[!is.na(gt) & gt == "protein_coding", , drop = FALSE]
}

# cohort -> (clinics table, membership restricted to ever-smokers)
cohorts <- list(
  wholeset        = list(clin = clinics_whole,  restrict = ever_samples),
  stage1          = list(clin = clinics_stage1, restrict = ever_samples),
  intersection    = list(clin = clinics_stage1, restrict = intersect(intersection_samples, ever_samples)),
  nonintersection = list(clin = clinics_whole,  restrict = intersect(nonintersection_samples, ever_samples))
)

cat("ever-smokers (labeled):", length(ever_samples), "\n")
summary_rows <- list()
for (set_name in names(deseq_sets)) {
  deseq_dir <- deseq_sets[[set_name]]
  cat("\n==== RNA-seq set:", set_name, "====\n")
  mod20 <- compute_mod20(deseq_dir)
  expr_files <- list(fpkm = file.path(deseq_dir, "counts.fpkm.csv"),
                     normalized = file.path(deseq_dir, "counts.norm.csv"))
  for (data_type in names(expr_files)) {
    expr <- read_expr(expr_files[[data_type]]); expr_coding <- coding_filter(expr)
    for (coh in names(cohorts)) {
      cc <- cohorts[[coh]]
      make_model_sets(expr,        cc$clin, mod20, file.path(out_root, set_name, coh, data_type),            cc$restrict)
      make_model_sets(expr_coding, cc$clin, mod20, file.path(out_root, set_name, paste0(coh, "_coding"), data_type), cc$restrict)
    }
    for (coh in unlist(lapply(names(cohorts), function(c) c(c, paste0(c, "_coding"))))) {
      for (task in c("123vs0", "12vs0", "3vs0")) {
        d <- readRDS(file.path(out_root, set_name, coh, data_type, paste0("fpkm_", task, ".RDS")))
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          set = set_name, cohort = coh, data = task, type = data_type,
          No.Recur = sum(d$recurrence == 1), No.NonRecu = sum(d$recurrence == 0), Total = nrow(d))
      }
    }
    cat(sprintf("  %s / %s: done\n", set_name, data_type))
  }
}
summary_table <- do.call(rbind, summary_rows)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(out_root, "counts_summary_smoking.csv"), row.names = FALSE)
cat("\nWrote", file.path(out_root, "counts_summary_smoking.csv"), "\n")
print(summary_table[summary_table$set == "current" & summary_table$type == "fpkm", ])
