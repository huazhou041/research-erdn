rm(list = ls())

# ---------------------------------------------------------------------------
# Add the `nonintersection` cohort = wholeset (150) MINUS the `intersection` cohort (58)
# = 92 samples, built in exactly the same format as wholeset / stage1 / intersection:
#   data/model_data/{prev,current}/{nonintersection,nonintersection_coding}/{fpkm,normalized}/
#   {clinics,module20,module20_clinics,fpkm,fpkm_clinics}_{123vs0,12vs0,3vs0}.RDS
#
# IMPORTANT CAVEAT (recorded in the summary csv):
#   `intersection` is NOT the full Nature Communications overlap -- it is
#   (old NatComm normal RNA-seq) AND (stage-1), which is why it is only 58.
#   Full NatComm contact = union of
#     (a) has a normal-tissue "-N" column in the old NatComm expression matrix  -> 66
#     (b) flagged "Nature" in Mapping_file_surgery_data Nature communication overlap -> 75
#   union = 76.  So subtracting only the 58 leaves 76 - 58 = 18 Nature Communications
#   patients still inside this cohort. It is therefore NOT a clean NatComm-free set.
#   (A leakage-free variant, `novel` = 150 - 76 = 74, is built by
#    build_modeldata_novel_202607.R.)
#
# Everything else (make_model_sets, module20 score, coding filter, task splits) is
# identical to build_modeldata_prev_current_202607.R.
# ---------------------------------------------------------------------------

base_dir <- normalizePath(".", mustWork = TRUE)
out_root <- file.path(base_dir, "data", "model_data")

metadata_file <- "raw_data/metadata_model_202606.csv"
gene_type_file <- "raw_data/gene_type.csv"
module_gene_file <- "raw_data/246samples.coding10000.modules.dim10.k20.genes.csv"
old_module_score_file <- "raw_data/246samples.coding10000.modules.dim10.k20.scores.csv"
nc_mapping_file <- "raw_data/Mapping_file_surgery_data Nature communication overlap.csv"

deseq_sets <- list(
  prev = "raw_data/prev",
  current = "raw_data/current"
)

read_expr <- function(path) {
  expr <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(expr)[1] <- "gene"
  expr
}
drop_missing_rows <- function(x) x[stats::complete.cases(x), , drop = FALSE]
extract_numeric_id <- function(x) suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))

make_model_sets <- function(expr, clinics, mod20, outdir, restrict_samples = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  expr_t <- data.frame(t(expr[, -1, drop = FALSE]), check.names = FALSE)
  names(expr_t) <- expr$gene
  expr_t$Study_ID <- rownames(expr_t)

  common_samples <- intersect(clinics$Study_ID, expr_t$Study_ID)
  if (!is.null(restrict_samples)) {
    common_samples <- intersect(restrict_samples, common_samples)
  }
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

metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
clinic_cols <- c("Study_ID", "recurrence", "TTP", "ProgType", "Tissue_Sex", "Tissue_Age",
                 "Tissue_race", "Tissue_smoke_status", "Tissue_Tumor_size")
clinics_whole <- metadata[!is.na(metadata$recurrence), clinic_cols]
rownames(clinics_whole) <- clinics_whole$Study_ID
clinics_stage1 <- metadata[metadata$stage_lungca %in% c("IA1", "IA2", "IB") & !is.na(metadata$recurrence), clinic_cols]

module_genes <- read.csv(module_gene_file, check.names = FALSE, stringsAsFactors = FALSE)
module20_genes <- module_genes[module_genes[, 2] == 20, 1]
gene_type_map <- read.csv(gene_type_file, check.names = FALSE, stringsAsFactors = FALSE)

# --- reproduce `intersection` exactly as the main builder defines it ---------
old_scores <- read.csv(old_module_score_file, stringsAsFactors = FALSE, check.names = FALSE)
old_ids <- extract_numeric_id(names(old_scores)[grep("-N", names(old_scores))])
current_numeric_ids <- extract_numeric_id(clinics_stage1$Study_ID)
intersection_numeric_ids <- intersect(old_ids, current_numeric_ids)
intersection_samples <- clinics_stage1$Study_ID[match(intersection_numeric_ids, current_numeric_ids)]
intersection_samples <- intersection_samples[!is.na(intersection_samples)]

# the requested cohort: wholeset minus intersection
nonintersection_samples <- setdiff(clinics_whole$Study_ID, intersection_samples)

# --- full NatComm set, to quantify the residual contamination ----------------
whole_ids <- extract_numeric_id(clinics_whole$Study_ID)
in_a <- clinics_whole$Study_ID[whole_ids %in% old_ids]              # 66
nc_map <- read.csv(nc_mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
flag_col <- grep("^Nature Comm", names(nc_map), value = TRUE)[1]
id_col   <- grep("Segal", names(nc_map), value = TRUE)[1]
flag_ids <- extract_numeric_id(nc_map[[id_col]][trimws(nc_map[[flag_col]]) == "Nature"])
in_b <- clinics_whole$Study_ID[whole_ids %in% flag_ids]             # 75
natcomm_samples <- union(in_a, in_b)                               # 76

cat(sprintf("intersection            = %d\n", length(intersection_samples)))
cat(sprintf("nonintersection (150-58)= %d\n", length(nonintersection_samples)))
cat(sprintf("full NatComm (union)    = %d  [(a) old RNA-seq %d ; (b) 'Nature' flag %d]\n",
            length(natcomm_samples), length(in_a), length(in_b)))
cat(sprintf(">>> NatComm patients STILL INSIDE nonintersection: %d\n",
            length(intersect(nonintersection_samples, natcomm_samples))))
stopifnot(length(nonintersection_samples) == 92)

compute_mod20 <- function(deseq_dir) {
  vst <- read_expr(file.path(deseq_dir, "counts.vst.csv"))
  vst_mat <- as.matrix(vst[, -1, drop = FALSE])
  rownames(vst_mat) <- vst$gene
  scaled_expr <- t(apply(vst_mat, 1, scale))
  mod20 <- data.frame(
    module20 = colMeans(scaled_expr[rownames(scaled_expr) %in% module20_genes, , drop = FALSE], na.rm = TRUE))
  rownames(mod20) <- colnames(vst_mat)
  mod20
}
coding_filter <- function(expr) {
  gene_type <- gene_type_map$gene_type[match(expr$gene, gene_type_map$gene_name)]
  expr[!is.na(gene_type) & gene_type == "protein_coding", , drop = FALSE]
}

NOTE <- paste0("nonintersection = wholeset(150) - intersection(58). WARNING: `intersection` is ",
               "(old NatComm normal RNA-seq) AND (stage-1), NOT the full NatComm overlap. Full NatComm ",
               "contact = 76 (union of 66 with old NatComm normal RNA-seq and 75 flagged 'Nature'). ",
               "Therefore 18 Nature Communications patients REMAIN in this cohort - it is NOT NatComm-free. ",
               "For a leakage-free set use cohort `novel` (n=74, 0 NatComm).")

summary_rows <- list()
for (set_name in names(deseq_sets)) {
  deseq_dir <- deseq_sets[[set_name]]
  cat("\n==== RNA-seq set:", set_name, "====\n")
  mod20 <- compute_mod20(deseq_dir)
  expr_files <- list(fpkm = file.path(deseq_dir, "counts.fpkm.csv"),
                     normalized = file.path(deseq_dir, "counts.norm.csv"))
  for (data_type in names(expr_files)) {
    expr <- read_expr(expr_files[[data_type]])
    expr_coding <- coding_filter(expr)
    make_model_sets(expr, clinics_whole, mod20,
                    file.path(out_root, set_name, "nonintersection", data_type),
                    restrict_samples = nonintersection_samples)
    make_model_sets(expr_coding, clinics_whole, mod20,
                    file.path(out_root, set_name, "nonintersection_coding", data_type),
                    restrict_samples = nonintersection_samples)
    for (cohort in c("nonintersection", "nonintersection_coding")) {
      for (task in c("123vs0", "12vs0", "3vs0")) {
        d <- readRDS(file.path(out_root, set_name, cohort, data_type, paste0("fpkm_", task, ".RDS")))
        ids <- rownames(d)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          set = set_name, cohort = cohort, data = task, type = data_type,
          No.Recur = sum(d$recurrence == 1), No.NonRecu = sum(d$recurrence == 0), Total = nrow(d),
          n_NatComm_remaining = length(intersect(ids, natcomm_samples)),
          pct_NatComm_remaining = round(100 * length(intersect(ids, natcomm_samples)) / nrow(d), 1),
          natcomm_free = FALSE,
          note = NOTE,
          stringsAsFactors = FALSE)
      }
    }
    cat(sprintf("  %s / %s: done\n", set_name, data_type))
  }
}
summary_table <- do.call(rbind, summary_rows)
out_csv <- file.path(out_root, "counts_summary_nonintersection.csv")
write.csv(summary_table, out_csv, row.names = FALSE)
cat("\nWROTE SUMMARY CSV: ", normalizePath(out_csv), "\n\n")
print(summary_table[, c("set","cohort","data","type","No.Recur","No.NonRecu","Total",
                        "n_NatComm_remaining","pct_NatComm_remaining")])
