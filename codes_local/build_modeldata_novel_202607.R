rm(list = ls())

# ---------------------------------------------------------------------------
# Add the `novel` cohort: the current/prev samples that have NO contact with the
# Nature Communications cohort. This is the independent set needed to test whether
# NatComm-derived features (notably Tissue_Module20, whose 170 genes were defined on
# the old 246-sample cohort) hold up on patients that did not contribute to them.
#
# Membership is the wholeset MINUS the union of two definitions of "in NatComm":
#   (a) has a normal-tissue ("-N") column in the old NatComm expression matrix  -> 66
#   (b) flagged "Nature" in Mapping_file_surgery_data Nature communication overlap -> 75
# The two are NOT nested (1 patient in (a) only, 10 in (b) only), so the union (76)
# is the conservative exclusion. novel = 150 - 76 = 74.
#
# NOTE this is deliberately NOT "wholeset minus intersection". `intersection` is only 58
# because it is additionally restricted to stage-1; subtracting it would leave 92 samples
# that still include 18 NatComm patients, defeating the purpose.
#
# Everything else (make_model_sets, module20 score, coding filter, task splits) is
# byte-identical to build_modeldata_prev_current_202607.R so the new cohort is directly
# comparable to wholeset/stage1/intersection.
# Output: data/model_data/{prev,current}/{novel,novel_coding}/{fpkm,normalized}/*.RDS
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

# --- helpers (identical to the main builder) --------------------------------
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

# --- metadata ---------------------------------------------------------------
metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
clinic_cols <- c("Study_ID", "recurrence", "TTP", "ProgType", "Tissue_Sex", "Tissue_Age",
                 "Tissue_race", "Tissue_smoke_status", "Tissue_Tumor_size")
clinics_whole <- metadata[!is.na(metadata$recurrence), clinic_cols]
rownames(clinics_whole) <- clinics_whole$Study_ID

module_genes <- read.csv(module_gene_file, check.names = FALSE, stringsAsFactors = FALSE)
module20_genes <- module_genes[module_genes[, 2] == 20, 1]
gene_type_map <- read.csv(gene_type_file, check.names = FALSE, stringsAsFactors = FALSE)

# --- define the NatComm exclusion set ---------------------------------------
whole_ids <- extract_numeric_id(clinics_whole$Study_ID)

# (a) patients carrying a normal-tissue column in the old NatComm expression matrix
old_scores <- read.csv(old_module_score_file, stringsAsFactors = FALSE, check.names = FALSE)
old_ids <- extract_numeric_id(names(old_scores)[grep("-N", names(old_scores))])
in_a <- clinics_whole$Study_ID[whole_ids %in% old_ids]

# (b) patients flagged "Nature" in the mapping file. read.csv mangles the header
# ("Nature Comm overlap " has a trailing space), so locate the columns by pattern.
nc_map <- read.csv(nc_mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
flag_col <- grep("^Nature Comm", names(nc_map), value = TRUE)[1]
id_col   <- grep("Segal", names(nc_map), value = TRUE)[1]
stopifnot(!is.na(flag_col), !is.na(id_col))
flag_ids <- extract_numeric_id(nc_map[[id_col]][trimws(nc_map[[flag_col]]) == "Nature"])
in_b <- clinics_whole$Study_ID[whole_ids %in% flag_ids]

natcomm_samples <- union(in_a, in_b)
novel_samples <- setdiff(clinics_whole$Study_ID, natcomm_samples)
cat(sprintf("NatComm contact: (a) old normal RNA-seq = %d ; (b) 'Nature' flag = %d ; union = %d\n",
            length(in_a), length(in_b), length(natcomm_samples)))
cat(sprintf("NOVEL cohort = %d samples (wholeset %d - NatComm %d)\n",
            length(novel_samples), nrow(clinics_whole), length(natcomm_samples)))
stopifnot(length(novel_samples) == 74)

compute_mod20 <- function(deseq_dir) {
  vst <- read_expr(file.path(deseq_dir, "counts.vst.csv"))
  vst_mat <- as.matrix(vst[, -1, drop = FALSE])
  rownames(vst_mat) <- vst$gene
  scaled_expr <- t(apply(vst_mat, 1, scale))
  hit <- rownames(scaled_expr) %in% module20_genes
  cat(sprintf("    module20 genes found in this VST: %d / %d\n", sum(hit), length(module20_genes)))
  mod20 <- data.frame(module20 = colMeans(scaled_expr[hit, , drop = FALSE], na.rm = TRUE))
  rownames(mod20) <- colnames(vst_mat)
  mod20
}
coding_filter <- function(expr) {
  gene_type <- gene_type_map$gene_type[match(expr$gene, gene_type_map$gene_name)]
  expr[!is.na(gene_type) & gene_type == "protein_coding", , drop = FALSE]
}

# --- build ------------------------------------------------------------------
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
                    file.path(out_root, set_name, "novel", data_type),
                    restrict_samples = novel_samples)
    make_model_sets(expr_coding, clinics_whole, mod20,
                    file.path(out_root, set_name, "novel_coding", data_type),
                    restrict_samples = novel_samples)
    for (cohort in c("novel", "novel_coding")) {
      for (task in c("123vs0", "12vs0", "3vs0")) {
        d <- readRDS(file.path(out_root, set_name, cohort, data_type, paste0("fpkm_", task, ".RDS")))
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          set = set_name, cohort = cohort, data = task, type = data_type,
          No.Recur = sum(d$recurrence == 1), No.NonRecu = sum(d$recurrence == 0), Total = nrow(d))
      }
    }
    cat(sprintf("  %s / %s: done\n", set_name, data_type))
  }
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(out_root, "counts_summary_novel.csv"), row.names = FALSE)
cat("\n")
print(summary_table)
