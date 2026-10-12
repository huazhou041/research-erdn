rm(list = ls())

# ---------------------------------------------------------------------------
# Build the 3 model-data cohorts (wholeset / stage1 / intersection) for TWO
# RNA-seq processings of the SAME 168 samples:
#   prev    = erdn_202502/deseq2  (previous RNA-seq run)
#   current = erdn_202606/deseq2  (current  RNA-seq run)
# Output: erdn_202607/data/model_data/{prev,current}/{wholeset,stage1,intersection}/{fpkm,normalized}/*.RDS
# Logic mirrors erdn_202606/codes/build_modeldata_*.R (module20 + all feature varieties).
# ---------------------------------------------------------------------------

base_dir <- normalizePath(".", mustWork = TRUE)
out_root <- file.path(base_dir, "data", "model_data")

# Shared inputs (metadata + references are identical for both RNA-seq runs) ----
metadata_file <- "raw_data/metadata_model_202606.csv"
gene_type_file <- "raw_data/gene_type.csv"
module_gene_file <- "raw_data/246samples.coding10000.modules.dim10.k20.genes.csv"
old_normalized_file <- "raw_data/old_counts.normalized.csv"
old_module_score_file <- "raw_data/246samples.coding10000.modules.dim10.k20.scores.csv"

# The two RNA-seq processings to run through the same pipeline -----------------
deseq_sets <- list(
  prev = "raw_data/prev",
  current = "raw_data/current"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
read_expr <- function(path) {
  expr <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(expr)[1] <- "gene"
  expr
}

drop_missing_rows <- function(x) {
  x[stats::complete.cases(x), , drop = FALSE]
}

extract_numeric_id <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9]", "", x)))
}

safe_cor_row <- function(gene, old_values, new_values) {
  keep <- is.finite(old_values) & is.finite(new_values)
  if (sum(keep) < 3) {
    return(data.frame(gene = gene, correlation = NA_real_, pval = NA_real_))
  }
  ct <- suppressWarnings(cor.test(old_values[keep], new_values[keep]))
  data.frame(gene = gene, correlation = unname(ct$estimate), pval = ct$p.value)
}

# Write the full set of feature varieties x subsets for one expression table.
# `restrict_samples` (optional) limits to the intersection sample set.
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

  task_defs <- list(
    "12vs0" = c(0, 1, 2),
    "3vs0" = c(0, 3)
  )
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

# ---------------------------------------------------------------------------
# Shared metadata + clinics views
# ---------------------------------------------------------------------------
metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
clinic_cols <- c("Study_ID", "recurrence", "TTP", "ProgType", "Tissue_Sex", "Tissue_Age",
                 "Tissue_race", "Tissue_smoke_status", "Tissue_Tumor_size")

clinics_whole <- metadata[!is.na(metadata$recurrence), clinic_cols]
rownames(clinics_whole) <- clinics_whole$Study_ID

clinics_stage1 <- metadata[metadata$stage_lungca %in% c("IA1", "IA2", "IB") & !is.na(metadata$recurrence), clinic_cols]
rownames(clinics_stage1) <- clinics_stage1$Study_ID

module_genes <- read.csv(module_gene_file, check.names = FALSE, stringsAsFactors = FALSE)
module20_genes <- module_genes[module_genes[, 2] == 20, 1]
gene_type_map <- read.csv(gene_type_file, check.names = FALSE, stringsAsFactors = FALSE)

# Old-cohort module-20 scores + numeric IDs used to define the intersection ---
old_scores <- read.csv(old_module_score_file, stringsAsFactors = FALSE, check.names = FALSE)
old_normal_scores <- old_scores[, grep("-N", names(old_scores)), drop = FALSE]
old_mod20 <- data.frame(module20 = as.numeric(old_normal_scores[20, ]))
rownames(old_mod20) <- names(old_normal_scores)
old_numeric_ids <- extract_numeric_id(rownames(old_mod20))

old_expr <- read.csv(old_normalized_file, check.names = FALSE, stringsAsFactors = FALSE)
old_expr <- old_expr[, c(1, grep("-N", names(old_expr))), drop = FALSE]
names(old_expr) <- gsub("[-|NYU|CH|S]", "", names(old_expr))

current_numeric_ids <- extract_numeric_id(clinics_stage1$Study_ID)
intersection_numeric_ids <- intersect(old_numeric_ids, current_numeric_ids)
intersection_samples <- clinics_stage1$Study_ID[match(intersection_numeric_ids, current_numeric_ids)]
intersection_samples <- intersection_samples[!is.na(intersection_samples)]

# ---------------------------------------------------------------------------
# Per RNA-seq-set module20 score (from that set's own VST)
# ---------------------------------------------------------------------------
compute_mod20 <- function(deseq_dir) {
  vst <- read_expr(file.path(deseq_dir, "counts.vst.csv"))
  vst_mat <- as.matrix(vst[, -1, drop = FALSE])
  rownames(vst_mat) <- vst$gene
  scaled_expr <- t(apply(vst_mat, 1, scale))
  mod20 <- data.frame(
    module20 = colMeans(scaled_expr[rownames(scaled_expr) %in% module20_genes, , drop = FALSE], na.rm = TRUE)
  )
  rownames(mod20) <- colnames(vst_mat)
  mod20
}

coding_filter <- function(expr) {
  gene_type <- gene_type_map$gene_type[match(expr$gene, gene_type_map$gene_name)]
  expr[!is.na(gene_type) & gene_type == "protein_coding", , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Drive both RNA-seq sets through the three cohorts
# ---------------------------------------------------------------------------
summary_rows <- list()
cor_rows_all <- list()

for (set_name in names(deseq_sets)) {
  deseq_dir <- deseq_sets[[set_name]]
  cat("\n==== RNA-seq set:", set_name, "(", deseq_dir, ") ====\n")
  mod20 <- compute_mod20(deseq_dir)

  expr_files <- list(
    fpkm = file.path(deseq_dir, "counts.fpkm.csv"),
    normalized = file.path(deseq_dir, "counts.norm.csv")
  )

  for (data_type in names(expr_files)) {
    expr <- read_expr(expr_files[[data_type]])
    expr_coding <- coding_filter(expr)

    # all-genes vs protein-coding done separately for wholeset and stage1;
    # intersection is protein-coding (+ old-cohort intersection) by definition.
    # 1) wholeset        : no stage filter, all genes
    make_model_sets(expr, clinics_whole, mod20,
                    file.path(out_root, set_name, "wholeset", data_type))
    # 2) wholeset_coding : no stage filter, protein-coding
    make_model_sets(expr_coding, clinics_whole, mod20,
                    file.path(out_root, set_name, "wholeset_coding", data_type))
    # 3) stage1          : stage IA1/IA2/IB, all genes
    make_model_sets(expr, clinics_stage1, mod20,
                    file.path(out_root, set_name, "stage1", data_type))
    # 4) stage1_coding   : stage IA1/IA2/IB, protein-coding
    make_model_sets(expr_coding, clinics_stage1, mod20,
                    file.path(out_root, set_name, "stage1_coding", data_type))
    # 5) intersection        : stage I, all genes, intersected with old 246 cohort
    make_model_sets(expr, clinics_stage1, mod20,
                    file.path(out_root, set_name, "intersection", data_type),
                    restrict_samples = intersection_samples)
    # 6) intersection_coding : stage I, protein-coding, intersected with old 246 cohort
    make_model_sets(expr_coding, clinics_stage1, mod20,
                    file.path(out_root, set_name, "intersection_coding", data_type),
                    restrict_samples = intersection_samples)

    # Counts summary per cohort/subset
    for (cohort in c("wholeset", "wholeset_coding", "stage1", "stage1_coding",
                     "intersection", "intersection_coding")) {
      for (task in c("123vs0", "12vs0", "3vs0")) {
        d <- readRDS(file.path(out_root, set_name, cohort, data_type, paste0("fpkm_", task, ".RDS")))
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          set = set_name, cohort = cohort, data = task, type = data_type,
          No.Recur = sum(d$recurrence == 1), No.NonRecu = sum(d$recurrence == 0), Total = nrow(d)
        )
      }
    }

    # Optional QC: correlation of top-MAD genes vs old manuscript data (intersection, normalized only)
    if (data_type == "normalized") {
      expr_t <- data.frame(t(expr_coding[, -1, drop = FALSE]), check.names = FALSE)
      names(expr_t) <- expr_coding$gene
      rownames(expr_t) <- names(expr_coding)[-1]
      common_samples <- intersect(intersection_samples, rownames(expr_t))
      topdata <- expr_t[match(common_samples, rownames(expr_t)), , drop = FALSE]
      ix <- apply(topdata, 2, mad, na.rm = TRUE)
      topgenes <- names(topdata)[order(ix, decreasing = TRUE)[seq_len(min(50, length(ix)))]]
      current_numeric <- extract_numeric_id(common_samples)
      old_common <- intersect(intersection_numeric_ids, current_numeric)

      cor_tab <- do.call(rbind, lapply(topgenes, function(gene) {
        old_values <- as.numeric(old_expr[match(gene, old_expr[, 1]), match(old_common, names(old_expr))])
        new_values <- as.numeric(topdata[match(old_common, current_numeric), match(gene, names(topdata))])
        safe_cor_row(gene, old_values, new_values)
      }))
      old_module_values <- old_mod20[match(old_common, old_numeric_ids), 1]
      new_module_values <- mod20[match(common_samples[match(old_common, current_numeric)], rownames(mod20)), 1]
      cor_tab <- rbind(cor_tab, safe_cor_row("Module20", old_module_values, new_module_values))
      cor_tab$set <- set_name
      cor_rows_all[[length(cor_rows_all) + 1]] <- cor_tab
    }
  }
}

summary_table <- do.call(rbind, summary_rows)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table, file.path(out_root, "counts_summary.csv"), row.names = FALSE)

if (length(cor_rows_all) > 0) {
  cor_table <- do.call(rbind, cor_rows_all)
  write.csv(cor_table, file.path(out_root, "correlation_summary_intersection.csv"), row.names = FALSE)
}

cat("\nWrote model data under:", out_root, "\n")
cat("Intersection samples:", length(intersection_samples), "\n\n")
print(summary_table)
