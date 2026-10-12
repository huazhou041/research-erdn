# Top-MAD gene overlap (nature_com Normal vs prev vs current) at several
# varcut thresholds x cohort (wholeset/stage1/intersection, each with a
# protein-coding-only or including-non-coding variant -- 6 total, matching
# the dashboard's own coh_lv exactly). nature_com's PATIENT set is fixed
# (flat) across all 6 cohort choices -- it has no wholeset/stage1/
# intersection subdivision of its own -- but its GENE universe still
# respects the coding/non-coding choice, so each cohort's 3-way comparison
# stays apples-to-apples on gene universe.
#
# nature_com patient set: the 123 patients in
# manuscript/figures/data/models/data/normalized/matched/control/
# control_fpkm_123vs0.RDS -- the "matched" (Normal+Tumor paired, NOT
# rm5years-restricted) cohort's 123vs0 subset. Confirmed exact match to a
# number the user recalled ("123 or 126") from the original nature_com
# training set -- traced by checking nrow() across every cohort x subset
# RDS variant; this is the only exact match found.
#
# nature_com's raw expression source is still the corrected 270-sample
# DESeq2 run (manuscript/ml/data/DGE-DESeq2-unstr-270samples-2groups/
# counts.normalized.csv) -- confirmed via manuscript/ml/getdata.R as the
# actual file that trained the AUC models (Section 2/3), NOT the
# revise_202302 246-sample re-export used in the earlier, uncorrected
# version of this script.
rm(list = ls())
suppressMessages(library(tidyverse))

# CONFIG (override via env). GENECORR_RAW must hold the regenerated expression
# inputs; the original ephemeral session scratchpad no longer exists.
proj <- Sys.getenv("ERDN_PROJ", ".")
outdir <- Sys.getenv("GENECORR_OUT", file.path(proj, "htmls"))
gc_raw <- Sys.getenv("GENECORR_RAW", file.path(proj, "provenance_inputs", "genecorr_raw"))
stopifnot(dir.exists(gc_raw))

gene_type_file <- Sys.getenv("GENE_TYPE_FILE", "raw_data/gene_type.csv")
gene_type_map <- read.csv(gene_type_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

read_expr <- function(path) { x <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE); names(x)[1] <- "gene"; x }
coding_only <- function(expr) { gt <- gene_type_map$gene_type[match(expr$gene, gene_type_map$gene_name)]; expr[!is.na(gt) & gt == "protein_coding", , drop = FALSE] }

mad_rank <- function(expr) {
  m <- as.matrix(expr[, -1]); storage.mode(m) <- "double"
  m <- log2(m + 1e-4)
  v <- apply(m, 1, mad, na.rm = TRUE)
  names(v) <- expr$gene
  sort(v, decreasing = TRUE)
}

## ---- nature_com: fixed 123-patient set, corrected 270-sample expression source ----
nature_com <- read_expr(Sys.getenv("NATURECOM_NORM_FILE", "raw_data/naturecom_counts.normalized.csv"))
nature_com_N <- nature_com[, c(1, grep("-N$", names(nature_com)))]
names(nature_com_N) <- gsub("NYU|-N$", "", names(nature_com_N))

nc123_patients <- rownames(readRDS(Sys.getenv("NATURECOM_123_RDS", "raw_data/naturecom_control_fpkm_123vs0.RDS")))
nc_present <- intersect(names(nature_com_N)[-1], nc123_patients)
cat(sprintf("nature_com 'matched/123vs0' target patient list: %d; found in raw 270-sample source: %d\n",
            length(nc123_patients), length(nc_present)))
nature_com_fixed <- nature_com_N[, c("gene", nc_present)]

## ---- prev/current: full 168-sample matrices, then subset per cohort ----
prev_full    <- read_expr("raw_data/prev/counts.norm.csv")
current_full <- read_expr("raw_data/current/counts.norm.csv")
names(prev_full)    <- gsub("\\.NL$", "", names(prev_full))
names(current_full) <- gsub("\\.NL$", "", names(current_full))
cat(sprintf("prev: %d samples | current: %d samples (full, pre-cohort-filter)\n", ncol(prev_full) - 1, ncol(current_full) - 1))

## cohort membership = rownames of that base cohort's own fpkm_123vs0.RDS
## (broadest-inclusion subset; "_coding" variants share the same patients as
## their plain counterpart -- only the gene-type filter differs, per Section 5).
cohort_patients <- function(set, base_cohort) {
  rds <- readRDS(file.path(gc_raw, set, base_cohort, "fpkm", "fpkm_123vs0.RDS"))
  gsub("\\.NL$", "", rownames(rds))
}

coh_lv <- c("wholeset", "wholeset_coding", "stage1", "stage1_coding", "intersection", "intersection_coding")
varcuts <- c(50, 100, 200, 500, 1000)
venn_rows <- list()
overlap_genes <- list()

for (coh in coh_lv) {
  is_coding <- grepl("_coding$", coh)
  base_coh <- sub("_coding$", "", coh)

  prev_pat <- cohort_patients("prev", base_coh)
  cur_pat  <- cohort_patients("current", base_coh)
  prev_sub    <- prev_full[,    c("gene", intersect(names(prev_full)[-1], prev_pat))]
  current_sub <- current_full[, c("gene", intersect(names(current_full)[-1], cur_pat))]
  nc_sub <- nature_com_fixed

  if (is_coding) { prev_sub <- coding_only(prev_sub); current_sub <- coding_only(current_sub); nc_sub <- coding_only(nc_sub) }

  cat(sprintf("\ncohort=%s (coding_only=%s): nature_com n=%d genes=%d | prev n=%d genes=%d | current n=%d genes=%d\n",
              coh, is_coding, ncol(nc_sub) - 1, nrow(nc_sub), ncol(prev_sub) - 1, nrow(prev_sub), ncol(current_sub) - 1, nrow(current_sub)))

  nc_mad   <- mad_rank(nc_sub)
  prev_mad <- mad_rank(prev_sub)
  cur_mad  <- mad_rank(current_sub)

  for (vc in varcuts) {
    A <- names(nc_mad)[seq_len(min(vc, length(nc_mad)))]
    B <- names(prev_mad)[seq_len(min(vc, length(prev_mad)))]
    C <- names(cur_mad)[seq_len(min(vc, length(cur_mad)))]

    abc <- Reduce(intersect, list(A, B, C))
    ab  <- setdiff(intersect(A, B), C)
    ac  <- setdiff(intersect(A, C), B)
    bc  <- setdiff(intersect(B, C), A)
    a_only <- setdiff(A, union(B, C))
    b_only <- setdiff(B, union(A, C))
    c_only <- setdiff(C, union(A, B))

    key <- paste(coh, vc)
    venn_rows[[key]] <- data.frame(
      cohort = coh, varcut = vc,
      nature_com_only = length(a_only), prev_only = length(b_only), current_only = length(c_only),
      nature_com_prev = length(ab), nature_com_current = length(ac), prev_current = length(bc),
      all_three = length(abc), stringsAsFactors = FALSE)
    overlap_genes[[key]] <- data.frame(cohort = coh, varcut = vc, gene = abc, stringsAsFactors = FALSE)
    cat(sprintf("  v%d: nc_only=%d prev_only=%d cur_only=%d nc&prev=%d nc&cur=%d prev&cur=%d all3=%d\n",
                vc, length(a_only), length(b_only), length(c_only), length(ab), length(ac), length(bc), length(abc)))
  }
}

venn_tab <- do.call(rbind, venn_rows)
write.csv(venn_tab, file.path(outdir, "venn_madgenes.csv"), row.names = FALSE)
overlap_tab <- do.call(rbind, overlap_genes)
write.csv(overlap_tab, file.path(outdir, "venn_madgenes_all3_genelist.csv"), row.names = FALSE)
cat("\nwrote venn_madgenes.csv and venn_madgenes_all3_genelist.csv (", nrow(venn_tab), "rows )\n")
print(venn_tab)
