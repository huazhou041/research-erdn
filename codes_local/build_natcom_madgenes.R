rm(list = ls())
# nature_com top-MAD gene ranking, computed IDENTICALLY to the dashboard Venn
# (make_venn_data.R mad_rank): per-gene MAD of log2(normalized+1e-4) over the
# nature_com NORMAL-tissue samples, restricted to the FIXED 123 nature_com training
# patients, ranked descending. This is the exact "nature_com top-MAD gene list" the
# report cites; it is the FIXED external feature set consumed by tune_drug.R
# fs_method=natcommad. Output: raw_data/natcom_mad_ranked.csv (gene,mad).
base_dir <- normalizePath(".", mustWork = TRUE)
RAW <- Sys.getenv("ERDN_RAW", file.path(base_dir, "raw_data"))
norm_file <- Sys.getenv("NATURECOM_NORM_FILE", file.path(RAW, "naturecom_counts.normalized.csv"))
nc123_rds <- Sys.getenv("NATURECOM_123_RDS",   file.path(RAW, "naturecom_control_fpkm_123vs0.RDS"))
out <- file.path(RAW, "natcom_mad_ranked.csv")
stopifnot(file.exists(norm_file), file.exists(nc123_rds))

read_expr <- function(path) { x <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE); names(x)[1] <- "gene"; x }

# nature_com Normal-tissue expression, restricted to the fixed 123-patient set
nc <- read_expr(norm_file)
nc_N <- nc[, c(1, grep("-N$", names(nc)))]
names(nc_N) <- gsub("NYU|-N$", "", names(nc_N))
nc123 <- rownames(readRDS(nc123_rds))
present <- intersect(names(nc_N)[-1], nc123)
cat(sprintf("nature_com fixed patient set: target %d, found %d in the 270-sample source\n", length(nc123), length(present)))
stopifnot(length(present) > 0)
fixed <- nc_N[, c("gene", present)]

# mad_rank IDENTICAL to make_venn_data.R: log2(x+1e-4) then per-gene MAD, sort desc
m <- as.matrix(fixed[, -1]); storage.mode(m) <- "double"
m <- log2(m + 1e-4)
v <- apply(m, 1, mad, na.rm = TRUE)
res <- data.frame(gene = fixed$gene, mad = v, stringsAsFactors = FALSE)
res <- res[is.finite(res$mad), ]
res <- res[order(res$mad, decreasing = TRUE), ]
res <- res[!duplicated(res$gene), ]        # unique gene symbols, keep highest-MAD
write.csv(res, out, row.names = FALSE)
cat(sprintf("wrote %s : %d ranked genes (top: %s)\n", out, nrow(res), paste(head(res$gene, 10), collapse = ", ")))
