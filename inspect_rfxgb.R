setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
# an rf and an xgb rds
for (pat in c("res_rfxgb/.*_rf_v200.rds$","res_rfxgb/.*_xgb_v200.rds$")) {
  f <- list.files("res_rfxgb", pattern=sub("res_rfxgb/","",pat), full.names=TRUE)[1]
  if (is.na(f)) { cat("no match for", pat, "\n"); next }
  x <- readRDS(f); cat("\n===", basename(f), "=== class:", class(x)[1], "\n")
  if (is.list(x)) print(utils::head(names(x), 20))
  # try to find importance / feature names
  cn <- tryCatch(dimnames(x)[[1]], error=function(e) NULL)
}
# GLM nonzero example from the natcommad model
g <- readRDS("res/current_wholeset_normalized_fpkm_123vs0_glm_natcommad_v200.rds")
b <- as.matrix(coef(g, s="lambda.min")); nz <- rownames(b)[b[,1]!=0]
cat("\nGLM natcommad wholeset/123vs0/v200: input features kept(nonzero incl intercept):", length(nz), "\n")
cat("first nonzero:", paste(head(setdiff(nz,"(Intercept)"),8), collapse=", "), "\n")
