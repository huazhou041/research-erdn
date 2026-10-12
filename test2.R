setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(glmnet)})
g <- readRDS("res/current_wholeset_normalized_fpkm_123vs0_glm_mad_v200.rds")
cat("GLM class:", class(g)[1], "\n")
b <- as.matrix(coef(g, s="lambda.min")); cat("coef nrow:", nrow(b), " rowname sample:", paste(head(rownames(b),6),collapse=","), "\n")
cat("nonzero genes:", length(setdiff(rownames(b)[b[,1]!=0], "(Intercept)")), " of offered", nrow(b)-1, "\n\n")
x <- readRDS("res_rfxgb/current_wholeset_normalized_fpkm_123vs0_xgb_v200.rds")
cat("XGB obj names:", paste(names(x),collapse=","), "\n")
cat("features class:", class(x$features), " len:", length(x$features), " head:", paste(head(x$features,6),collapse=","), "\n")
cat("model class:", class(x$model)[1], "\n")
imp <- tryCatch(xgboost::xgb.importance(model=x$model), error=function(e) conditionMessage(e))
if (is.data.frame(imp)) { cat("xgb.importance Feature head:", paste(head(imp$Feature,6),collapse=","), "\n"); 
  # map f-index to features if needed
  if (all(grepl("^f[0-9]+$", head(imp$Feature)))) { idx <- as.integer(sub("f","",imp$Feature))+1; cat("mapped genes:", paste(head(x$features[idx],6),collapse=","), "\n") } }
r <- readRDS("res_rfxgb/current_wholeset_normalized_fpkm_123vs0_rf_v200.rds")
cat("\nRF features head:", paste(head(r$features,6),collapse=","), " len:", length(r$features), "\n")
cat("RF model class:", class(r$model)[1], "\n")
