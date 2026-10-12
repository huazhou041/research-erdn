setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
suppressMessages({library(glmnet)})
CLIN <- c("(Intercept)","Tissue_Sex","Tissue_Age","Tissue_race","Tissue_smoke_status","Tissue_Tumor_size","Tissue_Module20")
glm_used <- function(f){
  g <- readRDS(f)
  b <- tryCatch(as.matrix(coef(g, s="lambda.min")), error=function(e) NULL)
  if (is.null(b)) return(list(err="coef failed", class=class(g)[1]))
  nz <- b[b[,1]!=0,,drop=FALSE]
  gn <- setdiff(rownames(nz), CLIN)
  list(offered=length(setdiff(rownames(b),CLIN)), used=length(gn), top=utils::head(gn,6))
}
rfxgb_used <- function(f){
  x <- readRDS(f)
  feat <- setdiff(as.character(x$features), CLIN)
  imp <- tryCatch({
    m <- x$model
    if (inherits(m,"xgb.Booster")) { xgboost::xgb.importance(model=m)$Feature }
    else if (inherits(m,"randomForest")) { rownames(m$importance)[order(-m$importance[,ncol(m$importance)])] }
    else NULL
  }, error=function(e) paste("imp-err:",conditionMessage(e)))
  list(offered=length(feat), imp_head=if(is.character(imp)) utils::head(imp,6) else NA)
}
cat("--- GLM samples ---\n")
for (f in head(list.files("res", pattern="_glm_(mad|cor|natcommad)_v200.rds$", full.names=TRUE), 4)) { cat(basename(f),"\n"); print(glm_used(f)) }
cat("--- RF/XGB samples ---\n")
for (f in c(head(list.files("res_rfxgb", pattern="_rf_v200.rds$", full.names=TRUE),2), head(list.files("res_rfxgb", pattern="_xgb_v200.rds$", full.names=TRUE),2))) { cat(basename(f),"\n"); print(rfxgb_used(f)) }
