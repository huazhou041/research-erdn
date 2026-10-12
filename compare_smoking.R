# Compare regenerated data/model_data_smoking (new Sept builder) vs preserved
# data/model_data_smoking.aug10bak (old row-subset). Confirms lineage change did
# NOT alter the data (so existing tune results remain valid).
setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
NEW <- "data/model_data_smoking"; OLD <- "/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_20260709/data/model_data_smoking"
cohorts <- c("wholeset","wholeset_coding","stage1","stage1_coding",
             "intersection","intersection_coding","nonintersection","nonintersection_coding")
files <- c("fpkm_123vs0.RDS","fpkm_clinics_123vs0.RDS","clinics_123vs0.RDS",
           "module20_123vs0.RDS","fpkm_12vs0.RDS","fpkm_3vs0.RDS")
align_eq <- function(a, b) {
  if (!setequal(rownames(a), rownames(b))) return(list(ok=FALSE, why="rownames differ"))
  if (!setequal(colnames(a), colnames(b))) return(list(ok=FALSE, why="colnames differ"))
  a2 <- a[sort(rownames(a)), sort(colnames(a)), drop=FALSE]
  b2 <- b[sort(rownames(b)), sort(colnames(b)), drop=FALSE]
  eq <- isTRUE(all.equal(a2, b2, check.attributes=FALSE))
  list(ok=eq, why=if(eq) "identical" else "VALUES DIFFER")
}
nfail <- 0; nok <- 0
for (s in c("prev","current")) for (co in cohorts) for (f in files) {
  pn <- file.path(NEW,s,co,"fpkm",f); po <- file.path(OLD,s,co,"fpkm",f)
  if (!file.exists(pn) || !file.exists(po)) { cat(sprintf("MISSING %s/%s/%s new=%s old=%s\n",s,co,f,file.exists(pn),file.exists(po))); nfail<-nfail+1; next }
  r <- align_eq(readRDS(pn), readRDS(po))
  if (r$ok) nok<-nok+1 else { nfail<-nfail+1; cat(sprintf("DIFF %s/%s/%s : %s\n",s,co,f,r$why)) }
}
cat(sprintf("\nCOMPARE: %d identical, %d differing/missing (of %d checks)\n", nok, nfail, nok+nfail))
if (nfail==0) cat("ALL_IDENTICAL: regeneration reproduced the data exactly (lineage fixed, data unchanged)\n")
