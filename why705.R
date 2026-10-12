setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
d  <- readRDS("data/model_data/current/stage1_coding/normalized/fpkm_12vs0.RDS")
cl <- readRDS("data/model_data/current/stage1_coding/normalized/fpkm_clinics_12vs0.RDS")
cat("n=",nrow(d)," events(recur=1)=",sum(d$recurrence==1)," non=",sum(d$recurrence==0)," genes=",ncol(d)-1,"\n")
# top-200 natcommad genes present in this coding cohort + how many are sex-linked
nr <- read.csv("raw_data/natcom_mad_ranked.csv", stringsAsFactors=FALSE)
ranked <- as.character(nr[[ if("gene"%in%names(nr))"gene" else names(nr)[1] ]])
present <- ranked[ranked %in% colnames(d)]
sel200  <- head(present, 200)
ygenes <- grep("Y[0-9]?$|^RPS4Y|^DDX3Y|^UTY$|^USP9Y|^KDM5D|^TXLNGY|^EIF1AY|^NLGN4Y|^ZFY|^PRKY|^UTX|XIST", sel200, value=TRUE)
cat("top-200 present used:", length(sel200), " | first 12:", paste(head(sel200,12),collapse=", "),"\n")
cat("obvious sex-linked among top-200:", length(ygenes), "->", paste(head(ygenes,15),collapse=", "),"\n")
# is SEX predictive of recurrence in THIS cohort? (the confound natcommad could ride)
if ("Tissue_Sex" %in% names(cl)) { cat("\nsex x recurrence:\n"); print(table(sex=cl$Tissue_Sex, recur=cl$recurrence)) }
# cross-method AUC for this exact cell
r <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors=FALSE)
sub <- r[r$set=="current"&r$cohort=="stage1_coding"&r$datatype=="normalized"&r$featureset=="fpkm"&r$subset=="12vs0", c("fsmethod","varcut","mean_auc","sd_auc","n_auc")]
cat("\n=== all methods, current/stage1_coding/normalized/fpkm/12vs0 ===\n"); print(sub[order(sub$fsmethod,sub$varcut),], row.names=FALSE)
