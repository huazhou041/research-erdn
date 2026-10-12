setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
# 1) the actual sample/event count for current/novel/3vs0 (the modeled data)
f <- "data/model_data/current/novel/normalized/fpkm_3vs0.RDS"
d <- readRDS(f)
cat("=== current/novel/normalized fpkm_3vs0.RDS ===\n")
cat("n samples:", nrow(d), " n genes(features):", ncol(d)-1, "\n")
cat("recurrence table:\n"); print(table(d$recurrence))
cat("events (recurrence==1):", sum(d$recurrence==1), " non-events:", sum(d$recurrence==0), "\n")

# 2) natcommad AUCs for this cell across ALL varcuts + its sd and n_auc
r <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors=FALSE)
sub <- r[r$set=="current" & r$cohort=="novel" & r$datatype=="normalized" &
         r$featureset=="fpkm" & r$subset=="3vs0", c("fsmethod","varcut","mean_auc","sd_auc","n_auc")]
sub <- sub[order(sub$fsmethod, sub$varcut),]
cat("\n=== all fsmethods, current/novel/normalized/fpkm/3vs0, by varcut ===\n")
print(sub, row.names=FALSE)

# 3) for contrast: same natcommad model in a WELL-POWERED cell (wholeset 123vs0)
cat("\n=== natcommad in a well-powered cell (current/wholeset/normalized/fpkm/123vs0) ===\n")
w <- r[r$set=="current" & r$cohort=="wholeset" & r$datatype=="normalized" &
       r$featureset=="fpkm" & r$subset=="123vs0" & r$fsmethod=="natcommad", c("varcut","mean_auc","sd_auc","n_auc")]
print(w[order(w$varcut),], row.names=FALSE)
