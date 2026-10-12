setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
d <- read.csv("results_fsmethod_comparison.csv", stringsAsFactors=FALSE)
# natcommad only has fpkm & fpkm_clinics featuresets; compare like-for-like
d <- d[d$featureset %in% c("fpkm","fpkm_clinics"), ]
cat("=== overall mean AUC by fsmethod (fpkm/fpkm_clinics, all cohorts/tasks/varcuts) ===\n")
print(round(tapply(d$mean_auc, d$fsmethod, mean, na.rm=TRUE), 4))
cat("\n=== best (max mean_auc) per fsmethod ===\n")
print(round(tapply(d$mean_auc, d$fsmethod, max, na.rm=TRUE), 4))
# head-to-head: natcommad vs mad on identical keys
K <- c("set","cohort","datatype","featureset","subset","varcut")
nc <- d[d$fsmethod=="natcommad",]; md <- d[d$fsmethod=="mad",]
kc <- do.call(paste,c(nc[K],sep="|")); km <- do.call(paste,c(md[K],sep="|"))
common <- intersect(kc,km)
nci <- nc[match(common,kc),]; mdi <- md[match(common,km),]
delta <- nci$mean_auc - mdi$mean_auc
cat(sprintf("\n=== natcommad vs mad, matched %d configs ===\n", length(common)))
cat(sprintf("mean delta (natcommad-mad)=%.4f | natcommad better in %.1f%% | mean natcommad=%.4f mean mad=%.4f\n",
    mean(delta,na.rm=TRUE), 100*mean(delta>0,na.rm=TRUE), mean(nci$mean_auc,na.rm=TRUE), mean(mdi$mean_auc,na.rm=TRUE)))
# well-powered headline: current set, 123vs0
cat("\n=== current set, 123vs0, top mean AUC per fsmethod ===\n")
d2 <- d[d$set=="current" & d$subset=="123vs0",]
print(round(tapply(d2$mean_auc, d2$fsmethod, max, na.rm=TRUE),4))
