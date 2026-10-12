import csv, io
f="roc_foldavg_202609.csv"
rows=list(csv.DictReader(io.open(f,encoding="utf-8")))
print("rows:", len(rows))
below=[r for r in rows if float(r['auc_vavg'])<0.5]
print("vertical-avg AUC < 0.5 (should be rare/none for real signal):", len(below))
# compare a few to headline AUC
import collections
# load headline AUCs (GLM)
glm={}
for pop,fn in [("full","results_fsmethod_comparison.csv"),("smoke","results_fsmethod_comparison_smoking.csv")]:
    for r in csv.DictReader(io.open(fn,encoding="utf-8")):
        if r['featureset'] in ('fpkm','fpkm_clinics'):
            k=pop+"|"+"_".join([r['set'],r['cohort'],r['datatype'],r['featureset'],r['subset'],"glm_"+r['fsmethod'],"v"+r['varcut']])
            glm[k]=float(r['mean_auc'])
famap={ r['pop']+"|"+r['key']: float(r['auc_vavg']) for r in rows }
import random
random.seed(1)
shared=[k for k in famap if k in glm]
print("shared GLM configs:", len(shared))
diffs=[abs(famap[k]-glm[k]) for k in shared]
print("mean|vavg-headline|=%.4f  max=%.4f"%(sum(diffs)/len(diffs), max(diffs)))
# the previously-flagged novel/3vs0 natcommad
for k in famap:
    if "novel_normalized_fpkm_3vs0_glm_natcommad_v1000" in k:
        print("novel/3vs0 natcommad v1000: vavg AUC=", famap[k], " headline=", glm.get(k))
