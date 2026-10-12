import csv, io, sys
f="results_roc_curves_fsmethods.csv"
rows=list(csv.DictReader(io.open(f,encoding="utf-8")))
print("columns:", list(rows[0].keys()))
def area(fpr,tpr):
    a=0.0
    for i in range(1,len(fpr)): a+=(fpr[i]-fpr[i-1])*(tpr[i]+tpr[i-1])/2.0
    return a
# AUC column: 'mean' holds the per-fold-averaged AUC in these files
examples=[]
for r in rows:
    try:
        fp=[float(x) for x in r['fpr'].split('|')]; tp=[float(x) for x in r['tpr'].split('|')]
    except Exception: continue
    if len(fp)<3: continue
    ar=abs(area(fp,tp))  # raw curve area (orientation as drawn)
    auc=float(r.get('mean') or r.get('auc') or 'nan')
    if auc>=0.8 and ar<0.5:
        examples.append((auc,ar,r['set'],r['cohort'],r['subset'],r['featureset'],r['datatype'],r['model'],r['varcut'],r.get('n_pos_pred'),r.get('n_neg_pred')))
examples.sort(key=lambda x:-x[0])
print("configs with per-fold AUC>=0.8 BUT raw pooled curve area<0.5:", len(examples))
for e in examples[:8]:
    print(f"  AUC(perfold)={e[0]:.3f}  pooledArea={e[1]:.3f}  {e[2]}/{e[3]}/{e[4]}/{e[5]}/{e[6]} {e[7]} v{e[8]}  pos={e[9]} neg={e[10]}")
