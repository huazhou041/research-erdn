import json, csv, io, random, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
H="htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking.html"
s=io.open(H,encoding="utf-8").read()
COLS=["set","cohort","datatype","featureset","subset","model","varcut","mean","sd","n","lo","hi"]
def extract_arrays(txt, key="LONG_DATA:"):
    out=[]
    i=0
    while True:
        j=txt.find(key,i)
        if j<0: break
        k=txt.find("[",j); depth=0; e=k
        while e<len(txt):
            if txt[e]=="[":depth+=1
            elif txt[e]=="]":
                depth-=1
                if depth==0: break
            e+=1
        out.append(json.loads(txt[k:e+1])); i=e+1
    return out
arrs=extract_arrays(s)
print("LONG_DATA blobs found:", len(arrs), "(expect 2: full, smoke)")
full=arrs[0]; smoke=arrs[1]
def to_map(arr):
    m={}
    for row in arr:
        d=dict(zip(COLS,row)); key=(d["set"],d["cohort"],d["datatype"],d["featureset"],d["subset"],d["model"],str(d["varcut"]))
        m[key]=d["mean"]
    return m
fm=to_map(full)
print("full rows:", len(full), " distinct models:", sorted(set(r[COLS.index('model')] for r in full)))
print("natcommad rows in full LONG_DATA:", sum(1 for r in full if r[COLS.index('model')]=='glm_natcommad'))
print("distinct cohorts in full:", len(set(r[COLS.index('cohort')] for r in full)))

def load_csv(fn): 
    return list(csv.DictReader(io.open(fn,encoding="utf-8")))
def check(agg_fn, model_from, meancol, label, n=25):
    rows=load_csv(agg_fn); random.seed(3); samp=random.sample(rows,min(n,len(rows)))
    maxd=0; miss=0; bad=0
    for r in samp:
        model=model_from(r)
        key=(r["set"],r["cohort"],r["datatype"],r["featureset"],r["subset"],model,str(r["varcut"]))
        if key not in fm: miss+=1; continue
        d=abs(float(fm[key])-float(r[meancol])); maxd=max(maxd,d)
        if d>2e-3: bad+=1; print("   MISMATCH",key,"html",fm[key],"agg",r[meancol])
    print(f"[{'PASS' if bad==0 and miss==0 else '**FAIL**'}] {label}: {len(samp)} sampled | max|html-agg|={maxd:.2e} | not-found-in-html={miss} | mismatch={bad}")
check("results_fsmethod_comparison.csv", lambda r:"glm_"+r["fsmethod"], "mean_auc", "GLM full -> HTML full")
check("rf_xgb_auc_202607.csv", lambda r:r["model"]+"_mad", "mean", "RFXGB full -> HTML full")
# smoke side
fm_s=to_map(smoke)
def check_s(agg_fn, model_from, meancol, label, n=20):
    rows=load_csv(agg_fn); random.seed(5); samp=random.sample(rows,min(n,len(rows)))
    maxd=0; miss=0; bad=0
    for r in samp:
        key=(r["set"],r["cohort"],r["datatype"],r["featureset"],r["subset"],model_from(r),str(r["varcut"]))
        if key not in fm_s: miss+=1; continue
        d=abs(float(fm_s[key])-float(r[meancol])); maxd=max(maxd,d)
        if d>2e-3: bad+=1
    print(f"[{'PASS' if bad==0 and miss==0 else '**FAIL**'}] {label}: {len(samp)} sampled | max={maxd:.2e} | not-found={miss} | mismatch={bad}")
check_s("results_fsmethod_comparison_smoking.csv", lambda r:"glm_"+r["fsmethod"], "mean_auc", "GLM smoke -> HTML smoke")
check_s("rf_xgb_auc_smoking.csv", lambda r:r["model"]+"_mad", "mean", "RFXGB smoke -> HTML smoke")
print("AUDIT_HTML_DONE")
