import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T = "codes_20260713/dashboard_template_fsmethods_cmp.html"
B = "codes_20260713/make_dashboard_fsmethods_cmp.R"

# ---------- template ----------
s = io.open(T, encoding="utf-8").read(); orig = s

# A) GENELISTS const right after LONG_COLS
anchorA = 'const LONG_COLS = ["set","cohort","datatype","featureset","subset","model","varcut","mean","sd","n","lo","hi"];'
if "const GENELISTS" in s:
    print("A const: already present")
elif anchorA in s:
    s = s.replace(anchorA, anchorA + "\nconst GENELISTS = __GENELISTS__;", 1); print("A const: added")
else:
    print("A const: ANCHOR NOT FOUND")

# C1) explorer cols: add Genes column
anchorC1 = "{key:'ci',label:'AUC CI'}];"
newC1 = "{key:'ci',label:'AUC CI'},{key:'__genes',label:'Genes'}];"
if newC1 in s: print("C1 cols: already")
elif anchorC1 in s: s = s.replace(anchorC1, newC1, 1); print("C1 cols: added")
else: print("C1 cols: ANCHOR NOT FOUND")

# C2) explorer body row: add gene-download button cell
anchorC2 = "<td>${r.ci}</td></tr>`"
newC2 = "<td>${r.ci}</td><td><button class=\"geneBtn\" data-gk=\"${geneKeyFor(r)}\" title=\"Download the genes this model used\">genes</button></td></tr>`"
if "geneBtn" in s and anchorC2 not in s: print("C2 row: already")
elif anchorC2 in s: s = s.replace(anchorC2, newC2, 1); print("C2 row: added")
else: print("C2 row: ANCHOR NOT FOUND")

# F) helper + download functions (append before init())
if "function geneKeyFor(" in s:
    print("F funcs: already present")
else:
    funcs = r'''function geneKeyFor(r){
  var m=r.model, km=(m.indexOf('glm_')===0)?m:m.split('_')[0];
  return (useSmoking()?'smoke':'full')+'|'+[r.set,r.cohort,r.datatype,r.featureset,r.subset,km,'v'+r.varcut].join('_');
}
function downloadGenesByKey(key){
  var e=(typeof GENELISTS!=='undefined')?GENELISTS[key]:null;
  if(!e||!e.g||!e.g.length){ alert('No downloadable gene list for this model.\n\nIt either uses a non-gene feature set (clinics / module20 / module20_clinics) or the model selected no genes. Full gene lists for all gene-based models are in the companion file model_gene_lists_202609.csv.gz.'); return; }
  var vlab=e.m==='glm_coef'?'glm_coefficient':(e.m==='xgb_gain'?'xgb_gain':'rf_importance');
  var lines=['rank,gene,'+vlab];
  e.g.forEach(function(gv,i){ lines.push((i+1)+','+gv[0]+','+gv[1]); });
  var blob=new Blob([lines.join('\n')],{type:'text/csv'});
  var a=document.createElement('a'); a.href=URL.createObjectURL(blob); a.download='genes_'+key.replace('|','_')+'.csv'; a.click(); URL.revokeObjectURL(a.href);
}
init();
'''
    s = s.replace("init();\n</script>", funcs + "</script>", 1)
    print("F funcs: inserted" if "function geneKeyFor(" in s else "F funcs: INSERT FAILED")

# G) delegated click listener on the explorer table (added once at init)
anchorG = "document.getElementById('fDownload').addEventListener('click', downloadExplorer);"
newG = anchorG + "\n  document.getElementById('explorerTable').addEventListener('click', function(ev){ var b=ev.target.closest && ev.target.closest('.geneBtn'); if(b){ ev.stopPropagation(); downloadGenesByKey(b.dataset.gk); } });"
if "explorerTable').addEventListener('click'" in s: print("G listener: already")
elif anchorG in s: s = s.replace(anchorG, newG, 1); print("G listener: added")
else: print("G listener: ANCHOR NOT FOUND")

if s != orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")

# ---------- builder: inject __GENELISTS__ ----------
r = io.open(B, encoding="utf-8").read(); rorig = r
anchorR = 'tpl <- sub1(tpl, "__VENN_DATA__", venn_json)'
inj = anchorR + '\ngenelists_json <- paste(readLines("model_gene_lists_compact.json", warn = FALSE), collapse = "")\ntpl <- sub1(tpl, "__GENELISTS__", genelists_json)'
if "__GENELISTS__" in r: print("builder: already injects")
elif anchorR in r: r = r.replace(anchorR, inj, 1); print("builder: injection added")
else: print("builder: ANCHOR NOT FOUND")
if r != rorig: io.open(B,"w",encoding="utf-8").write(r); print("WROTE builder")

for probe in ["const GENELISTS", "{key:'__genes'", "geneBtn", "function geneKeyFor(", "explorerTable').addEventListener('click'"]:
    print(("OK  " if probe in s else "MISS ")+probe)
print(("OK  __GENELISTS__ in builder" if "__GENELISTS__" in r else "MISS builder inject"))
