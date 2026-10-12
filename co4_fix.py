import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"
s=io.open(T,encoding="utf-8").read(); orig=s
reps=[
 ("'at cutoff <b>'+p.toFixed(2)+'</b>: sensitivity <b>'+sens.toFixed(3)+'</b>, specificity <b>'+spec.toFixed(3)+'</b>'",
  "'Sensitivity <b>'+sens.toFixed(3)+'</b> &nbsp;|&nbsp; Specificity <b>'+spec.toFixed(3)+'</b> &nbsp;&middot; at probability cutoff '+p.toFixed(2)"),
 ("'at specificity <b>'+spec.toFixed(2)+'</b>: sensitivity <b>'+sens.toFixed(3)+'</b> (FPR '+fprv.toFixed(2)+')'",
  "'Sensitivity <b>'+sens.toFixed(3)+'</b> &nbsp;|&nbsp; Specificity <b>'+spec.toFixed(2)+'</b> &nbsp;&middot; FPR '+fprv.toFixed(2)"),
 # add both to marker hovers
 ("var cur={x:[fpr[k]],y:[tp[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'}};",
  "var cur={x:[fpr[k]],y:[tp[k]],mode:'markers',name:'cutoff '+p.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'sensitivity '+sens.toFixed(3)+'<br>specificity '+spec.toFixed(3)+'<extra></extra>'};"),
 ("var cur={x:[fprv],y:[sens],mode:'markers',name:'spec '+spec.toFixed(2),marker:{size:12,color:'#e67e22'}};",
  "var cur={x:[fprv],y:[sens],mode:'markers',name:'spec '+spec.toFixed(2),marker:{size:12,color:'#e67e22'},hovertemplate:'sensitivity '+sens.toFixed(3)+'<br>specificity '+spec.toFixed(2)+'<extra></extra>'};"),
]
for old,new in reps:
    c=s.count(old); s=s.replace(old,new); print(("OK %d"%c if c else "MISS"), old[:48])
if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE")
