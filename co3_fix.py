import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"
s=io.open(T,encoding="utf-8").read(); orig=s
old="margin:{t:20},legend:{orientation:'h'}"
new="margin:{t:30,b:70,l:60,r:20},legend:{orientation:'h',y:1.15,yanchor:'bottom',x:0}"
n=s.count(old)
s=s.replace(old,new)
print("replaced", n, "cutoff-chart layouts")
if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")
print("remaining old:", s.count(old), " new present:", s.count(new))
