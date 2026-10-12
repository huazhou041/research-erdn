import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"
s=io.open(T,encoding="utf-8").read(); orig=s
if "coSlider').addEventListener('input'" in s:
    print("listeners: already present")
else:
    anchor="  document.getElementById('rocReset').addEventListener"
    add=("  document.getElementById('coSlider').addEventListener('input', renderCutoff);\n"
         "  document.getElementById('coCurve').addEventListener('change', renderCutoff);\n"
         "  ['rocSet','rocCohort','rocDatatype','rocFeatureset','rocSubset','rocModel','rocVarcut'].forEach(function(id){ document.getElementById(id).addEventListener('change', populateCutoff); });\n"
         "  document.getElementById('rocSearch').addEventListener('input', populateCutoff);\n")
    if anchor in s: s=s.replace(anchor, add+anchor, 1); print("listeners: ADDED")
    else: print("listeners: ANCHOR NOT FOUND")
if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE")
print("verify listener present:", "coSlider').addEventListener('input'" in s)
