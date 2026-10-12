import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"; B="codes_20260713/make_dashboard_fsmethods_cmp.R"
s=io.open(T,encoding="utf-8").read(); orig=s
co=io.open("co2_funcs.js",encoding="utf-8").read()

# A) const FOLDAVG after const GENELISTS
if "const FOLDAVG" in s: print("FOLDAVG const: already")
elif "const GENELISTS = __GENELISTS__;" in s:
    s=s.replace("const GENELISTS = __GENELISTS__;","const GENELISTS = __GENELISTS__;\nconst FOLDAVG = __FOLDAVG__;",1); print("FOLDAVG const: added")
else: print("FOLDAVG const: ANCHOR NOT FOUND")

# B) Basis toggle before the Curve selector (inside cutoff section)
if 'id="coMode"' in s: print("coMode select: already")
elif '<label>Curve <select id="coCurve"></select></label>' in s:
    sel='<label>Basis <select id="coMode"><option value="pooled">Pooled (probability cutoff)</option><option value="vavg">Per-fold vertical-averaged</option></select></label>\n    <label>Curve <select id="coCurve"></select></label>'
    s=s.replace('<label>Curve <select id="coCurve"></select></label>', sel, 1); print("coMode select: added")
else: print("coMode select: ANCHOR NOT FOUND")

# C) relabel slider span
if 'id="coSlideLab"' in s: print("slider label: already")
elif '<label>Cutoff <input type="range" id="coSlider"' in s:
    s=s.replace('<label>Cutoff <input type="range" id="coSlider"','<label><span id="coSlideLab">Cutoff (probability)</span> <input type="range" id="coSlider"',1); print("slider label: added")
else: print("slider label: ANCHOR NOT FOUND")

# D) override renderCutoff (append before init so it wins)
if "function coMode()" in s: print("override funcs: already")
else:
    s=s.replace("init();\n</script>", co+"init();\n</script>",1)
    print("override funcs: inserted" if "function coMode()" in s else "override funcs: INSERT FAILED")

# E) wire coMode change
if "getElementById('coMode').addEventListener" in s: print("coMode wiring: already")
elif "document.getElementById('coSlider').addEventListener('input', renderCutoff);" in s:
    s=s.replace("document.getElementById('coSlider').addEventListener('input', renderCutoff);",
                "document.getElementById('coSlider').addEventListener('input', renderCutoff);\n  document.getElementById('coMode').addEventListener('change', renderCutoff);",1); print("coMode wiring: added")
else: print("coMode wiring: ANCHOR NOT FOUND")

if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")

# F) builder: inject __FOLDAVG__
r=io.open(B,encoding="utf-8").read(); rorig=r
if "__FOLDAVG__" in r: print("builder: already")
elif 'tpl <- sub1(tpl, "__GENELISTS__", genelists_json)' in r:
    add='tpl <- sub1(tpl, "__GENELISTS__", genelists_json)\nfoldavg_json <- paste(readLines("roc_foldavg_compact.json", warn = FALSE), collapse = "")\ntpl <- sub1(tpl, "__FOLDAVG__", foldavg_json)'
    r=r.replace('tpl <- sub1(tpl, "__GENELISTS__", genelists_json)', add, 1); print("builder: injection added")
else: print("builder: ANCHOR NOT FOUND")
if r!=rorig: io.open(B,"w",encoding="utf-8").write(r); print("WROTE builder")

for probe in ["const FOLDAVG",'id="coMode"','id="coSlideLab"',"function coMode()","getElementById('coMode').addEventListener"]:
    print(("OK  " if probe in s else "MISS ")+probe)
print(("OK  __FOLDAVG__ builder" if "__FOLDAVG__" in r else "MISS builder"))
