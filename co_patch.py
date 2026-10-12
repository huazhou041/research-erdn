import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T = "codes_20260713/dashboard_template_fsmethods_cmp.html"
s = io.open(T, encoding="utf-8").read(); orig = s
co = io.open("co_funcs.js", encoding="utf-8").read()

# 1) nav link after ROC curves
if '#cutoff' not in s:
    a='<a href="#roc">ROC curves</a>'
    if a in s: s=s.replace(a, a+'\n  <a href="#cutoff">Cutoff</a>',1); print("nav: added")
    else: print("nav: ANCHOR NOT FOUND")
else: print("nav: already")

# 2) section before genevenn
if 'id="cutoff"' not in s:
    sec='''<section id="cutoff">
  <h2>Operating point &amp; optimal cutoff</h2>
  <p class="sub">Pick one ROC curve (from the filters above), then slide the probability cutoff to read sensitivity and specificity at that threshold. The red star marks the <b>Youden-optimal</b> cutoff (maximizes sensitivity + specificity − 1). These are pooled cross-validated predictions, so this is an illustrative operating point, not a validated deployment threshold. Follows the full / ever-smoker toggle.</p>
  <div class="controls">
    <label>Curve <select id="coCurve"></select></label>
    <label>Cutoff <input type="range" id="coSlider" min="0" max="1" step="0.01" value="0.5" style="vertical-align:middle;width:200px"> <span id="coCut">0.50</span></label>
    <span id="coReadout" class="sub" style="margin:0"></span>
  </div>
  <div id="coChart" style="height:430px"></div>
  <div id="coOpt" class="sub" style="margin-top:8px"></div>
</section>
<section id="genevenn">'''
    if '<section id="genevenn">' in s: s=s.replace('<section id="genevenn">', sec, 1); print("section: added")
    else: print("section: ANCHOR NOT FOUND")
else: print("section: already")

# 3) POP_SECTIONS add cutoff
if "'roc','cutoff'" in s: print("POP_SECTIONS: already")
elif "'roc','explorer','natcommad']" in s:
    s=s.replace("'roc','explorer','natcommad']","'roc','cutoff','explorer','natcommad']",1); print("POP_SECTIONS: added")
else: print("POP_SECTIONS: ANCHOR NOT FOUND")

# 4) functions before init()
if "function renderCutoff(" in s: print("funcs: already")
else:
    s=s.replace("init();\n</script>", co+"init();\n</script>",1)
    print("funcs: inserted" if "function renderCutoff(" in s else "funcs: INSERT FAILED")

# 5) renderAll: also populate cutoff
if "renderNatcommad(); populateCutoff();" in s: print("renderAll: already")
elif "renderExplorer(); renderNatcommad();" in s:
    s=s.replace("renderExplorer(); renderNatcommad();","renderExplorer(); renderNatcommad(); populateCutoff();",1); print("renderAll: added")
else: print("renderAll: ANCHOR NOT FOUND")

# 6) init wiring: slider + select + ROC-filter repopulate + plotly-ready
if "getElementById('coSlider')" not in s:
    anchor="  document.getElementById('rocReset').addEventListener"
    add=("  document.getElementById('coSlider').addEventListener('input', renderCutoff);\n"
         "  document.getElementById('coCurve').addEventListener('change', renderCutoff);\n"
         "  ['rocSet','rocCohort','rocDatatype','rocFeatureset','rocSubset','rocModel','rocVarcut'].forEach(function(id){ document.getElementById(id).addEventListener('change', populateCutoff); });\n"
         "  document.getElementById('rocSearch').addEventListener('input', populateCutoff);\n")
    if anchor in s: s=s.replace(anchor, add+anchor, 1); print("init wiring: added")
    else: print("init wiring: ANCHOR NOT FOUND")
else: print("init wiring: already")

# 7) plotly-ready: render cutoff when Plotly becomes available
if "renderGeneCorr(); populateCutoff();" in s: print("plotly-ready: already")
elif "renderROC(); renderGeneCorr(); });" in s:
    s=s.replace("renderROC(); renderGeneCorr(); });","renderROC(); renderGeneCorr(); populateCutoff(); });",1); print("plotly-ready: added")
else: print("plotly-ready: ANCHOR NOT FOUND (will still render via renderAll)")

if s != orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")
for probe in ['id="cutoff"','function renderCutoff(',"'roc','cutoff'","getElementById('coSlider')",'#cutoff']:
    print(("OK  " if probe in s else "MISS ")+probe)
