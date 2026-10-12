import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T = "codes_20260713/dashboard_template_fsmethods_cmp.html"
with io.open(T, encoding="utf-8") as f: s = f.read()
orig = s

if 'id="natcommad"' in s:
    print("section already present -- skipping HTML/nav insert")
else:
    # 1) nav link before Explorer
    s = s.replace('<a href="#explorer">Explorer</a>',
                  '<a href="#natcommad">ncMAD table</a>\n  <a href="#explorer">Explorer</a>', 1)
    # 2) section after venn (which ends right before the explorer section)
    section = '''</section>
<section id="natcommad">
  <h2>nature_com top-MAD (GLM+ncMAD): AUC results table</h2>
  <p class="sub">GLM trained on the <b>fixed nature_com top-MAD gene list</b> (the same ranking shown in the Venn above), transferred to our prev/current data &mdash; the &ldquo;start from the nature_com paper&rsquo;s top genes&rdquo; arm. Gene-based feature sets only (fpkm, fpkm_clinics). AUC = 100&times;10 fold-averaged, fixed direction. This table follows the FULL vs EVER-SMOKER toggle. Filter below, or use Download CSV to export the shown rows.</p>
  <div class="controls">
    <label>Set <select id="ncSet"><option value="all">both</option><option value="prev">prev</option><option value="current" selected>current</option></select></label>
    <label>Task <select id="ncTask"><option value="all">all</option><option value="123vs0" selected>123vs0</option><option value="12vs0">12vs0</option><option value="3vs0">3vs0</option></select></label>
    <label>Feature set <select id="ncFs"><option value="all" selected>both</option><option value="fpkm">fpkm</option><option value="fpkm_clinics">fpkm_clinics</option></select></label>
    <label>Sort <select id="ncSort"><option value="mean" selected>AUC (high&rarr;low)</option><option value="cohort">cohort</option></select></label>
    <span id="ncSummary" class="sub" style="margin:0"></span>
    <button id="ncDownload">Download CSV</button>
  </div>
  <div class="tablewrap" style="max-height:480px;overflow:auto"><table id="ncTable"><thead><tr><th>set</th><th>cohort</th><th>datatype</th><th>feature set</th><th>task</th><th>top-N genes</th><th>mean AUC</th><th>sd</th><th>n</th></tr></thead><tbody></tbody></table></div>
</section>
<section id="explorer">'''
    s = s.replace('</section>\n<section id="explorer">', section, 1)

# 3) JS functions (before renderAll)
if 'function renderNatcommad' not in s:
    js = '''function ncRows() {
  const setF=document.getElementById('ncSet').value, taskF=document.getElementById('ncTask').value,
        fsF=document.getElementById('ncFs').value, sortF=document.getElementById('ncSort').value;
  let rows = LONG.filter(r => r.model==='glm_natcommad'
     && (setF==='all'||r.set===setF) && (taskF==='all'||r.subset===taskF) && (fsF==='all'||r.featureset===fsF));
  if (sortF==='mean') rows.sort((a,b)=> (b.mean==null?-1:b.mean)-(a.mean==null?-1:a.mean));
  else rows.sort((a,b)=> String(a.cohort).localeCompare(String(b.cohort)) || (b.mean==null?-1:b.mean)-(a.mean==null?-1:a.mean));
  return rows;
}
function renderNatcommad() {
  if (!document.getElementById('ncTable')) return;
  const rows = ncRows();
  document.querySelector('#ncTable tbody').innerHTML = rows.map(r =>
    `<tr><td>${r.set}</td><td>${r.cohort}</td><td>${r.datatype}</td><td>${r.featureset}</td><td>${r.subset}</td><td>${r.varcut}</td><td>${r.mean==null?'':r.mean.toFixed(3)}</td><td>${r.sd==null?'':r.sd.toFixed(3)}</td><td>${r.n==null?'':r.n}</td></tr>`).join('');
  const best = rows.length ? rows.reduce((m,r)=> ((r.mean==null?-1:r.mean)>(m.mean==null?-1:m.mean)?r:m)) : null;
  document.getElementById('ncSummary').textContent = `${rows.length} configs \\u00b7 ${useSmoking()?'ever-smoker':'full'} cohort` + (best? ` \\u00b7 best AUC ${best.mean.toFixed(3)} (${best.set}/${best.cohort}/${best.featureset}/top-${best.varcut})`:'');
}
function downloadNatcommad() {
  const rows = ncRows();
  const cols = ['set','cohort','datatype','featureset','subset','varcut','mean','sd','n'];
  const hdr = ['set','cohort','datatype','featureset','task','top_n_genes','mean_auc','sd_auc','n_auc'];
  const csv = [hdr.join(',')].concat(rows.map(r => cols.map(c => `"${String(r[c]==null?'':r[c]).replace(/"/g,'""')}"`).join(','))).join('\\n');
  const blob = new Blob([csv],{type:'text/csv'});
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob);
  a.download = `natcommad_${useSmoking()?'eversmoker':'full'}.csv`; a.click(); URL.revokeObjectURL(a.href);
}
function renderAll() {'''
    s = s.replace('function renderAll() {', js, 1)

# 4) call renderNatcommad inside renderAll
if 'renderExplorer(); renderNatcommad();' not in s:
    s = s.replace('explorerPage=0; renderExplorer();\n}', 'explorerPage=0; renderExplorer(); renderNatcommad();\n}', 1)

# 5) add to POP_SECTIONS
if "'explorer','natcommad'" not in s:
    s = s.replace("'roc','explorer'];", "'roc','explorer','natcommad'];", 1)

# 6) wire listeners in init
if 'ncDownload' not in s.split('function downloadNatcommad')[-1]:  # only in the wiring area
    pass
if "getElementById('ncDownload').addEventListener" not in s:
    s = s.replace("  document.getElementById('vnCohort').addEventListener('change', renderVenn);\n  renderVenn();",
                  "  document.getElementById('vnCohort').addEventListener('change', renderVenn);\n  renderVenn();\n  ['ncSet','ncTask','ncFs','ncSort'].forEach(id => document.getElementById(id).addEventListener('change', renderNatcommad));\n  document.getElementById('ncDownload').addEventListener('click', downloadNatcommad);", 1)

if s != orig:
    with io.open(T, "w", encoding="utf-8") as f: f.write(s)
    print("PATCHED", T)
else:
    print("no change (already patched)")
# verify anchors took
for probe in ['id="natcommad"', 'function renderNatcommad', "'roc','explorer','natcommad'", 'renderExplorer(); renderNatcommad();', "getElementById('ncDownload')"]:
    print(("OK  " if probe in s else "MISS ") + probe)
