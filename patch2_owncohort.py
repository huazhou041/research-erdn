import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T = "codes_20260713/dashboard_template_fsmethods_cmp.html"
s = io.open(T, encoding="utf-8").read(); orig = s
ov = io.open("nc_override.js", encoding="utf-8").read()

# A) table header: add own-MAD + delta columns (replace the mean/sd/n tail)
old_h = "<th>top-N genes</th><th>mean AUC</th><th>sd</th><th>n</th>"
new_h = "<th>top-N genes</th><th>AUC own-cohort MAD</th><th>AUC nature_com MAD</th><th>&Delta; (nc&minus;own)</th><th>n</th>"
if old_h in s: s = s.replace(old_h, new_h, 1); print("header: replaced")
elif new_h in s: print("header: already updated")
else: print("header: ANCHOR NOT FOUND")

# B) add a delta sort option
old_o = '<option value="cohort">cohort</option></select>'
new_o = '<option value="cohort">cohort</option><option value="delta">&Delta; (nc&minus;own)</option></select>'
if 'value="delta"' in s: print("sort opt: already present")
elif old_o in s: s = s.replace(old_o, new_o, 1); print("sort opt: added")
else: print("sort opt: ANCHOR NOT FOUND")

# C) description tweak
old_p = "Filter below, or use Download CSV to export the shown rows."
new_p = "Each row also shows the AUC from that cohort&rsquo;s OWN top-MAD genes (GLM+MAD) and the difference, for a direct comparison. Filter below, or Download CSV."
if old_p in s: s = s.replace(old_p, new_p, 1); print("desc: updated")
else: print("desc: (anchor not found, skipping)")

# D) append override functions just before init() so they win (hoisting: last decl wins)
if "function ncRowsWithOwn()" in s:
    print("override JS: already present")
else:
    s = s.replace("init();\n</script>", ov + "\ninit();\n</script>", 1)
    print("override JS: inserted" if ov.split('(')[0] in s else "override JS: INSERT FAILED")

if s != orig:
    io.open(T, "w", encoding="utf-8").write(s); print("WROTE", T)
else:
    print("no change")
for probe in [new_h, 'value="delta"', "function ncRowsWithOwn()", "auc_own_cohort_mad"]:
    print(("OK  " if probe in s else "MISS ") + probe[:40])
