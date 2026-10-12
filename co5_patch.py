import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"; B="codes_20260713/make_dashboard_fsmethods_cmp.R"
s=io.open(T,encoding="utf-8").read(); orig=s
co=io.open("co5_funcs.js",encoding="utf-8").read()

# A) const INSAMPLE after const FOLDAVG
if "const INSAMPLE" in s: print("INSAMPLE const: already")
elif "const FOLDAVG = __FOLDAVG__;" in s:
    s=s.replace("const FOLDAVG = __FOLDAVG__;","const FOLDAVG = __FOLDAVG__;\nconst INSAMPLE = __INSAMPLE__;",1); print("INSAMPLE const: added")
else: print("INSAMPLE const: ANCHOR NOT FOUND")

# B) third option in coMode select
if 'value="insample"' in s: print("insample option: already")
elif '<option value="vavg">Per-fold vertical-averaged</option>' in s:
    s=s.replace('<option value="vavg">Per-fold vertical-averaged</option>',
                '<option value="vavg">Per-fold vertical-averaged</option><option value="insample">Full-data model (in-sample)</option>',1); print("insample option: added")
else: print("insample option: ANCHOR NOT FOUND")

# C) override renderCutoff v3 (append before init so it wins)
if "three bases -- pooled CV" in s: print("renderCutoff v3: already")
else:
    s=s.replace("init();\n</script>", co+"init();\n</script>",1)
    print("renderCutoff v3: inserted" if "three bases -- pooled CV" in s else "renderCutoff v3: INSERT FAILED")

if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")

# D) builder inject __INSAMPLE__
r=io.open(B,encoding="utf-8").read(); rorig=r
if "__INSAMPLE__" in r: print("builder: already")
elif 'tpl <- sub1(tpl, "__FOLDAVG__", foldavg_json)' in r:
    add='tpl <- sub1(tpl, "__FOLDAVG__", foldavg_json)\ninsample_json <- paste(readLines("roc_insample_compact.json", warn = FALSE), collapse = "")\ntpl <- sub1(tpl, "__INSAMPLE__", insample_json)'
    r=r.replace('tpl <- sub1(tpl, "__FOLDAVG__", foldavg_json)', add, 1); print("builder: injection added")
else: print("builder: ANCHOR NOT FOUND")
if r!=rorig: io.open(B,"w",encoding="utf-8").write(r); print("WROTE builder")

for probe in ["const INSAMPLE",'value="insample"',"three bases -- pooled CV"]:
    print(("OK  " if probe in s else "MISS ")+probe)
print(("OK  __INSAMPLE__ builder" if "__INSAMPLE__" in r else "MISS builder"))
