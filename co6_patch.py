import io, os
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
T="codes_20260713/dashboard_template_fsmethods_cmp.html"; B="codes_20260713/make_dashboard_fsmethods_cmp.R"
s=io.open(T,encoding="utf-8").read(); orig=s
co=io.open("co6_funcs.js",encoding="utf-8").read()

# A) const CALIB after const INSAMPLE
if "const CALIB" in s: print("CALIB const: already")
elif "const INSAMPLE = __INSAMPLE__;" in s:
    s=s.replace("const INSAMPLE = __INSAMPLE__;","const INSAMPLE = __INSAMPLE__;\nconst CALIB = __CALIB__;",1); print("CALIB const: added")
else: print("CALIB const: ANCHOR NOT FOUND")

# B) reorder + add calib option as default (first)
old_opts='<option value="pooled">Pooled (probability cutoff)</option><option value="vavg">Per-fold vertical-averaged</option><option value="insample">Full-data model (in-sample)</option>'
new_opts=('<option value="calib" selected>Calibrated CV (deployable, recommended)</option>'
          '<option value="vavg">Per-fold vertical-averaged (honest trade-off)</option>'
          '<option value="pooled">Pooled CV (raw probability)</option>'
          '<option value="insample">Full-data model (in-sample, optimistic)</option>')
if 'value="calib"' in s: print("options: already")
elif old_opts in s: s=s.replace(old_opts,new_opts,1); print("options: reordered + calib default")
else: print("options: ANCHOR NOT FOUND")

# C) override renderCutoff v4
if "four bases -- calibrated CV" in s: print("renderCutoff v4: already")
else:
    s=s.replace("init();\n</script>", co+"init();\n</script>",1)
    print("renderCutoff v4: inserted" if "four bases -- calibrated CV" in s else "renderCutoff v4: INSERT FAILED")

if s!=orig: io.open(T,"w",encoding="utf-8").write(s); print("WROTE template")

# D) builder inject __CALIB__
r=io.open(B,encoding="utf-8").read(); rorig=r
if "__CALIB__" in r: print("builder: already")
elif 'tpl <- sub1(tpl, "__INSAMPLE__", insample_json)' in r:
    add='tpl <- sub1(tpl, "__INSAMPLE__", insample_json)\ncalib_json <- paste(readLines("roc_calib_compact.json", warn = FALSE), collapse = "")\ntpl <- sub1(tpl, "__CALIB__", calib_json)'
    r=r.replace('tpl <- sub1(tpl, "__INSAMPLE__", insample_json)', add, 1); print("builder: injection added")
else: print("builder: ANCHOR NOT FOUND")
if r!=rorig: io.open(B,"w",encoding="utf-8").write(r); print("WROTE builder")

for probe in ["const CALIB",'value="calib" selected',"four bases -- calibrated CV"]:
    print(("OK  " if probe in s else "MISS ")+probe)
print(("OK  __CALIB__ builder" if "__CALIB__" in r else "MISS builder"))
