import io, os, sys
os.chdir("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")

def edit(path, repls):
    with io.open(path, encoding="utf-8") as f: s = f.read()
    orig = s; log = []
    for old, new, allf in repls:
        if new.split("glm_natcommad")[0] and "glm_natcommad" in new and s.count("glm_natcommad") and old not in s and new.replace(old,"") in s:
            log.append("  already wired (skip)"); continue
        if old not in s:
            log.append("  ANCHOR NOT FOUND: %r" % old[:40]); continue
        n = s.count(old) if allf else 1
        s = s.replace(old, new) if allf else s.replace(old, new, 1)
        log.append("  replaced %dx: %r" % (n, old[:40]))
    if s != orig:
        with io.open(path, "w", encoding="utf-8") as f: f.write(s)
    print(path, "CHANGED" if s != orig else "unchanged")
    for l in log: print(l)

builder_old = 'mod_lv <- c("glm_mad","glm_cor","glm_mad_nosexmt","rf_mad","xgb_mad")'
builder_new = 'mod_lv <- c("glm_mad","glm_cor","glm_mad_nosexmt","rf_mad","xgb_mad","glm_natcommad")'
for b in ["codes_20260713/make_dashboard_fsmethods_cmp.R",
          "codes_20260713/make_dashboard_fsmethods.R",
          "codes_20260713/make_dashboard_fsmethods_smoking.R"]:
    edit(b, [(builder_old, builder_new, False)])

tpl_repls = [
    ("xgb_mad:'XGB+MAD'}", "xgb_mad:'XGB+MAD', glm_natcommad:'GLM+ncMAD'}", False),
    ("xgb_mad:'#d68910'",  "xgb_mad:'#d68910', glm_natcommad:'#17a2b8'",     True),
]
for t in ["codes_20260713/dashboard_template_fsmethods_cmp.html",
          "codes_20260713/dashboard_template_fsmethods.html"]:
    edit(t, tpl_repls)
