import glob, re, os, statistics as st
pat = re.compile(r"^(prev|current)_(wholeset_coding|wholeset|stage1_coding|stage1|intersection_coding|intersection)_(normalized|fpkm)_(clinics|fpkm_clinics|fpkm|module20_clinics|module20)_(12vs0|3vs0|123vs0)_glm_(mad_nosexmt|mad|cor)_v(\d+)_cvSummary\.csv$")
out = ["\t".join("set cohort datatype dataname subset fs varcut n meanAUC sdAUC bootAUC".split())]
for f in sorted(glob.glob("res/*_cvSummary.csv")):
    m = pat.match(os.path.basename(f))
    if not m:
        continue
    aucs, boots = [], {}
    with open(f) as fh:
        hdr = fh.readline().strip().split(",")
        try:
            ia, ib = hdr.index("AUC"), hdr.index("bootid")
        except ValueError:
            continue
        for line in fh:
            p = line.strip().split(",")
            if len(p) <= max(ia, ib):
                continue
            try:
                a = float(p[ia])
            except ValueError:
                continue
            aucs.append(a)
            boots.setdefault(p[ib], []).append(a)
    if not aucs:
        continue
    bm = [st.mean(v) for v in boots.values()]
    out.append("\t".join(list(m.groups()) + [str(len(aucs)), f"{st.mean(aucs):.4f}", f"{st.pstdev(aucs):.4f}", f"{st.mean(bm):.4f}"]))
print("\n".join(out))
