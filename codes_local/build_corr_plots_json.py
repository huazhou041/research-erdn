import base64, json, os, re

# CONFIG (override via env). CORR_PNG_DIR must hold the correlation-heatmap PNGs;
# the original ephemeral session scratchpad no longer exists, so set this to a
# persistent location holding regenerated inputs before re-running.
PROJ = os.environ.get("ERDN_PROJ", ".")
SCRATCH = os.environ.get("CORR_PNG_DIR", os.path.join(PROJ, "provenance_inputs", "corrpng"))
OUT = os.environ.get("CORR_PLOTS_OUT", os.path.join(PROJ, "htmls", "corr_plots.json"))

records = []

def b64_of(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode("ascii")

# prev_qc dropped from the dashboard per instruction -- internal prev-only check,
# not a comparison against nature_com/current. PNGs left in scratch, just not embedded.

# current_naturecom: corr_natureCom_{Normal|Tumor}_vs_{prev|current}_{full|topMAD}_{pearson|spear}[_diagonal].png
pat2 = re.compile(r"^corr_natureCom_(Normal|Tumor)_vs_(prev|current)_(full|topMAD)_(pearson|spear)(_diagonal)?$")
d2 = os.path.join(SCRATCH, "current_naturecom")
for fn in sorted(os.listdir(d2)):
    if not fn.endswith(".png"): continue
    stem = fn[:-4]
    m = pat2.match(stem)
    if not m:
        print("NO MATCH (current_naturecom):", fn); continue
    tissue, compare, geneset, method, diag = m.group(1), m.group(2), m.group(3), m.group(4), bool(m.group(5))
    if tissue == "Tumor":
        continue  # dropped per instruction -- Normal only, matches the rest of this dashboard's nature_com scope
    # folded directly into "source" -- one of the 3 named comparisons, not a sub-selector
    src = "prev_vs_naturecom" if compare == "prev" else "current_vs_naturecom"
    records.append({
        "source": src, "filter": None,
        "geneset": geneset, "method": method, "diagonal": diag,
        "label": f"{compare} vs nature_com (Normal) - {geneset} genes - {method} - {'diagonal' if diag else 'raw'}",
        "img": b64_of(os.path.join(d2, fn))
    })

# prev_vs_current (direct): corr_prev_vs_current_{full|topMAD}_{pearson|spear}[_diagonal].png
pat3 = re.compile(r"^corr_prev_vs_current_(full|topMAD)_(pearson|spear)(_diagonal)?$")
d3 = os.path.join(SCRATCH, "prev_vs_current")
if os.path.isdir(d3):
    for fn in sorted(os.listdir(d3)):
        if not fn.endswith(".png"): continue
        stem = fn[:-4]
        m = pat3.match(stem)
        if not m:
            print("NO MATCH (prev_vs_current):", fn); continue
        geneset, method, diag = m.group(1), m.group(2), bool(m.group(3))
        records.append({
            "source": "prev_vs_current", "filter": None, "tissue": None, "compare": None,
            "geneset": geneset, "method": method, "diagonal": diag,
            "label": f"prev vs current (direct) - {geneset} genes - {method} - {'diagonal' if diag else 'raw'}",
            "img": b64_of(os.path.join(d3, fn))
        })
else:
    print("prev_vs_current dir not found yet, skipping")

with open(OUT, "w") as fh:
    json.dump(records, fh)

print("wrote", OUT, "with", len(records), "records")
total_bytes = sum(len(r["img"]) for r in records)
print("total base64 payload: {:.2f} MB".format(total_bytes / 1e6))
