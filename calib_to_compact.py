#!/usr/bin/env python3
# roc_calib_202609.csv -> roc_calib_compact.json  {"pop|key":[auc,[fpr...],[tpr...]]}
import csv, json, sys
src = "roc_calib_202609.csv"; dst = "roc_calib_compact.json"
out = {}
with open(src, newline="") as f:
    for row in csv.DictReader(f):
        k = row["pop"] + "|" + row["key"]
        fpr = [round(float(x), 4) for x in row["fpr"].split("|")]
        tpr = [round(float(x), 4) for x in row["tpr"].split("|")]
        out[k] = [round(float(row["auc_calib"]), 4), fpr, tpr]
with open(dst, "w") as f:
    json.dump(out, f, separators=(",", ":"))
print("wrote %s with %d configs" % (dst, len(out)))
