#!/usr/bin/env python3
# Build a filtered variant of the full-vs-smoking dashboard builder:
#   - keep datatype == "normalized" only (drop fpkm datatype everywhere)
#   - drop task subset "3vs0" everywhere
# Writes a NEW builder (…_trim.R) -> NEW html, leaving the main report untouched.
import io, sys
SRC = "codes_20260713/make_dashboard_fsmethods_cmp.R"
DST = "codes_20260713/make_dashboard_fsmethods_cmp_trim.R"
s = io.open(SRC, encoding="utf-8").read()
orig = s

# 1) insert a read.csv wrapper right after the library load (shadows base read.csv
#    for every CSV the builder reads; filters only where the columns exist)
anchor = "suppressPackageStartupMessages({ library(jsonlite) })"
wrapper = anchor + """

## ---- TRIM VARIANT: normalized-only, 3vs0 removed -------------------------
## Shadow read.csv so EVERY data frame with a datatype/subset column is filtered
## identically at the source. datatype=="normalized" only; subset!="3vs0".
.read_csv_base <- utils::read.csv
read.csv <- function(file, ...) {
  df <- .read_csv_base(file, ...)
  if (is.data.frame(df)) {
    if ("datatype" %in% names(df)) df <- df[df$datatype == "normalized", , drop = FALSE]
    if ("subset"   %in% names(df)) df <- df[df$subset != "3vs0",       , drop = FALSE]
  }
  df
}
## -------------------------------------------------------------------------"""
assert anchor in s, "library anchor not found"
assert "TRIM VARIANT" not in s
s = s.replace(anchor, wrapper, 1)

# 2) trim the menu vectors (drive the UI selectors + META)
s = s.replace('dm_lv  <- c("fpkm","normalized")', 'dm_lv  <- c("normalized")', 1)
s = s.replace('cmp_lv <- c("123vs0","12vs0","3vs0")', 'cmp_lv <- c("123vs0","12vs0")', 1)

# 3) mk_cmp_table hardcoded task list
s = s.replace('  tasks <- c("123vs0","12vs0","3vs0")', '  tasks <- c("123vs0","12vs0")', 1)

# 4) drop the now-stale 3vs0 clause in the cmp-table caption (if present)
s = s.replace(
  ", and <span style=\"color:#b8860b\">&#9888;</span> 3vs0 cells have &le;12 recurrence events",
  "", 1)

# 5) new output filename (leave the main report untouched)
s = s.replace(
  'out_file <- "htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking.html"',
  'out_file <- "htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking_normOnly_no3vs0.html"', 1)

changed = sum(1 for a,b in [
  ('dm_lv  <- c("normalized")', s),
  ('cmp_lv <- c("123vs0","12vs0")', s),
  ('  tasks <- c("123vs0","12vs0")', s),
  ('normOnly_no3vs0.html', s),
  ('TRIM VARIANT', s),
] if a in b)
io.open(DST, "w", encoding="utf-8").write(s)
print("wrote", DST, "| edits applied:", changed, "/5")
for probe in ['dm_lv  <- c("normalized")','cmp_lv <- c("123vs0","12vs0")',
              '  tasks <- c("123vs0","12vs0")','normOnly_no3vs0.html','TRIM VARIANT']:
    print(("OK   " if probe in s else "MISS ")+probe)
