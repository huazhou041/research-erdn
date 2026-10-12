cohorts_re  <- "nonintersection_coding|nonintersection|wholeset_coding|wholeset|stage1_coding|stage1|intersection_coding|intersection"
subsets_re  <- "123vs0|12vs0|3vs0"
fsmethod_re <- "mad_nosexmt|mad|cor"

parse_name <- function(path) {
  bn <- sub("_cvSummary\\.csv$", "", basename(path))
  varcut <- sub(".*_v", "", bn)
  rest   <- sub("_v[0-9]+$", "", bn)
  fsmethod <- regmatches(rest, regexpr(paste0("(", fsmethod_re, ")$"), rest))
  cat("  fsmethod=", fsmethod, "\n")
  rest     <- sub(paste0("_", fsmethod, "$"), "", rest)
  model <- sub(".*_", "", rest)
  rest  <- sub("_glm$", "", rest)
  subset <- regmatches(rest, regexpr(paste0("(", subsets_re, ")$"), rest))
  cat("  subset=", subset, "\n")
  rest   <- sub(paste0("_", subset, "$"), "", rest)
  set <- sub("_.*", "", rest)
  r2  <- sub("^[^_]+_", "", rest)
  cat("  r2=", r2, "\n")
  cohort <- regmatches(r2, regexpr(paste0("^(", cohorts_re, ")"), r2))
  cat("  cohort=[", cohort, "] len=", length(cohort), "\n")
  r3     <- sub(paste0("^", cohort, "_"), "", r2)
  cat("  r3=", r3, "\n")
  datatype   <- sub("_.*", "", r3)
  featureset <- sub("^[^_]+_", "", r3)
  cat("  datatype=", datatype, " featureset=", featureset, "\n")
  data.frame(set, cohort, datatype, featureset, subset, model, fsmethod, varcut,
             stringsAsFactors = FALSE)
}

tests <- c(
  "res/current_nonintersection_coding_fpkm_clinics_123vs0_glm_cor_v1000_cvSummary.csv",
  "res/prev_wholeset_normalized_fpkm_clinics_12vs0_glm_mad_nosexmt_v200_cvSummary.csv"
)
for (t in tests) {
  cat("\nFILE:", basename(t), "\n")
  m <- tryCatch(parse_name(t), error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(m)) {
    bad <- any(sapply(m, function(x) length(x) == 0 || x == ""))
    cat("  -> PARSED, would_be_skipped=", bad, "\n"); print(m)
  }
}
# what do the real nonintersection files actually look like?
f <- list.files("res", pattern="nonintersection.*_cvSummary\\.csv$")
cat("\n=== distinct nonintersection name shapes (first 6) ===\n")
print(head(f, 6))
