setwd("/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609")
f <- "res/current_wholeset_normalized_fpkm_123vs0_glm_natcommad_v200.rds"
x <- readRDS(f)
cat("class:", class(x), "\n")
cat("is list:", is.list(x), " length:", length(x), "\n")
if (is.list(x)) { cat("names:\n"); print(names(x)); for (nm in names(x)) cat(sprintf("  %-20s class=%s len/dim=%s\n", nm, class(x[[nm]])[1], paste(if(!is.null(dim(x[[nm]]))) dim(x[[nm]]) else length(x[[nm]]), collapse="x"))) }
# look for any gene-name vectors
find_genes <- function(o, path="") {
  if (is.character(o) && length(o)>1 && any(grepl("^[A-Z0-9]", o))) cat("  char vec at", path, "len", length(o), "ex:", paste(head(o,4),collapse=","), "\n")
  if (is.list(o)) for (nm in names(o)) find_genes(o[[nm]], paste0(path,"$",nm))
}
find_genes(x, "x")
