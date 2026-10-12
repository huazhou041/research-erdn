# For each dataset (set,cohort,datatype,dataname,subset): n, pos, nfeat,
# max |cor(feature, recurrence)|, which feature, and # features passing BH q<0.05.
sets    <- c("prev","current")
cohorts <- c("wholeset","wholeset_coding","stage1","stage1_coding","intersection","intersection_coding")
dtypes  <- c("normalized","fpkm")
subs    <- c("12vs0","3vs0","123vs0")
out <- list()
stat <- function(X, y) {
  X <- as.matrix(X)
  keep <- apply(X, 2, function(g) { g <- as.numeric(g); is.finite(sd(g, na.rm=TRUE)) && sd(g, na.rm=TRUE) > 0 })
  X <- X[, keep, drop=FALSE]
  if (!ncol(X)) return(list(mx=NA_real_, ft=NA_character_, nsig=NA_integer_, nf=0L))
  r <- suppressWarnings(as.vector(cor(X, y, use="pairwise.complete.obs")))
  r[!is.finite(r)] <- 0
  a <- abs(r); i <- which.max(a); n <- length(y)
  tt <- a*sqrt((n-2)/pmax(1-pmin(a,0.999999)^2, 1e-12))
  p  <- 2*pt(-abs(tt), df=n-2)
  q  <- p.adjust(p, method="BH")
  list(mx=a[i], ft=colnames(X)[i], nsig=sum(q<0.05), nf=ncol(X))
}
for (s in sets) for (co in cohorts) for (dt in dtypes) for (su in subs) {
  d <- file.path("data/model_data", s, co, dt)
  f_g <- file.path(d, paste0("fpkm_",su,".RDS"))
  f_c <- file.path(d, paste0("clinics_",su,".RDS"))
  f_m <- file.path(d, paste0("module20_",su,".RDS"))
  if (!all(file.exists(f_g,f_c,f_m))) next
  G <- readRDS(f_g); C <- readRDS(f_c); M <- readRDS(f_m)
  y <- as.numeric(as.character(G$recurrence))
  gen <- G[, setdiff(colnames(G),"recurrence"), drop=FALSE]
  cl  <- C[, setdiff(colnames(C),"recurrence"), drop=FALSE]
  mo  <- M[, setdiff(colnames(M),"recurrence"), drop=FALSE]
  cl <- cl[rownames(gen), , drop=FALSE]; mo <- mo[rownames(gen), , drop=FALSE]
  fsets <- list(
    fpkm             = gen,
    fpkm_clinics     = cbind(cl, gen),
    clinics          = cl,
    module20         = mo,
    module20_clinics = cbind(cl, mo))
  for (nm in names(fsets)) {
    st <- stat(fsets[[nm]], y)
    out[[length(out)+1]] <- data.frame(set=s, cohort=co, datatype=dt, dataname=nm, subset=su,
      n=nrow(gen), pos=sum(y==1), neg=sum(y==0), nfeat=st$nf,
      max_abs_cor=round(st$mx,4), top_feature=st$ft, n_sig_BH05=st$nsig, stringsAsFactors=FALSE)
  }
  cat("done:", s, co, dt, su, "\n")
}
write.csv(do.call(rbind,out), "maxcor_by_dataset.csv", row.names=FALSE)
