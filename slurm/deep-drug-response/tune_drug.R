#!/usr/bin/env Rscript
# inputs    -------------------------------------------------------------------
library(argparser)

p = arg_parser("Run training pipeline")
p = add_argument(p, "drug", help = "Drug to learn")
p = add_argument(p, "matrix", help = "Matrix to train on (RDS file)")
p = add_argument(p, "hyperparams", help = "json file with hyper parameters")
p = add_argument(p, "prefix", help = "prefix to write files")
p = add_argument(p, "--max_models", default = 15L,
                 help = "Max number of models to compute")
p = add_argument(p, "--varcut", default = 100,
                 help = "Percent of genes to keep based on variance.")
p = add_argument(p, "--nboot", default = 50,
                 help = "Number of bootstrap.")
p = add_argument(p, "--learner", default = "DL",
                 help = "Algorithm to use for learning (DL, RF, GLM)")
p = add_argument(p, "--test", help = "Path to clinical data")
p = add_argument(p, "--kfold", default = 5L, help = "Number of CV-splits")
p = add_argument(p, "--max_mem", default = "64G", help = "H2O parameter")
p = add_argument(p, "--cores", default = "Slurm",
                 help = "Method to determine cores with future::availableCores")
p = add_argument(p, "--fs_method", default = "mad",
                 help = "Gene prefilter method: mad (unsupervised variance, original), cor (supervised: abs correlation with outcome, computed fold-locally), or mad_nosexmt (mad, but excluding sex-chromosome and mitochondrial genes first)")
p = add_argument(p, "--exclude_genes", default = NA,
                 help = "Optional: path to a newline-delimited gene list to exclude before ranking (legacy; mad_nosexmt now derives the set inline)")
p = add_argument(p, "--exclude_sexmt_annot", default = NA,
                 help = "Path to a gene annotation CSV (gene_name,chr,...) used by fs_method=mad_nosexmt to derive sex-chromosome (chrX/chrY) genes INLINE; mitochondrial genes are embedded in this script (annotation carries no chrM). No precomputed exclude file needed.")
p = add_argument(p, "--rank_genes_file", default = NA,
                 help = "For fs_method=natcommad: CSV with a 'gene' column pre-ranked by nature_com MAD (descending). The top -v genes present in the data are used as a FIXED external feature set for every fold (transfer of the nature_com MAD gene list to our data; no fold-local ranking, no outcome use).")

argv = parse_args(p)
# argv = parse_args(p, c("Bortezomib", "bortezomib_train.rds", "hyperparam_grid_dl.json",
#                        "debug", "--max_models", "2", "--varcut", "50"))

message("Running with the following arguments")
# str(argv)

drug = argv$drug
mat_file = argv$matrix
lift_null_file = argv$lift
hyperparam_file = argv$hyperparams
prefix = argv$prefix
Nmodels = argv$max_models
varcut = argv$varcut
learner = argv$learner
kfold = argv$kfold
kboot=argv$nboot
max_mem = argv$max_mem
test_file = argv$test
cores_method = argv$cores
fs_method = argv$fs_method
exclude_genes_file = argv$exclude_genes
exclude_sexmt_annot = argv$exclude_sexmt_annot
rank_genes_file = argv$rank_genes_file

stopifnot(learner %in% c("DL", "RF", "GLM"))
# GUARD: the DL/RF branches select the best inner-CV model via glmnet-only fields
# (mod$cvm / mod$lambda.min, ~lines 348/428), which do not exist for h2o/randomForest
# objects -- so those paths are broken. Use the dedicated scripts instead.
if (learner %in% c("DL", "RF"))
  stop("tune_drug.R --learner ", learner,
       " is unsupported (glmnet-only model selection). Use tune_drug_rf.R, ",
       "tune_drug_dl.R, or tune_drug_xgb.R instead.")
stopifnot(fs_method %in% c("mad", "cor", "mad_nosexmt", "natcommad"))

# fs_method=natcommad: load the nature_com MAD-ranked gene list ONCE (descending
# MAD over nature_com normal-tissue samples; built by build_natcom_madgenes.R).
# This is a FIXED external feature set applied to our data -- the top -v of these
# genes that are present are used for every CV fold, so it neither ranks on our
# data nor sees the outcome (a clean transfer of nature_com's gene list).
natcom_ranked = character(0)
if (fs_method == "natcommad") {
    stopifnot(!is.na(rank_genes_file), file.exists(rank_genes_file))
    .nr = read.csv(rank_genes_file, check.names = FALSE, stringsAsFactors = FALSE)
    gcol_nr = if ("gene" %in% names(.nr)) "gene" else names(.nr)[1]
    natcom_ranked = as.character(.nr[[gcol_nr]])            # already in MAD-descending order
    message(sprintf("natcommad: loaded %d nature_com MAD-ranked genes from %s", length(natcom_ranked), rank_genes_file))
}

# Mitochondrial genes to exclude for fs_method=mad_nosexmt, embedded INLINE because
# the gene annotation carries no chrM (only chrX/chrY). Covers the standard 37 MT
# genes in both the "MT-" and bare-alias naming conventions plus the MTRNR2L
# nuclear-mitochondrial pseudogene family (86 names total). Sex-chromosome genes are
# derived from the annotation's chr column at run time; union of the two reproduces
# the former exclude_sex_mt_genes.txt exactly.
MT_GENES = c(
  "ATP6","ATP8","COX1","COX2","COX3","CYTB",
  "MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB",
  "MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6","MT-RNR1","MT-RNR2",
  "MT-TA","MT-TC","MT-TD","MT-TE","MT-TF","MT-TG","MT-TH","MT-TI","MT-TK","MT-TL1","MT-TL2",
  "MT-TM","MT-TN","MT-TP","MT-TQ","MT-TR","MT-TS1","MT-TS2","MT-TT","MT-TV","MT-TW","MT-TY",
  "MTRNR2L1","MTRNR2L11","MTRNR2L12","MTRNR2L13","MTRNR2L2","MTRNR2L3","MTRNR2L4",
  "MTRNR2L5","MTRNR2L6","MTRNR2L7","MTRNR2L8","MTRNR2L9",
  "ND1","ND2","ND3","ND4","ND4L","ND5","ND6","RNR1","RNR2",
  "TA","TC","TD","TE","TF","TG","TH","TI","TK","TL1","TL2","TM","TN","TP","TQ","TR","TS1","TS2","TT","TV","TW","TY")

exclude_gene_set = character(0)
if (fs_method == "mad_nosexmt") {
    if (!is.na(exclude_sexmt_annot)) {
        stopifnot(file.exists(exclude_sexmt_annot))
        ann = read.csv(exclude_sexmt_annot, check.names = FALSE, stringsAsFactors = FALSE)
        chrcol = if ("chr" %in% names(ann)) "chr" else grep("chr|chrom|seqname", names(ann), ignore.case = TRUE, value = TRUE)[1]
        gcol   = if ("gene_name" %in% names(ann)) "gene_name" else grep("gene_name|^gene$|symbol", names(ann), ignore.case = TRUE, value = TRUE)[1]
        stopifnot(!is.na(chrcol), !is.na(gcol))
        sex_genes = unique(ann[[gcol]][ann[[chrcol]] %in% c("chrX", "chrY", "X", "Y")])
        exclude_gene_set = union(sex_genes, MT_GENES)
        message(sprintf("mad_nosexmt: excluding %d genes INLINE (%d sex-chromosome from annotation chr column + %d mitochondrial embedded) before MAD ranking",
                        length(exclude_gene_set), length(sex_genes), length(MT_GENES)))
    } else if (!is.na(exclude_genes_file)) {
        stopifnot(file.exists(exclude_genes_file))                 # legacy fallback
        exclude_gene_set = readLines(exclude_genes_file)
        message(sprintf("mad_nosexmt: excluding %d genes from legacy file %s", length(exclude_gene_set), exclude_genes_file))
    } else {
        stop("fs_method=mad_nosexmt requires --exclude_sexmt_annot <gene_annotation.csv> (preferred) or --exclude_genes <file>")
    }
}

Ncores = future::availableCores(method = cores_method)
options(mc.cores = Ncores)
cat('*****************\n')
cat(learner,' : ' ,mat_file,'\n')
cat('*****************\n')

# Library    -------------------------------------------------------------------

suppressPackageStartupMessages({
    library(tidyverse)
    library(parallel)
    library(glmnet)
    library(glmnetUtils)
    library(h2o)
	library(pROC)
  library(caret)
})

all_equal = function(x, y) {
    if (length(y) == 1)
        y = rep(y, length(x))
    !is.character(all.equal(x, y, check.attributes = FALSE))
}

is_scaled = function(mat) {
    # result of `stat::scale`
    if (!is.null(attr(mat, "scaled:center")) )
        return(TRUE)
    # zero-mean
    cmeans = colMeans(mat, na.rm = TRUE)
    if (!all_equal(cmeans, 0))
        return(FALSE)
    # unit variance
    csd = apply(mat, 2, sd)
    if (!all_equal(csd, 1))
        return(FALSE)

    return(TRUE)
}

rbind_data_frame = function(...)
    # simple wrapper
    rbind.data.frame(..., make.row.names = FALSE, stringsAsFactors = FALSE)

delist_col = function(x, collapse = '|') {
    if (!is.list(x))
        return(x)
    vapply(x, paste, '', sep = '', collapse = collapse)
}

delist_df = function(x, collapse = '|')
    as.data.frame(lapply(x, delist_col, collapse = collapse))

FS_top = function(mat, top) {
    cand_genes = genes
    if (fs_method == "natcommad") {
        # FIXED external feature set: the top-`top` nature_com MAD-ranked genes that
        # are present in this dataset, in nature_com's MAD order. Identical for every
        # fold (no ranking on `mat`, no outcome use) -- a transfer of nature_com's gene
        # list. Clinical `tissues` columns are kept as in the other methods.
        ranked_present = natcom_ranked[natcom_ranked %in% genes]
        sel = head(ranked_present, top)
        if (length(sel) > 1) return(c(tissues, sel)) else return(tissues)
    }
    if (fs_method == "mad_nosexmt") {
        # Same unsupervised MAD ranking as the original, but with sex-chromosome and
        # mitochondrial genes removed from the candidate pool first -- these were found
        # to dominate the top-MAD list purely due to patient sex / technical artifacts,
        # unrelated to the recurrence outcome (see Section 7 diagnosis).
        cand_genes = setdiff(genes, exclude_gene_set)
    }
    if (fs_method == "cor") {
        # Supervised: rank genes by |correlation| with the outcome. `mat` here is
        # whatever the caller already restricted (a single CV fold's training rows,
        # or the whole dataset for the final model) -- same fold-local discipline as
        # the mad() branch below, just a different (outcome-aware) ranking statistic.
        y = suppressWarnings(as.numeric(as.character(mat[, drug])))
        ix = apply(mat[, cand_genes], 2, function(g)
            abs(suppressWarnings(cor(as.numeric(g), y, use = "pairwise.complete.obs"))))
        ix[is.na(ix)] = 0
    } else {
        ix = apply(mat[, cand_genes], 2, mad, na.rm = TRUE)
    }
    if (length(ix)<top)
	top=length(ix)
    ix = order(ix, decreasing = TRUE)[1:top]
    genes = cand_genes[ix]
    #if (!is.na(genes))
if (length(genes)>1)
    c(tissues,genes)
	else
	tissues
}

expand_hyperparams_h2o = function(hyper_params) {
    hyper_params = split(hyper_params, seq(nrow(hyper_params)))
    hyper_params = lapply(hyper_params, as.list)
    names(hyper_params) = NULL
    for (i in seq(length(hyper_params))) {
        for (j in seq(length(hyper_params[[i]]))) {
            if (is.list(hyper_params[[i]][[j]]))
                hyper_params[[i]][[j]] = unlist(hyper_params[[i]][[j]],
                                                recursive = FALSE,
                                                use.names = FALSE)
        }
    }
    hyper_params
}

validate_h2o = function(model, newdata) {
    y_hat = h2o.performance(model, newdata)
    data.frame(MSE = h2o.mse(y_hat), AUC = h2o.auc(y_hat))
}

safe_auc = function(y, prob) {
  y = as.numeric(y)
  if (length(unique(y[is.finite(y)])) < 2) {
    return(NA_real_)
  }
  # direction="<": predictor (prob) assumed LOWER for recurrence=0 than recurrence=1.
  # Without this, pROC's default direction="auto" silently flips per-call based on
  # whichever group has the higher median in that fold -- inflating AUC when the
  # model is actually anti-correlated with truth in that fold.
  as.numeric(pROC::auc(y, as.numeric(prob), quiet = TRUE, direction = "<"))
}

validate_glmnet = function(model, newdata, s = 'lambda.min') {
  idrug = which(colnames(newdata) == drug)
  #//y_hat = predict(model, as(newdata[, -idrug],'dgCMatrix'), s = s)
  y_hat = predict(model, newdata[, -idrug], s = s)
  dy = y_hat - newdata[, idrug]
  data.frame(lambda = model[[s]],  # hacky
             MSE = mean(dy * dy), AUC = safe_auc(newdata[, idrug], y_hat))
}

# Preprocessing -----------------------------------------------------------------

message(paste0('Sourcing ', hyperparam_file, "..."))
hyper_params <- jsonlite::fromJSON(hyperparam_file,
                                   simplifyVector = TRUE,
                                   simplifyDataFrame = FALSE,
                                   simplifyMatrix = FALSE)

message(paste0('Reading ', mat_file, '...'))
#mat = as.matrix(readRDS(mat_file))
mat = readRDS(mat_file)
mat$recurrence=as.factor(mat$recurrence) ##hua zhou 2020.03.13
tissues = grep("^Tissue_", colnames(mat), value = TRUE)
if (any(grepl("_GE", colnames(mat)))) {
    genes = grep("_GE$", colnames(mat), value = TRUE)
} else {
    genes = setdiff(colnames(mat), c(drug, tissues))
}

message('Subseting matrix...')
mat = mat[!is.na(mat[, drug]),  ]
origmat=mat

set.seed(123)

cvtotal=vector("list", kboot)

cvtotal2 = vector("list", kboot)

###doing kboot boostraps
for (tt in 1:kboot) {
# cat('doing ',tt,'\n') 
  # k-fold cross-validation
  mat=origmat
#mat=mat[sample(1:nrow(mat),nrow(mat),replace=F),]
fold_id = cut(1:nrow(mat), breaks = kfold, labels = FALSE)

flds <- createFolds(mat[,1], k = kfold, list = TRUE, returnTrain = FALSE)
for (i in 1:kfold) {
  fold_id[flds[[i]]]=i
}
fold_size = sum(fold_id != 1L)

message(sprintf('Selecting genes based on top %f genes cutoff', varcut))
#top = floor(varcut * length(genes) / 100)
top = varcut
features = mclapply(1:kfold, function(i) FS_top(mat[fold_id != i, ], top))
features[[kfold + 1L]] = FS_top(mat, top)
if (length(features[[1]])<100)
cat(features[[1]])

message('Setting up learning algorithm...')
if (learner == "DL") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.deeplearning
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else if (learner == "RF") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.randomForest
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else {
    # GLM
    validate = validate_glmnet
    convert = function(x, ...) {
		x[,1]=(x[,1]==1)+0
			as.matrix(x)}
# learn = function(training_frame, y, x, alpha) {
#        cv.glmnet(training_frame[, x], training_frame[, y],
#                  parallel = TRUE, alpha = alpha, family = "gaussian")
#    }
  #  learn = function(training_frame, y, x, alpha) {
  #      cv.glmnet(training_frame[, x], as.vector(training_frame[, y]),
  #                parallel = TRUE, alpha = alpha, family = "binomial") }
###hua added, prevent failing in lambda default
learn=function(training_frame, y, x, alpha) {
  yv = as.vector(training_frame[, y])
  # cv.glmnet dies with "one binomial class has 1 or 0 observations" on rare-event data
  # (e.g. novel/3vs0: 4 events). Two causes, both handled here:
  #   1. nfolds=10 is more folds than the rare class has members.
  #   2. cv.glmnet assigns inner folds AT RANDOM -- it does not stratify -- so even with
  #      nfolds=3 and 3 events a fold can end up with zero. We therefore pass an explicit
  #      class-stratified foldid so every inner fold is guaranteed both classes.
  # nfolds=10 unstratified is still attempted FIRST, so datasets that already ran are
  # bit-for-bit unaffected; this only engages where the old code crashed.
  strat_foldid = function(yv, nf) {
    fid = integer(length(yv))
    for (cl in unique(yv)) {
      idx = which(yv == cl)
      fid[idx] = sample(rep_len(seq_len(nf), length(idx)))
    }
    fid
  }
  nmin = min(table(yv))
  nf = max(3L, min(10L, as.integer(nmin)))
  result=tryCatch( {
    cv.glmnet(training_frame[, x], yv,
              parallel = TRUE, alpha = alpha, family = "binomial") },
    error=function(err) {
      cat('something wrong')
      fid = strat_foldid(yv, nf)
      tryCatch({
        cv.glmnet(training_frame[, x], yv, parallel = TRUE, alpha = alpha,
                  foldid = fid, family = "binomial")
      }, error=function(err2) {
        cv.glmnet(training_frame[, x], yv, parallel = TRUE, alpha = alpha, foldid = fid,
                  lambda=exp(seq(log(0.001), log(5), length.out=100)), family = "binomial")
      })
    })
  return(result)
}



    hyper_params = hyper_params
    hyper_params_df = data.frame(model_id = 1:length(hyper_params),
                                 alpha = hyper_params)
    saveModel = function(x, nam) saveRDS(x, paste0(nam, ".rds"))
    Nmodels = length(hyper_params)
}

message('Running cross-validation')

cvSummary = vector("list", kfold)

cvSummary2=data.frame(matrix(NaN,0,4))
names(cvSummary2)=c('subj','foldid','bootid','prob')
for (i in seq_len(kfold)) {
    kmat = mat[, c(drug, features[[i]]), drop = FALSE]
    ix = which(fold_id == i)
    train_hex = convert(kmat[-ix, , drop = FALSE])
    valid_hex = convert(kmat[ix, , drop = FALSE])
 if (length(features[[i]])>1 || learner!='GLM') {
    kmodel = lapply(hyper_params, function(params) {
                    x = features[[i]]
                    if ('mtries' %in% names(params)) {# special case
                        params[['mtries']] = round(params[['mtries']] * length(x))
			  if (params[['mtries']]==0)
                          params[['mtries']]=-1  ##hua added to block mtries=0 error
                    }
                    do.call(learn, c(list(y = drug, x = x, training_frame = train_hex), params))
                  })
    cvms = sapply(kmodel, function(mod) mod$cvm[mod$lambda == mod$lambda.min])
    imod = which.min(cvms)
    kmodel = kmodel[imod]
    
    
    ####get pred prob
    
    
    pred_t=predict(kmodel[[1]],valid_hex[,-which(colnames(valid_hex) == drug)],type='response')
    
    pred_t2=data.frame(subj=rownames(pred_t),foldid=i,bootid=tt,prob=pred_t[,1])
    cvSummary2=rbind.data.frame(cvSummary2,pred_t2)
    
    
    
    kvalid = lapply(kmodel, validate, newdata = valid_hex)
    kvalid = do.call(rbind, kvalid)
    kvalid$fold_id = i
    # kvalid$model_id = imod
    tmp_df = kvalid
    # tmp_df = cbind(data.frame(fold_id = i),
    #                hyper_params_df[imod],
    #                kvalid)
    # rownames(tmp_df) = NULL
    cvSummary[[i]] = tmp_df } else {
    kmodel=glm(as.formula(paste0(drug,'~',features[[i]])),family='binomial',data=data.frame(train_hex))
    cauc=safe_auc(valid_hex[,1], predict(kmodel,newdata=data.frame(valid_hex)))

    kvalid = data.frame(model_id=1,MSE=0,AUC=cauc)
    tmp_df = cbind(data.frame(fold_id = rep(i, 1)),
                   hyper_params_df[1,],
                   kvalid)
       pred_t=predict(kmodel,newdata=data.frame(valid_hex),type='response')
    
    pred_t2=data.frame(subj=names(pred_t),foldid=i,bootid=tt,prob=pred_t)
    cvSummary2=rbind.data.frame(cvSummary2,pred_t2)
    rownames(tmp_df) = NULL
    cvSummary[[i]] = kvalid }

}



# Write results -----------------------------------------------------------
cvSummary = do.call(rbind, cvSummary)
cvSummary$bootid=tt
# cvSummary = delist_df(cvSummary)
# cvSummary = cvSummary[order(cvSummary[["AUC"]],decreasing = T), ]
cvtotal[[tt]]=cvSummary
cvtotal2[[tt]]=cvSummary2
#write_csv(cvSummary, paste0(prefix, "_cv.csv"))
#stop()
# 
}

cvtotal = do.call(rbind, cvtotal)
cvtotal2 = do.call(rbind, cvtotal2)
write_csv(cvtotal, paste0(prefix, "_cvSummary.csv"))

write_csv(cvtotal2, paste0(prefix, "_pred_prob.csv"))



# Train best model on all of the data
#best_model = cvSummary$model_id[nrow(cvSummary)]
train_hex = convert(mat)
params = hyper_params[[1]]
x = features[[kfold + 1L]]
if ('mtries' %in% names(params))
  params[['mtries']] = floor(params[['mtries']] * length(x))
if (length(x)>1 || learner != 'GLM') {
  kmodel = lapply(hyper_params, function(params) {
    x = features[[kfold + 1L]]  ##hua 2026: use whole-data feature set for final model
    if ('mtries' %in% names(params)) {# special case
      params[['mtries']] = round(params[['mtries']] * length(x))
      if (params[['mtries']]==0)
        params[['mtries']]=-1  ##hua added to block mtries=0 error
    }
    do.call(learn, c(list(y = drug, x = x, training_frame = train_hex), params))
  })
  cvms = sapply(kmodel, function(mod) mod$cvm[mod$lambda == mod$lambda.min])
  imod = which.min(cvms)
  model=kmodel[[imod]] } else
    model=glm(as.formula(paste0(drug,'~',x)),family='binomial',data=data.frame(train_hex))
saveModel(model, prefix)







if (!is.na(test_file)) {
    test = readRDS(test_file)
    patient_id = rownames(test)
    test = convert(test[, x, drop = FALSE])
    y = as.vector(predict(model, test))
    df = tibble(patient_id = patient_id, ic50 = y)
    write_csv(df, paste0(prefix, "_clinicalPrediction.csv"))
}

if (learner %in% c("DL", "RF"))
    h2o.shutdown(FALSE)


cat('done\n',date())

