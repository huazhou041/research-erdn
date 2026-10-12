#!/usr/bin/env python3
# Python/sklearn MLP replacement for the crashing h2o deeplearning arm.
# No JVM -> no port-54321 crashes. Mirrors the RF/GLM CV discipline:
#   * 10-fold stratified CV x nboot repeats
#   * fold-local MAD feature selection (top-N genes; unsupervised, no leakage)
#   * per-fold StandardScaler fit on train only
#   * modest, regularized MLP (small hidden layer + strong L2 + early stopping)
#     because n~90 with p~20k is the overfit-prone regime
#   * POOLED out-of-fold AUC per repeat (honest estimator)
#   * matched PERMUTATION NULL in the same run (shuffle y, identical pipeline)
import sys, numpy as np, pandas as pd
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
# Input: a CSV whose first column is 'recurrence' (0/1) and the rest are genes.

def load_csv(path):
    return pd.read_csv(path)

def pooled_oof_auc(X, y, ntop=200, k=10, seed=0):
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    prob = np.full(len(y), np.nan)
    for tr, te in skf.split(X, y):
        mad = np.median(np.abs(X[tr] - np.median(X[tr],axis=0)), axis=0)
        sel = np.argsort(mad)[::-1][:ntop]
        sc = StandardScaler().fit(X[tr][:,sel])
        Xtr, Xte = sc.transform(X[tr][:,sel]), sc.transform(X[te][:,sel])
        clf = MLPClassifier(hidden_layer_sizes=(32,), alpha=1.0, activation='relu',
                            max_iter=300, early_stopping=True, random_state=seed)
        clf.fit(Xtr, y[tr])
        prob[te] = clf.predict_proba(Xte)[:, list(clf.classes_).index(1)]
    return roc_auc_score(y, prob)

if __name__ == "__main__":
    csv, nperm = sys.argv[1], (int(sys.argv[2]) if len(sys.argv)>2 else 20)
    df = load_csv(csv)
    y = df['recurrence'].astype(float).astype(int).values
    X = df.drop(columns=['recurrence']).astype(float).values
    keep = X.std(axis=0) > 0; X = X[:, keep]
    # NOTE: no global MAD pre-filter here -- feature selection is done strictly
    # fold-locally inside pooled_oof_auc (top-N MAD on the training fold only), so
    # the held-out fold never influences which features are kept. (A former global
    # top-3000 pre-filter let validation-fold expression distributions leak into
    # selection; removed for a strictly fold-local estimate.)
    obs = np.mean([pooled_oof_auc(X, y, seed=s) for s in range(3)])
    rng = np.random.RandomState(1)
    null = np.array([pooled_oof_auc(X, rng.permutation(y), seed=s) for s in range(nperm)])
    p = (1+np.sum(null>=obs))/(nperm+1)
    print(f"{csv.split('/')[-1]}: n={len(y)} events={int(y.sum())}  OBSERVED={obs:.3f}  "
          f"NULLmean={null.mean():.3f} 95pct={np.quantile(null,0.95):.3f} "
          f"max={null.max():.3f}  p={p:.3f}", flush=True)
