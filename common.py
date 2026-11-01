import itertools
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def load():
    df = pd.read_csv(DATA / "activation_dataset.csv", na_values=["NA"])
    for raw, new in [("anneal_temp_C", "anneal_T"), ("implant_temp_C", "implant_T"),
                     ("anneal_time_min", "time"), ("activation_pct", "act")]:
        df[new] = pd.to_numeric(df[raw], errors="coerce")
    df.loc[(df.source_id == "S1_S2") & df.implant_T.isna(), "implant_T"] = 500.0
    return df


def sic_model_frame(df):
    sic = df[df.system == "SiC"]
    return sic.dropna(subset=["anneal_T", "implant_T", "act"]).reset_index(drop=True)


def zscore(X):
    X = np.asarray(X, float)
    return (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)


def design(X):
    return np.column_stack([np.ones(len(X)), X])


def ols(X, y):
    A = design(X)
    return np.linalg.lstsq(A, y, rcond=None)[0]


def r_squared(X, y):
    A = design(X)
    b = np.linalg.lstsq(A, y, rcond=None)[0]
    ss_res = np.sum((y - A @ b) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    n, p = A.shape
    return r2, 1 - (1 - r2) * (n - 1) / (n - p)


def cluster_enumeration(sources, X, y):
    """Exact source-level bootstrap: every ordered resample of the sources.

    Each ordered tuple of sources has equal probability, so percentiles over
    the enumerated fits are the exact bootstrap percentiles. Resamples with a
    single distinct source or a rank-deficient design are excluded and counted.
    """
    sources = np.asarray(sources)
    uniq = sorted(set(sources))
    rows = {s: np.where(sources == s)[0] for s in uniq}
    betas, n_single, n_rank = [], 0, 0
    for combo in itertools.product(uniq, repeat=len(uniq)):
        if len(set(combo)) < 2:
            n_single += 1
            continue
        idx = np.concatenate([rows[s] for s in combo])
        A = design(X[idx])
        if np.linalg.matrix_rank(A) < A.shape[1]:
            n_rank += 1
            continue
        betas.append(np.linalg.lstsq(A, y[idx], rcond=None)[0])
    return {"betas": np.array(betas), "n_total": len(uniq) ** len(uniq),
            "n_single_source": n_single, "n_rank_deficient": n_rank}


def summarize_enumeration(res, names):
    b = res["betas"]
    out = {"n_total_resamples": res["n_total"], "n_single_source": res["n_single_source"],
           "n_rank_deficient": res["n_rank_deficient"], "n_valid": int(len(b)), "terms": {}}
    for j, name in enumerate(names, start=1):
        v = b[:, j]
        lo, hi = np.percentile(v, [2.5, 97.5])
        out["terms"][name] = {"ci95": [float(lo), float(hi)], "min": float(v.min()),
                              "max": float(v.max()), "n_distinct": int(len(np.unique(np.round(v, 8)))),
                              "share_positive": float(np.mean(v > 0))}
    return out


def reml_random_intercept(X, y, groups):
    from scipy.optimize import minimize_scalar
    A = design(X)
    n, p = A.shape
    groups = np.asarray(groups)
    idx = [np.where(groups == g)[0] for g in sorted(set(groups))]

    def fit(logr):
        V = np.eye(n)
        for ii in idx:
            V[np.ix_(ii, ii)] += np.exp(logr)
        Vi = np.linalg.inv(V)
        M = A.T @ Vi @ A
        beta = np.linalg.solve(M, A.T @ Vi @ y)
        res = y - A @ beta
        s2 = float(res @ Vi @ res) / (n - p)
        ll = -0.5 * ((n - p) * np.log(s2) + np.linalg.slogdet(V)[1] + np.linalg.slogdet(M)[1])
        return beta, s2, M, ll

    opt = minimize_scalar(lambda t: -fit(t)[3], bounds=(-8, 6), method="bounded")
    beta, s2, M, _ = fit(opt.x)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(M)))
    return {"var_ratio_source_to_residual": float(np.exp(opt.x)), "beta": beta, "se": se,
            "n_groups": len(idx)}
