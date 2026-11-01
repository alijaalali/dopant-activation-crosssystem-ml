import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor

from common import DATA, FIG
from models import arrhenius_log_time, inv_logit, logit_act, reml_profile

warnings.filterwarnings("ignore")
db = pd.read_csv(DATA / "database.csv")
al = db[(db.system == "SiC") & (db.dopant == "Al")].copy()
al["log_time"] = np.log10(al.anneal_time_min)
al["log_conc"] = np.log10(al.implanted_conc_cm3)
al = al.dropna(subset=["activation_pct", "anneal_temp_C", "anneal_time_min", "implant_temp_C", "implanted_conc_cm3"])
al = al[(al.definition_class != "approximate") & (~al.two_step_anneal.fillna(False).astype(bool))].reset_index(drop=True)
al["y"] = logit_act(al.activation_pct)
EA_GRID = [0.5, 1, 2, 3, 4, 5, 6, 8, 10]
OUT = {"n_rows": len(al), "n_sources": int(al.source_id.nunique())}


def design(frame, kind, ea=None):
    if kind == "kinetic":
        te = arrhenius_log_time(frame.anneal_temp_C, frame.anneal_time_min, ea)
        return np.column_stack([frame.implant_temp_C, frame.log_conc, te])
    if kind == "linear":
        return np.column_stack([frame.implant_temp_C, frame.anneal_temp_C, frame.log_time, frame.log_conc])
    if kind == "linear_no_implant":
        return np.column_stack([frame.anneal_temp_C, frame.log_time, frame.log_conc])
    raise ValueError(kind)


def mixed_predict(train, test, kind, ea=None):
    Xtr, Xte = design(train, kind, ea), design(test, kind, ea)
    mu, sd = Xtr.mean(0), Xtr.std(0, ddof=1)
    sd[sd == 0] = 1
    A = np.column_stack([np.ones(len(train)), (Xtr - mu) / sd])
    g = pd.factorize(train.source_id)[0]
    if g.max() == 0:
        beta = np.linalg.lstsq(A, train.y.to_numpy(), rcond=None)[0]
    else:
        beta = reml_profile(A, train.y.to_numpy(), g)[0]
    return 100 * inv_logit(np.column_stack([np.ones(len(test)), (Xte - mu) / sd]) @ beta)


def kinetic_predict(train, test):
    best, best_err = None, np.inf
    for ea in EA_GRID:
        errs = []
        for s in train.source_id.unique():
            tr, te = train[train.source_id != s], train[train.source_id == s]
            if tr.source_id.nunique() < 2:
                continue
            errs.append(np.abs(mixed_predict(tr, te, "kinetic", ea) - te.activation_pct.to_numpy()))
        if errs:
            e = np.concatenate(errs).mean()
            if e < best_err:
                best, best_err = ea, e
    ea = best if best is not None else 2.0
    return mixed_predict(train, test, "kinetic", ea)


def mean_predict(train, test):
    return np.full(len(test), 100 * inv_logit(train.y.mean()))


def rf_predict(train, test):
    m = RandomForestRegressor(n_estimators=200, max_depth=4, random_state=0, min_samples_leaf=2)
    m.fit(design(train, "linear"), train.activation_pct)
    return m.predict(design(test, "linear"))


MODELS = {
    "training_mean": mean_predict,
    "linear_mixed": lambda a, b: mixed_predict(a, b, "linear"),
    "linear_mixed_without_implantT": lambda a, b: mixed_predict(a, b, "linear_no_implant"),
    "arrhenius_kinetic_mixed": kinetic_predict,
    "random_forest": rf_predict,
}


def evaluate(splits):
    res = {}
    for name, fn in MODELS.items():
        errs, preds, trues = [], [], []
        for tr, te in splits:
            p = fn(tr, te)
            errs.append(np.abs(p - te.activation_pct.to_numpy()))
            preds.append(p)
            trues.append(te.activation_pct.to_numpy())
        e = np.concatenate(errs)
        P, T = np.concatenate(preds), np.concatenate(trues)
        res[name] = {"mae_pp": float(e.mean()), "median_abs_error_pp": float(np.median(e)), "n_test": int(len(e)),
                     "spearman": float(stats.spearmanr(P, T)[0]) if np.std(P) > 0 else None}
    return res


def leave_one(col):
    return [(al[al[col] != v], al[al[col] == v]) for v in al[col].unique()]


OUT["leave_one_source_out"] = evaluate(leave_one("source_id"))
al["lab"] = al.lab_group.fillna(al.source_id)
OUT["leave_one_lab_out"] = evaluate(leave_one("lab"))
old, new = al[al.year <= 2005], al[al.year >= 2017]
OUT["temporal"] = {"train_years": "<=2005", "test_years": ">=2017", "n_train": len(old), "n_test": len(new),
                   "train_sources": sorted(old.source_id.unique()), "test_sources": sorted(new.source_id.unique()),
                   **evaluate([(old, new)])}
for k in ["leave_one_source_out", "leave_one_lab_out", "temporal"]:
    print(k)
    for m, v in OUT[k].items():
        if isinstance(v, dict):
            print("  ", m, round(v["mae_pp"], 1), v["spearman"])

prof = {}
for ea in EA_GRID:
    errs = []
    for s in al.source_id.unique():
        tr, te = al[al.source_id != s], al[al.source_id == s]
        errs.append(np.abs(mixed_predict(tr, te, "kinetic", ea) - te.activation_pct.to_numpy()))
    prof[ea] = float(np.concatenate(errs).mean())
base = {}
for kind in ["linear", "linear_no_implant"]:
    errs = []
    for s in al.source_id.unique():
        tr, te = al[al.source_id != s], al[al.source_id == s]
        errs.append(np.abs(mixed_predict(tr, te, kind) - te.activation_pct.to_numpy()))
    base[kind] = float(np.concatenate(errs).mean())
OUT["kinetic_Ea_profile_loso_mae_pp"] = prof
OUT["linear_loso_mae_pp"] = base
print("Ea profile", prof, base)
json.dump(OUT, open(DATA / "kinetic_validation_results.json", "w"), indent=2)

fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.plot(list(prof), list(prof.values()), "o-", label="Arrhenius-equivalent time")
ax.axhline(base["linear"], color="tab:red", linestyle="--", label="linear (anneal T, log time)")
ax.set_xlabel("Assumed activation energy $E_a$ (eV)")
ax.set_ylabel("Leave-one-source-out MAE (pp)")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIG / "kinetic_profile.png", dpi=150)
plt.savefig(FIG / "kinetic_profile.pdf")
print("DONE")
