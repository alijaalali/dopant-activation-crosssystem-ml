import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import DATA, FIG
from models import gibbs_mixed, inv_logit, logit_act, summarize

db = pd.read_csv(DATA / "database.csv")
al = db[(db.system == "SiC") & (db.dopant == "Al")].copy()
al["log_time"] = np.log10(al.anneal_time_min)
al["log_conc"] = np.log10(al.implanted_conc_cm3)
need = ["activation_pct", "anneal_temp_C", "anneal_time_min", "implant_temp_C", "implanted_conc_cm3"]
al = al.dropna(subset=need)
al = al[(al.definition_class != "approximate") & (~al.two_step_anneal.fillna(False).astype(bool))].reset_index(drop=True)
FEATS = ["implant_temp_C", "anneal_temp_C", "log_time", "log_conc"]
print("rows", len(al), "sources", al.source_id.nunique(), al.source_id.value_counts().to_dict())


def fit(frame, feats=FEATS, use_lab=True, beta_sd=5.0, seed=1, n_iter=8000, burn=3000):
    X = frame[feats].to_numpy(float)
    mu, sd = X.mean(0), X.std(0, ddof=1)
    Xz = (X - mu) / sd
    A = np.column_stack([np.ones(len(frame)), Xz])
    y = logit_act(frame.activation_pct)
    src = pd.factorize(frame.source_id)[0]
    groups = [src]
    if use_lab:
        groups.append(pd.factorize(frame.lab_group.fillna(frame.source_id))[0])
    dr = gibbs_mixed(A, y, groups, n_iter=n_iter, burn=burn, seed=seed, beta_sd=beta_sd)
    s = summarize(dr, ["intercept"] + feats)
    s["sd_source"] = float(np.sqrt(dr["s2u"][:, 0]).mean())
    s["sd_residual"] = float(np.sqrt(dr["s2e"]).mean())
    mid = np.zeros(len(feats))
    lo, hi = mid.copy(), mid.copy()
    pp = {}
    for j, f in enumerate(feats):
        a, b = np.zeros(A.shape[1]), np.zeros(A.shape[1])
        a[0] = b[0] = 1
        a[1 + j], b[1 + j] = 0.5, -0.5
        eff = 100 * (inv_logit(dr["beta"] @ a) - inv_logit(dr["beta"] @ b))
        pp[f] = {"mean": float(eff.mean()), "ci95": [float(np.percentile(eff, 2.5)), float(np.percentile(eff, 97.5))],
                 "p_positive": float(np.mean(eff > 0))}
    s["marginal_effect_pp_per_sd"] = pp
    s["feature_sd"] = {f: float(v) for f, v in zip(feats, sd)}
    return s, dr, (mu, sd)


OUT = {"n_rows": len(al), "n_sources": int(al.source_id.nunique()), "rows_per_source": al.source_id.value_counts().to_dict()}
main, draws, (mu, sd) = fit(al)
OUT["main"] = main
print("main", json.dumps({k: main[k] for k in ["implant_temp_C", "anneal_temp_C", "log_time", "log_conc"]}))
sens = {}
for label, frame, kw in [
        ("drop_Saks2004", al[al.source_id != "Saks2004"], {}),
        ("drop_Saks2001", al[al.source_id != "Saks2001_JAP"], {}),
        ("drop_ECS92_implantT_is_range", al[al.source_id != "Nipoti2019_ECSTrans"], {}),
        ("drop_Pernot", al[al.source_id != "Pernot2005_JAP"], {}),
        ("no_lab_effect", al, {"use_lab": False}),
        ("prior_sd_2", al, {"beta_sd": 2.0}),
        ("prior_sd_10", al, {"beta_sd": 10.0}),
        ("without_conc", al, {"feats": ["implant_temp_C", "anneal_temp_C", "log_time"]})]:
    s, _, _ = fit(frame, **kw)
    sens[label] = {"n_rows": len(frame), "n_sources": int(frame.source_id.nunique()), "implant_T": s["implant_temp_C"],
                   "marginal_pp": s["marginal_effect_pp_per_sd"]["implant_temp_C"]}
    print(label, len(frame), sens[label]["implant_T"]["mean"], sens[label]["implant_T"]["ci95"], sens[label]["implant_T"]["p_positive"])
OUT["sensitivity"] = sens
json.dump(OUT, open(DATA / "hierarchical_results.json", "w"), indent=2)

grid_i = np.linspace(25, 1000, 40)
grid_a = np.linspace(1500, 1950, 40)
fig, ax = plt.subplots(figsize=(6, 4.6))
tt, cc = 30.0, 1e20
Z = np.zeros((len(grid_a), len(grid_i)))
bm = draws["beta"].mean(0)
for ia, a in enumerate(grid_a):
    for ii, i in enumerate(grid_i):
        x = (np.array([i, a, np.log10(tt), np.log10(cc)]) - mu) / sd
        Z[ia, ii] = 100 * inv_logit(bm[0] + x @ bm[1:])
cs = ax.contourf(grid_i, grid_a, Z, levels=np.linspace(0, 100, 11), cmap="viridis")
pts = al.drop_duplicates(["implant_temp_C", "anneal_temp_C"])
ax.scatter(pts.implant_temp_C, pts.anneal_temp_C, c="w", edgecolor="k", s=25)
plt.colorbar(cs, label="Predicted activation (%)")
ax.set_xlabel("Implant temperature (°C)")
ax.set_ylabel("Anneal temperature (°C)")
ax.set_title("Posterior-mean activation, 30 min anneal, 1e20 cm$^{-3}$ Al\n(white points: conditions in the data)")
plt.tight_layout()
plt.savefig(FIG / "design_map.png", dpi=150)
plt.savefig(FIG / "design_map.pdf")
print("DONE")
