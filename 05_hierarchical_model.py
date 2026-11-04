import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import DATA, FIG
from models import FEATS, fit_hier, inv_logit, load_sic_al

al = load_sic_al(DATA / "database.csv")
print("rows", len(al), "sources", al.source_id.nunique(), al.source_id.value_counts().to_dict())
lab = al.lab_group.fillna(al.source_id)


def strip(res):
    return {k: v for k, v in res.items() if k != "_draws"}


OUT = {"n_rows": len(al), "n_sources": int(al.source_id.nunique()), "n_labs": int(lab.nunique()),
       "rows_per_source": al.source_id.value_counts().to_dict(),
       "sources_per_lab": al.groupby(lab).source_id.nunique().to_dict(),
       "rows_per_lab": lab.value_counts().to_dict()}
main = fit_hier(al)
draws = main["_draws"]
OUT["main"] = strip(main)
OUT["main"]["sd_source"] = main["sd"]["source"]["mean"]
OUT["main"]["sd_residual"] = main["sd"]["residual"]["mean"]
for k in FEATS:
    print(k, main["coef"][k]["mean"], main["coef"][k]["ci95"], main["coef"][k]["rhat"], main["coef"][k]["ess"])
print("sd", {k: (v["mean"], v["ci95"], v["rhat"], v["ess"]) for k, v in main["sd"].items()})
print("share", {k: (v["mean"], v["ci95"]) for k, v in main["share"].items()})

rng = np.random.default_rng(0)
prior_var = 1.0 / rng.gamma(2.0, 1.0, size=200000)
OUT["prior_sd_inverse_gamma_2_1"] = {"median": float(np.median(np.sqrt(prior_var))),
                                     "ci95": [float(np.percentile(np.sqrt(prior_var), 2.5)),
                                              float(np.percentile(np.sqrt(prior_var), 97.5))]}

factor_models = {}
for f in ["source", "lab"]:
    r = fit_hier(al, factors=f)
    factor_models[f] = strip(r)
    print("factors", f, r["sd"], r["share"]["between_groups_total"]["mean"])
OUT["factor_variants"] = factor_models

sens = {}
variants = [
    ("drop_Saks2004", al[al.source_id != "Saks2004"], {}),
    ("drop_Saks2001", al[al.source_id != "Saks2001_JAP"], {}),
    ("drop_ECS92_implantT_is_range", al[al.source_id != "Nipoti2019_ECSTrans"], {}),
    ("drop_Pernot", al[al.source_id != "Pernot2005_JAP"], {}),
    ("no_lab_effect", al, {"factors": "source"}),
    ("prior_sd_2", al, {"beta_sd": 2.0}),
    ("prior_sd_10", al, {"beta_sd": 10.0}),
    ("without_conc", al, {"feats": ["implant_temp_C", "anneal_temp_C", "log_time"]})]
for label, frame, kw in variants:
    r = fit_hier(frame.reset_index(drop=True), **kw)
    sens[label] = {"n_rows": len(frame), "n_sources": int(frame.source_id.nunique()),
                   "implant_T": r["coef"]["implant_temp_C"], "marginal_pp": r["marginal_effect_pp_per_sd"]["implant_temp_C"],
                   "share_between_groups": r["share"]["between_groups_total"]}
    print(label, len(frame), sens[label]["marginal_pp"])
OUT["sensitivity"] = sens

nos = al[al.source_id != "Saks2004"]
OUT["implant_T_range"] = {
    "all": {"min": float(al.implant_temp_C.min()), "max": float(al.implant_temp_C.max()),
            "distinct": sorted(float(v) for v in al.implant_temp_C.unique()), "sd": float(al.implant_temp_C.std(ddof=1))},
    "without_Saks2004": {"min": float(nos.implant_temp_C.min()), "max": float(nos.implant_temp_C.max()),
                         "distinct": sorted(float(v) for v in nos.implant_temp_C.unique()),
                         "sd": float(nos.implant_temp_C.std(ddof=1)),
                         "rows_per_level": {str(k): int(v) for k, v in nos.implant_temp_C.value_counts().items()},
                         "sources_per_level": {str(k): sorted(nos[nos.implant_temp_C == k].source_id.unique())
                                               for k in sorted(nos.implant_temp_C.unique())}},
    "per_source": {s: sorted(float(v) for v in g.implant_temp_C.unique()) for s, g in al.groupby("source_id")}}
json.dump(OUT, open(DATA / "hierarchical_results.json", "w"), indent=2)

mu, sd = draws["mu"], draws["sd"]
beta = draws["beta"].mean(0)
grid_i = np.linspace(25, 1000, 40)
grid_a = np.linspace(1250, 1950, 50)
Z = np.zeros((len(grid_a), len(grid_i)))
for ia, a in enumerate(grid_a):
    for ii, i in enumerate(grid_i):
        x = (np.array([i, a, np.log10(30.0), np.log10(1e20)]) - mu) / sd
        Z[ia, ii] = 100 * inv_logit(beta[0] + x @ beta[1:])
fig, ax = plt.subplots(figsize=(5.2, 4.0))
cs = ax.contourf(grid_i, grid_a, Z, levels=np.linspace(0, 100, 11), cmap="viridis")
pts = al.drop_duplicates(["implant_temp_C", "anneal_temp_C"])
ax.scatter(pts.implant_temp_C, pts.anneal_temp_C, c="w", edgecolor="k", s=25)
plt.colorbar(cs, label="Predicted activation (%)")
ax.set_xlabel("Implant temperature (°C)")
ax.set_ylabel("Anneal temperature (°C)")
ax.set_title("Posterior-mean activation, 30 min anneal, $10^{20}$ cm$^{-3}$ Al\n(white points: conditions in the data)", fontsize=9)
plt.tight_layout()
plt.savefig(FIG / "design_map.png", dpi=150)
plt.savefig(FIG / "design_map.pdf")
print("DONE")
