import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import DATA, FIG
from models import FEATS, fit_hier, inv_logit, load_sic_al

al = load_sic_al(DATA / "database.csv")
db = pd.read_csv(DATA / "database.csv")
OUT = {}


def brief(r):
    return {"n_rows": r["n_rows"], "n_sources": r["n_sources"], "n_labs": r["n_labs"],
            "implant_T_logit": r["coef"]["implant_temp_C"], "anneal_T_logit": r["coef"]["anneal_temp_C"],
            "implant_T_pp": r["marginal_effect_pp_per_sd"]["implant_temp_C"], "anneal_T_pp": r["marginal_effect_pp_per_sd"]["anneal_temp_C"],
            "sd": r["sd"], "share": r["share"]}


main = fit_hier(al)
draws = main["_draws"]
mu, sd, beta = draws["mu"], draws["sd"], draws["beta"]
OUT["baseline_covariates"] = {"feature_means": dict(zip(FEATS, map(float, mu))), "feature_sd": dict(zip(FEATS, map(float, sd))),
                              "implant_T_contrast_C": float(sd[0]), "anneal_T_contrast_C": float(sd[1]),
                              "mean_anneal_time_min": float(10 ** mu[2]), "mean_implanted_conc_cm3": float(10 ** mu[3])}
pp = {}
for j, f in enumerate(FEATS):
    a, b = np.zeros(5), np.zeros(5)
    a[0] = b[0] = 1
    a[1 + j], b[1 + j] = 0.5, -0.5
    per_draw = 100 * (inv_logit(beta @ a) - inv_logit(beta @ b))
    plug_in = 100 * (inv_logit(beta.mean(0) @ a) - inv_logit(beta.mean(0) @ b))
    pp[f] = {"posterior_mean_of_transformed_prediction": float(per_draw.mean()), "credible_interval95": [float(np.percentile(per_draw, 2.5)), float(np.percentile(per_draw, 97.5))],
             "transformation_of_posterior_mean_coefficients": float(plug_in)}
OUT["percentage_point_effect_definition"] = pp

classes = sorted(al.definition_class.unique())
cls = {}
for c in [None] + classes:
    frame = al if c is None else al[al.definition_class != c]
    frame = frame.reset_index(drop=True)
    entry = {}
    for fac in ["both", "source"]:
        r = fit_hier(frame, factors=fac)
        entry[fac] = brief(r)
    entry["rows_dropped"] = int(len(al) - len(frame))
    entry["sources_dropped"] = int(al.source_id.nunique() - frame.source_id.nunique())
    cls["all" if c is None else "without_" + c] = entry
    print("class", c, entry["both"]["implant_T_pp"]["mean"], entry["source"]["share"]["between_groups_total"]["mean"])
for label, drop in [("without_resistivity_and_cv", ["resistivity_adjusted", "cv_activated_fraction"])]:
    frame = al[~al.definition_class.isin(drop)].reset_index(drop=True)
    entry = {fac: brief(fit_hier(frame, factors=fac)) for fac in ["both", "source"]}
    entry["rows_dropped"] = int(len(al) - len(frame))
    entry["sources_dropped"] = int(al.source_id.nunique() - frame.source_id.nunique())
    cls[label] = entry
OUT["leave_one_definition_class_out"] = cls

priors = {"IG(2,1)": (2.0, 1.0), "IG(0.5,0.5)": (0.5, 0.5), "IG(1,0.1)": (1.0, 0.1), "IG(3,2)": (3.0, 2.0)}
tab = {}
for name, (a0, b0) in priors.items():
    row = {}
    for fac in ["source", "both"]:
        r = fit_hier(al, factors=fac, a0=a0, b0=b0)
        row[fac] = {"share_between_groups": r["share"]["between_groups_total"], "implant_T_pp": r["marginal_effect_pp_per_sd"]["implant_temp_C"],
                    "anneal_T_pp": r["marginal_effect_pp_per_sd"]["anneal_temp_C"], "sd": r["sd"]}
    tab[name] = row
    print("prior", name, row["source"]["share_between_groups"]["mean"], row["both"]["share_between_groups"]["mean"])
OUT["variance_share_by_prior_and_structure"] = tab

kv = json.load(open(DATA / "kinetic_validation_results.json"))
fold = []
for s in al.source_id.unique():
    te, tr = al[al.source_id == s], al[al.source_id != s]
    row = {"held_out_source": s, "lab": str(te.lab_group.iloc[0]), "definition_class": str(te.definition_class.iloc[0]), "n_test_rows": int(len(te)), "n_train_rows": int(len(tr)),
           "n_train_sources": int(tr.source_id.nunique())}
    for f in FEATS:
        row["test_mean_" + f], row["train_mean_" + f] = float(te[f].mean()), float(tr[f].mean())
    row["test_activation_mean_pp"], row["train_activation_mean_pp"] = float(te.activation_pct.mean()), float(tr.activation_pct.mean())
    for m, v in kv["leave_one_source_out"].items():
        row["mae_" + m] = v["per_unit_mae_pp"][s]
    fold.append(row)
OUT["loso_fold_table"] = fold
OUT["loso_summary"] = {m: {"pooled_observation_weighted_mae_pp": v["mae_pp"], "source_weighted_mae_pp": v["unit_mean_mae_pp"], "held_out_units": v["paired_vs_mean"]["units"],
                           "units_better_than_training_mean": v["paired_vs_mean"]["units_better_than_mean"]} for m, v in kv["leave_one_source_out"].items()}
OUT["lolo_summary"] = {m: {"pooled_observation_weighted_mae_pp": v["mae_pp"], "source_weighted_mae_pp": v["unit_mean_mae_pp"], "held_out_units": v["paired_vs_mean"]["units"],
                           "units_better_than_training_mean": v["paired_vs_mean"]["units_better_than_mean"]} for m, v in kv["leave_one_lab_out"].items()}
OUT["temporal_summary"] = {m: {"mae_pp": v["mae_pp"], "n_test_rows": v["n_test"], "test_sources": kv["temporal"]["test_sources"]} for m, v in kv["temporal"].items() if isinstance(v, dict) and "mae_pp" in v}
ea = [r["selected_Ea_eV"] for r in kv["selected_Ea_by_fold_loso"]]
OUT["selected_Ea_loso"] = {"values_eV": ea, "range_eV": [min(ea), max(ea)], "at_upper_grid_bound": int(sum(e == 10 for e in ea)), "folds": len(ea)}

m = db[(db.source_id == "Nipoti2022_MSF") & db.activation_pct.notna()].copy()
m["sample"] = m.notes.str.extract(r"Sample ([^;]+);", expand=False)
w = m.pivot_table(index=["sample", "implanted_conc_cm3", "anneal_temp_C", "anneal_time_min"], columns="method_variant", values="activation_pct", aggfunc="first").dropna().reset_index()
w["diff_pp"] = w["p"] - w["p-mu"]
w["mean_pp"] = (w["p"] + w["p-mu"]) / 2
OUT["method_pairs"] = {"table": w.to_dict("records"), "n": int(len(w)), "mean_diff_pp": float(w.diff_pp.mean()), "sd_diff_pp": float(w.diff_pp.std(ddof=1)),
                       "spearman_diff_vs_conc": [float(x) for x in stats.spearmanr(w.diff_pp, np.log10(w.implanted_conc_cm3))],
                       "spearman_diff_vs_mean_activation": [float(x) for x in stats.spearmanr(w.diff_pp, w.mean_pp)],
                       "spearman_diff_vs_p_activation": [float(x) for x in stats.spearmanr(w.diff_pp, w["p"])],
                       "distinct_concentrations": sorted(float(x) for x in w.implanted_conc_cm3.unique()),
                       "mean_diff_by_conc": {str(k): float(v) for k, v in w.groupby("implanted_conc_cm3").diff_pp.mean().items()}}
fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.2))
for _, r in w.iterrows():
    ax[0].plot([0, 1], [r["p"], r["p-mu"]], "o-", color="0.35", lw=1)
ax[0].set_xticks([0, 1])
ax[0].set_xticklabels(["method p", "method p-$\\mu$"])
ax[0].set_ylabel("Reported activation (%)")
ax[0].set_title("Same nine samples, two fitting methods", fontsize=9)
sc = ax[1].scatter(w.mean_pp, w.diff_pp, c=np.log10(w.implanted_conc_cm3), cmap="viridis", edgecolor="k")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel("Mean of the two reported values (%)")
ax[1].set_ylabel("Difference, p minus p-$\\mu$ (pp)")
plt.colorbar(sc, ax=ax[1], label="log$_{10}$ implanted conc. (cm$^{-3}$)")
plt.tight_layout()
plt.savefig(FIG / "method_pairs.png", dpi=150)
plt.savefig(FIG / "method_pairs.pdf")

json.dump(OUT, open(DATA / "revision_v2_results.json", "w"), indent=2)
print("DONE")
