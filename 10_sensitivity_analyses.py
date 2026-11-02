import json

import numpy as np
import pandas as pd

from common import DATA
from models import FEATS, fit_hier, load_sic_al

db = pd.read_csv(DATA / "database.csv")
al = load_sic_al(DATA / "database.csv")
OUT = {}


def brief(r):
    return {"n_rows": r["n_rows"], "n_sources": r["n_sources"], "n_labs": r["n_labs"],
            "implant_T_logit": r["coef"]["implant_temp_C"], "anneal_T_logit": r["coef"]["anneal_temp_C"],
            "implant_T_pp": r["marginal_effect_pp_per_sd"]["implant_temp_C"],
            "anneal_T_pp": r["marginal_effect_pp_per_sd"]["anneal_temp_C"],
            "sd": r["sd"], "share": r["share"]}


clip = {}
for eps in [0.02, 0.01, 0.005]:
    lo, hi = 100 * eps, 100 * (1 - eps)
    n_lo = int((al.activation_pct < lo).sum())
    n_hi = int((al.activation_pct > hi).sum())
    n_edge = int(((al.activation_pct <= lo) | (al.activation_pct >= hi)).sum())
    affected = al[(al.activation_pct < lo) | (al.activation_pct > hi)]
    entry = {"lower_pct": lo, "upper_pct": hi, "rows_below_lower": n_lo, "rows_above_upper": n_hi,
             "rows_at_or_beyond_bounds": n_edge, "affected_sources": sorted(affected.source_id.unique()),
             "affected_values": sorted(float(v) for v in affected.activation_pct)}
    full = fit_hier(al, eps=eps)
    entry["all_rows"] = brief(full)
    nos = fit_hier(al[al.source_id != "Saks2004"].reset_index(drop=True), eps=eps)
    entry["without_Saks2004"] = brief(nos)
    clip[str(eps)] = entry
    print("eps", eps, n_lo, n_hi, entry["all_rows"]["implant_T_pp"], entry["without_Saks2004"]["implant_T_pp"],
          entry["all_rows"]["share"]["between_groups_total"]["mean"])
OUT["clipping"] = clip

vp = {}
for a0, b0 in [(2.0, 1.0), (0.5, 0.5), (1.0, 0.1), (3.0, 2.0)]:
    r = fit_hier(al, a0=a0, b0=b0)
    rs = np.random.default_rng(0).gamma(a0, 1.0, 200000)
    pr = np.sqrt(b0 / rs)
    vp[f"IG({a0},{b0})"] = {**brief(r), "prior_sd_median": float(np.median(pr)), "prior_sd_ci95": [float(np.percentile(pr, 2.5)), float(np.percentile(pr, 97.5))]}
    print("var prior", a0, b0, vp[f"IG({a0},{b0})"]["share"]["between_groups_total"]["mean"], r["sd"]["source"]["mean"], r["sd"]["lab"]["mean"])
OUT["variance_prior_sensitivity"] = vp

tab = (db[db.system == "SiC"].groupby("definition_class", dropna=False)
       .agg(rows=("db_id", "count"), sources=("source_id", "nunique")).reset_index())
modeled = al.groupby("definition_class").agg(rows=("db_id", "count"), sources=("source_id", "nunique")).reset_index()
OUT["definition_classes_all_systems"] = (db.groupby(["system", "definition_class"], dropna=False)
                                         .size().reset_index(name="rows").to_dict("records"))
OUT["definition_classes_sic_all"] = tab.to_dict("records")
OUT["definition_classes_modeled"] = modeled.to_dict("records")

hall_like = ["fit_acceptor_fraction", "derived_ratio", "net_of_compensation", "stated_hall"]
subsets = {
    "fit_acceptor_fraction_only": al[al.definition_class == "fit_acceptor_fraction"],
    "hall_based_only": al[al.definition_class.isin(hall_like)],
    "room_temperature_ionized_fraction": al[al.definition_class.isin(["stated_hall", "derived_ratio", "cv_activated_fraction"])],
    "without_resistivity_referenced": al[al.definition_class != "resistivity_adjusted"],
}
restricted = {}
for name, fr in subsets.items():
    fr = fr.reset_index(drop=True)
    entry = {"n_rows": len(fr), "n_sources": int(fr.source_id.nunique()),
             "n_labs": int(fr.lab_group.fillna(fr.source_id).nunique()),
             "classes": fr.definition_class.value_counts().to_dict(),
             "rows_per_source": fr.source_id.value_counts().to_dict(),
             "implant_T_levels": {s: sorted(float(v) for v in g.implant_temp_C.unique()) for s, g in fr.groupby("source_id")},
             "activation_range_pct": [float(fr.activation_pct.min()), float(fr.activation_pct.max())],
             "mean_logit_activation_by_source": {s: float(np.mean(np.log(np.clip(g.activation_pct, 2, 98) / 100 / (1 - np.clip(g.activation_pct, 2, 98) / 100)))) for s, g in fr.groupby("source_id")}}
    nlev = fr.implant_temp_C.nunique()
    entry["implant_T_distinct_levels"] = int(nlev)
    within = fr.groupby("source_id").implant_temp_C.nunique()
    entry["sources_with_implant_T_contrast"] = int((within > 1).sum())
    if len(fr) >= 10 and fr.source_id.nunique() >= 3 and (within > 1).sum() >= 1:
        feats = ["implant_temp_C", "anneal_temp_C", "log_time"]
        r = fit_hier(fr, feats=feats, factors="source")
        entry["model"] = {"features": feats, **brief(r)}
        entry["status"] = "fitted (exploratory)"
    else:
        entry["status"] = "not modelled: too few independent sources or no within-source implant-temperature contrast"
    restricted[name] = entry
    print(name, entry["n_rows"], entry["n_sources"], entry["status"])
OUT["restricted_definition_subsets"] = restricted

main = fit_hier(al)
d = main["_draws"]
s2e, s2u = d["s2e"], d["s2u"]
tot = s2u.sum(axis=1) + s2e
meas = {}
for label, sd_m in [("pm5pp", 0.2), ("pm10pp", 0.4), ("factor1.5", 0.8), ("factor2", 1.2)]:
    s2m = sd_m ** 2
    other = np.maximum(s2e - s2m, 0.0)
    meas[label] = {
        "assumed_measurement_sd_logit": sd_m,
        "share_between_groups": float(np.mean(s2u.sum(axis=1) / tot)),
        "share_measurement": float(np.mean(np.minimum(s2m, s2e) / tot)),
        "share_other_within_source": float(np.mean(other / tot)),
        "share_measurement_ci95": [float(np.percentile(np.minimum(s2m, s2e) / tot, 2.5)),
                                   float(np.percentile(np.minimum(s2m, s2e) / tot, 97.5))],
        "prob_measurement_variance_exceeds_residual": float(np.mean(s2m > s2e))}
OUT["measurement_error_scenarios"] = {
    "note": "assumed measurement SD on the logit scale; 0.2 is about 5 points and 0.4 about 10 points near 50% activation; 0.8 and 1.2 correspond to a factor of about 1.5 and 2 in the carrier fraction",
    "share_between_groups_mean": float(np.mean(s2u.sum(axis=1) / tot)),
    "share_between_groups_ci95": [float(np.percentile(s2u.sum(axis=1) / tot, 2.5)), float(np.percentile(s2u.sum(axis=1) / tot, 97.5))],
    "upper_bound_if_all_residual_is_measurement": {
        "between_groups_variance_fraction_of_total": float(np.mean(s2u.sum(axis=1) / tot)),
        "interpretation": "the between-group share is a property of the model decomposition; measurement error can only be taken out of the residual part"},
    "residual_sd_posterior_mean": float(np.mean(np.sqrt(s2e))),
    "scenarios": meas}

def flow(frame):
    steps = []
    f = frame
    steps.append(("Al-implanted 4H-SiC rows in the database", len(f), f.source_id.nunique()))
    for col, why in [("activation_pct", "no numeric activation"), ("implant_temp_C", "no implant temperature"),
                     ("anneal_temp_C", "no anneal temperature"), ("anneal_time_min", "no anneal time"),
                     ("implanted_conc_cm3", "no implanted concentration")]:
        keep = f.dropna(subset=[col])
        steps.append((f"after removing rows with {why}", len(keep), keep.source_id.nunique(), len(f) - len(keep)))
        f = keep
    keep = f[f.definition_class != "approximate"]
    steps.append(("after removing approximate values", len(keep), keep.source_id.nunique(), len(f) - len(keep)))
    f = keep
    keep = f[~f.two_step_anneal.fillna(False).astype(bool)]
    steps.append(("after removing two-step anneals", len(keep), keep.source_id.nunique(), len(f) - len(keep)))
    return steps

sic_al = db[(db.system == "SiC") & (db.dopant == "Al")]
OUT["dataset_flow"] = {
    "database_rows": int(len(db)), "database_sources": int(db.source_id.nunique()),
    "by_origin": db.groupby("origin").agg(rows=("db_id", "count"), sources=("source_id", "nunique")).reset_index().to_dict("records"),
    "by_system": db.groupby("system").agg(rows=("db_id", "count"), sources=("source_id", "nunique")).reset_index().to_dict("records"),
    "sic_by_dopant": db[db.system == "SiC"].groupby("dopant").agg(rows=("db_id", "count"), sources=("source_id", "nunique")).reset_index().to_dict("records"),
    "sic_al_steps": [list(s) for s in flow(sic_al)]}
init = pd.read_csv(DATA / "activation_dataset.csv", na_values=["NA"])
sic0 = init[init.system == "SiC"].copy()
for c in ["anneal_temp_C", "implant_temp_C", "activation_pct"]:
    sic0[c] = pd.to_numeric(sic0[c], errors="coerce")
sic0.loc[(sic0.source_id == "S1_S2") & sic0.implant_temp_C.isna(), "implant_temp_C"] = 500.0
a = sic0.dropna(subset=["anneal_temp_C", "activation_pct"])
b = a.dropna(subset=["implant_temp_C"])
OUT["dataset_flow"]["initial_sic_steps"] = [
    ["SiC rows in the initial compilation", len(sic0), int(sic0.source_id.nunique())],
    ["with numeric anneal temperature and activation", len(a), int(a.source_id.nunique())],
    ["with implant temperature (multivariate model)", len(b), int(b.source_id.nunique())]]
for s in OUT["dataset_flow"]["sic_al_steps"]:
    print(s)
for s in OUT["dataset_flow"]["initial_sic_steps"]:
    print(s)
json.dump(OUT, open(DATA / "sensitivity_results.json", "w"), indent=2)
print("DONE")
