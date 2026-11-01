import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


RNG = np.random.RandomState(42)

df = pd.read_csv(DATA / "activation_dataset.csv")

sic = df[df["system"] == "SiC"].copy()
sic.loc[sic["source_id"] == "S1_S2", "implant_temp_C"] = 500
sic["implant_T"] = pd.to_numeric(sic["implant_temp_C"], errors="coerce")
sic["anneal_T"] = pd.to_numeric(sic["anneal_temp_C"], errors="coerce")
sic["act"] = pd.to_numeric(sic["activation_pct"], errors="coerce")
usable = sic.dropna(subset=["implant_T", "anneal_T", "act"]).copy()
print(f"SiC multivariate-usable n = {len(usable)}, sources = {sorted(usable['source_id'].unique())}")

X = usable[["anneal_T", "implant_T"]].to_numpy()
y = usable["act"].to_numpy()
sources = usable["source_id"].to_numpy()
unique_sources = sorted(set(sources))

scaler = StandardScaler()
Xs = scaler.fit_transform(X)
ols = LinearRegression().fit(Xs, y)
coef_anneal, coef_implant = ols.coef_
intercept = ols.intercept_
print(f"Point estimate (unchanged check): intercept={intercept:.3f}, anneal_T_std={coef_anneal:.3f}, implant_T_std={coef_implant:.3f}")

n_boot = 5000
boot_anneal, boot_implant = [], []
for _ in range(n_boot):
    resampled_sources = RNG.choice(unique_sources, size=len(unique_sources), replace=True)
    idx = np.concatenate([np.where(sources == s)[0] for s in resampled_sources])
    if len(set(sources[idx])) < 2:
        continue
    try:
        m = LinearRegression().fit(Xs[idx], y[idx])
        boot_anneal.append(m.coef_[0])
        boot_implant.append(m.coef_[1])
    except Exception:
        continue

ci_anneal = np.percentile(boot_anneal, [2.5, 97.5])
ci_implant = np.percentile(boot_implant, [2.5, 97.5])
print(f"Cluster-bootstrap CI (n_valid={len(boot_anneal)}): anneal_T={ci_anneal}, implant_T={ci_implant}")
print(f"SiC multivariate model: {len(usable)} rows / {len(unique_sources)} sources")

system_fittability = {}
for sysname in ["GaN", "Diamond", "GaAs", "InP", "Si"]:
    sub = df[df["system"] == sysname].copy()
    sub["anneal_T"] = pd.to_numeric(sub["anneal_temp_C"], errors="coerce")
    sub["act"] = pd.to_numeric(sub["activation_pct"], errors="coerce")
    us = sub.dropna(subset=["anneal_T", "act"])
    n_src = len(us["source_id"].unique())
    system_fittability[sysname] = {
        "n_total_rows": int(len(sub)),
        "n_usable_anneal_act": int(len(us)),
        "n_usable_sources": n_src,
        "fittable": bool(n_src >= 5 and len(us) >= 6),
    }
    print(f"{sysname}: total={len(sub)}, usable(anneal+act)={len(us)}, sources={n_src}, fittable={system_fittability[sysname]['fittable']}")

summary_rows = []
for sysname in ["SiC", "GaN", "Diamond", "GaAs", "InP", "Si"]:
    sub = df[df["system"] == sysname].copy()
    sub["act"] = pd.to_numeric(sub["activation_pct"], errors="coerce")
    sub["anneal_T_num"] = pd.to_numeric(sub["anneal_temp_C"], errors="coerce")
    anneal_vals = sub["anneal_T_num"].dropna()
    n_sources = len(sub["source_id"].unique())
    summary_rows.append({
        "system": sysname,
        "n_rows": int(len(sub)),
        "n_sources": int(n_sources),
        "best_activation_pct": float(sub["act"].max()) if sub["act"].notna().any() else None,
        "anneal_T_min": float(anneal_vals.min()) if len(anneal_vals) else None,
        "anneal_T_max": float(anneal_vals.max()) if len(anneal_vals) else None,
        "fittable_model": bool(sysname == "SiC"),
    })

for r in summary_rows:
    print(r)

results = {
    "sic_multivariate_model": {
        "n": len(usable),
        "n_sources": len(unique_sources),
        "sources": unique_sources,
        "point_estimate": {"intercept": float(intercept), "anneal_T_std": float(coef_anneal), "implant_T_std": float(coef_implant)},
        "cluster_bootstrap_ci": {
            "anneal_T_std": [float(ci_anneal[0]), float(ci_anneal[1])],
            "implant_T_std": [float(ci_implant[0]), float(ci_implant[1])],
            "n_valid_resamples": len(boot_anneal),
        },
        "implant_T_excludes_zero": bool(ci_implant[0] > 0 or ci_implant[1] < 0),
        "anneal_T_excludes_zero": bool(ci_anneal[0] > 0 or ci_anneal[1] < 0),
    },
    "system_fittability": system_fittability,
    "six_system_summary": summary_rows,
}
with open(DATA / "six_system_summary.json", "w") as f:
    json.dump(results, f, indent=2)
print("Saved data/six_system_summary.json")
