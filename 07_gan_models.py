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
gan = df[df["system"] == "GaN"].copy()
gan["anneal_T"] = pd.to_numeric(gan["anneal_temp_C"], errors="coerce")
gan["act"] = pd.to_numeric(gan["activation_pct"], errors="coerce")
gan["time"] = pd.to_numeric(gan["anneal_time_min"], errors="coerce")

usable = gan.dropna(subset=["anneal_T", "act"]).copy()
print(f"GaN univariate-usable n = {len(usable)}, sources = {sorted(usable['source_id'].unique())}")

def cluster_bootstrap_univariate(X, y, sources, n_boot=5000):
    unique_sources = sorted(set(sources))
    boot_coef = []
    for _ in range(n_boot):
        resampled = RNG.choice(unique_sources, size=len(unique_sources), replace=True)
        idx = np.concatenate([np.where(sources == s)[0] for s in resampled])
        if len(set(sources[idx])) < 2:
            continue
        Xb, yb = X[idx].reshape(-1, 1), y[idx]
        try:
            m = LinearRegression().fit(Xb, yb)
            boot_coef.append(m.coef_[0])
        except Exception:
            continue
    return np.array(boot_coef)

X = usable[["anneal_T"]].to_numpy()
y = usable["act"].to_numpy()
sources = usable["source_id"].to_numpy()
unique_sources = sorted(set(sources))

scaler = StandardScaler()
Xs = scaler.fit_transform(X)
ols = LinearRegression().fit(Xs, y)
point_coef = ols.coef_[0]
point_intercept = ols.intercept_
print(f"Model 1 (anneal_T only): point coef(std)={point_coef:.3f}, intercept={point_intercept:.3f}, n_sources={len(unique_sources)}")

boot1 = cluster_bootstrap_univariate(Xs[:, 0], y, sources)
ci1 = np.percentile(boot1, [2.5, 97.5]) if len(boot1) > 0 else [np.nan, np.nan]
print(f"Model 1 cluster-bootstrap CI (n_valid={len(boot1)}): anneal_T_std={ci1}")

loso_preds, loso_true = [], []
for s in unique_sources:
    train_idx = sources != s
    test_idx = sources == s
    if train_idx.sum() < 2:
        continue
    m = LinearRegression().fit(Xs[train_idx], y[train_idx])
    loso_preds.extend(m.predict(Xs[test_idx]))
    loso_true.extend(y[test_idx])
loso_mae = np.mean(np.abs(np.array(loso_preds) - np.array(loso_true))) if loso_preds else np.nan
print(f"Model 1 LOSO-CV (n_sources={len(unique_sources)}): MAE={loso_mae:.2f}")

usable_no_outlier = usable[usable["source_id"] != "Rocco2020_JAP"].copy()
Xo = usable_no_outlier[["anneal_T"]].to_numpy()
yo = usable_no_outlier["act"].to_numpy()
sources_o = usable_no_outlier["source_id"].to_numpy()
Xso = StandardScaler().fit_transform(Xo)
ols_o = LinearRegression().fit(Xso, yo)
boot1o = cluster_bootstrap_univariate(Xso[:, 0], yo, sources_o)
ci1o = np.percentile(boot1o, [2.5, 97.5]) if len(boot1o) > 0 else [np.nan, np.nan]
print(f"Model 1b (excl. Rocco2020 outlier, n={len(usable_no_outlier)}, sources={len(set(sources_o))}): coef={ols_o.coef_[0]:.3f}, CI={ci1o}")

usable2 = gan.dropna(subset=["anneal_T", "act", "time"]).copy()
print(f"\nModel 2 (anneal_T + time) usable n = {len(usable2)}, sources = {sorted(usable2['source_id'].unique())}")
X2 = usable2[["anneal_T", "time"]].to_numpy()
y2 = usable2["act"].to_numpy()
sources2 = usable2["source_id"].to_numpy()
unique_sources2 = sorted(set(sources2))
Xs2 = StandardScaler().fit_transform(X2)
ols2 = LinearRegression().fit(Xs2, y2)
r2_2 = ols2.score(Xs2, y2)
print(f"Model 2 point estimate: anneal_T_std={ols2.coef_[0]:.3f}, time_std={ols2.coef_[1]:.3f}, R2(in-sample)={r2_2:.3f}, n_sources={len(unique_sources2)}")

boot_anneal2, boot_time2 = [], []
for _ in range(5000):
    resampled = RNG.choice(unique_sources2, size=len(unique_sources2), replace=True)
    idx = np.concatenate([np.where(sources2 == s)[0] for s in resampled])
    if len(set(sources2[idx])) < 2:
        continue
    Xb, yb = Xs2[idx], y2[idx]
    try:
        m = LinearRegression().fit(Xb, yb)
        boot_anneal2.append(m.coef_[0])
        boot_time2.append(m.coef_[1])
    except Exception:
        continue
ci_anneal2 = np.percentile(boot_anneal2, [2.5, 97.5]) if boot_anneal2 else [np.nan, np.nan]
ci_time2 = np.percentile(boot_time2, [2.5, 97.5]) if boot_time2 else [np.nan, np.nan]
print(f"Model 2 cluster-bootstrap CI (n_valid={len(boot_anneal2)}): anneal_T_std={ci_anneal2}, time_std={ci_time2}")

sic_implant_ci = [0.0, 43.2]
sic_anneal_ci = None
try:
    with open(DATA / "sic_final_model_results.json") as f:
        sic_res = json.load(f)
    sic_anneal_ci = sic_res["cluster_bootstrap_ci"]["anneal_T_std"]
    sic_implant_ci = sic_res["cluster_bootstrap_ci"]["implant_T_std"]
    sic_anneal_point = sic_res["point_estimate"]["anneal_T_std"]
except Exception as e:
    print("Could not load SiC results:", e)
    sic_anneal_point = None

results = {
    "model1_univariate_annealT": {
        "n": len(usable), "n_sources": len(unique_sources), "sources": unique_sources,
        "point_coef_std": float(point_coef), "intercept": float(point_intercept),
        "cluster_bootstrap_ci": [float(ci1[0]), float(ci1[1])],
        "ci_excludes_zero": bool(ci1[0] > 0 or ci1[1] < 0),
        "loso_mae": float(loso_mae) if not np.isnan(loso_mae) else None,
        "n_valid_bootstrap": int(len(boot1)),
    },
    "model1b_sensitivity_excl_outlier": {
        "n": len(usable_no_outlier), "n_sources": len(set(sources_o)),
        "point_coef_std": float(ols_o.coef_[0]),
        "cluster_bootstrap_ci": [float(ci1o[0]), float(ci1o[1])],
        "ci_excludes_zero": bool(ci1o[0] > 0 or ci1o[1] < 0),
    },
    "model2_bivariate_annealT_time": {
        "n": len(usable2), "n_sources": len(unique_sources2), "sources": unique_sources2,
        "point_coef_annealT_std": float(ols2.coef_[0]), "point_coef_time_std": float(ols2.coef_[1]),
        "r2_insample": float(r2_2),
        "cluster_bootstrap_ci_annealT": [float(ci_anneal2[0]), float(ci_anneal2[1])],
        "cluster_bootstrap_ci_time": [float(ci_time2[0]), float(ci_time2[1])],
        "annealT_excludes_zero": bool(ci_anneal2[0] > 0 or ci_anneal2[1] < 0),
        "time_excludes_zero": bool(ci_time2[0] > 0 or ci_time2[1] < 0),
    },
    "cross_system_comparison": {
        "SiC_implant_T_ci": sic_implant_ci,
        "SiC_anneal_T_ci": sic_anneal_ci,
        "GaN_anneal_T_ci": [float(ci1[0]), float(ci1[1])],
        "note": "GaN sources rarely report implant temperature, so the GaN model uses anneal temperature (and anneal time) only.",
    },
}
with open(DATA / "gan_model_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("\nSaved data/gan_model_results.json")
