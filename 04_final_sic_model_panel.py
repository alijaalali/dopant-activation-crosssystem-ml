import json
import numpy as np
import pandas as pd
from sklearn.linear_model import (LinearRegression, Ridge, Lasso, ElasticNet,
                                    BayesianRidge, ARDRegression, HuberRegressor)
from sklearn.cross_decomposition import PLSRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.svm import SVR
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
print(f"n_sources = {len(unique_sources)}")

scaler = StandardScaler()
Xs = scaler.fit_transform(X)

ols = LinearRegression().fit(Xs, y)
coef_anneal, coef_implant = ols.coef_
intercept = ols.intercept_
print(f"Point estimate: intercept={intercept:.3f}, anneal_T_std={coef_anneal:.3f}, implant_T_std={coef_implant:.3f}")

n_boot = 5000
boot_anneal, boot_implant, boot_intercept = [], [], []
for _ in range(n_boot):
    resampled_sources = RNG.choice(unique_sources, size=len(unique_sources), replace=True)
    idx = np.concatenate([np.where(sources == s)[0] for s in resampled_sources])
    if len(set(sources[idx])) < 2:
        continue
    Xb, yb = Xs[idx], y[idx]
    try:
        m = LinearRegression().fit(Xb, yb)
        boot_anneal.append(m.coef_[0])
        boot_implant.append(m.coef_[1])
        boot_intercept.append(m.intercept_)
    except Exception:
        continue

ci_anneal = np.percentile(boot_anneal, [2.5, 97.5])
ci_implant = np.percentile(boot_implant, [2.5, 97.5])
print(f"Cluster-bootstrap CI (n_valid={len(boot_anneal)}): anneal_T={ci_anneal}, implant_T={ci_implant}")

loso_preds, loso_true = [], []
for s in unique_sources:
    train_idx = sources != s
    test_idx = sources == s
    if train_idx.sum() < 3:
        continue
    m = LinearRegression().fit(Xs[train_idx], y[train_idx])
    loso_preds.extend(m.predict(Xs[test_idx]))
    loso_true.extend(y[test_idx])
loso_mae = np.mean(np.abs(np.array(loso_preds) - np.array(loso_true)))
loso_rmse = np.sqrt(np.mean((np.array(loso_preds) - np.array(loso_true)) ** 2))
print(f"LOSO-CV (n_sources={len(unique_sources)}): MAE={loso_mae:.2f}, RMSE={loso_rmse:.2f}")

models = {
    "OLS": LinearRegression(),
    "Ridge(a=1)": Ridge(alpha=1.0),
    "Lasso(a=0.5)": Lasso(alpha=0.5, max_iter=5000),
    "ElasticNet": ElasticNet(alpha=0.5, l1_ratio=0.5, max_iter=5000),
    "BayesianRidge": BayesianRidge(),
    "ARDRegression": ARDRegression(),
    "HuberRegressor": HuberRegressor(max_iter=500),
    "PLSRegression(1comp)": PLSRegression(n_components=1),
    "DecisionTree(d=2)": DecisionTreeRegressor(max_depth=2, random_state=42),
    "DecisionTree(d=3)": DecisionTreeRegressor(max_depth=3, random_state=42),
    "RandomForest(n=50,d=3)": RandomForestRegressor(n_estimators=50, max_depth=3, random_state=42),
    "ExtraTrees(n=50,d=3)": ExtraTreesRegressor(n_estimators=50, max_depth=3, random_state=42),
    "AdaBoost(n=30)": AdaBoostRegressor(n_estimators=30, random_state=42),
    "KNN(k=3)": KNeighborsRegressor(n_neighbors=3),
    "KNN(k=5)": KNeighborsRegressor(n_neighbors=5),
    "GaussianProcess": GaussianProcessRegressor(random_state=42, normalize_y=True),
    "SVR(rbf)": SVR(kernel="rbf"),
}

panel_results = {}
n_perm = 300
for name, model in models.items():
    try:
        model.fit(Xs, y)
        if hasattr(model, "coef_"):
            coefs = np.ravel(model.coef_)
            eff_anneal_sign = np.sign(coefs[0]) if len(coefs) > 0 else np.nan
            eff_implant_sign = np.sign(coefs[1]) if len(coefs) > 1 else np.nan
        else:

            base_pred = model.predict(Xs)
            Xs_perm_anneal = Xs.copy(); Xs_perm_anneal[:, 0] = Xs[:, 0].mean()
            Xs_perm_implant = Xs.copy(); Xs_perm_implant[:, 1] = Xs[:, 1].mean()
            eff_anneal_sign = np.sign(np.mean(base_pred) - np.mean(model.predict(Xs_perm_anneal)))
            eff_implant_sign = np.sign(np.mean(base_pred) - np.mean(model.predict(Xs_perm_implant)))

        obs_implant_effect = model.predict(np.array([[0, 1], [0, -1]]) * Xs[:, 1].std())
        obs_implant_effect = obs_implant_effect[0] - obs_implant_effect[1]
        obs_anneal_effect_arr = model.predict(np.array([[1, 0], [-1, 0]]) * Xs[:, 0].std())
        obs_anneal_effect = obs_anneal_effect_arr[0] - obs_anneal_effect_arr[1]

        count_implant, count_anneal, n_valid = 0, 0, 0
        for _ in range(n_perm):
            yp = RNG.permutation(y)
            try:
                mp = type(model)(**model.get_params()) if hasattr(model, "get_params") else model
                mp.fit(Xs, yp)
                pi = mp.predict(np.array([[0, 1], [0, -1]]) * Xs[:, 1].std())
                pa = mp.predict(np.array([[1, 0], [-1, 0]]) * Xs[:, 0].std())
                if abs(pi[0] - pi[1]) >= abs(obs_implant_effect):
                    count_implant += 1
                if abs(pa[0] - pa[1]) >= abs(obs_anneal_effect):
                    count_anneal += 1
                n_valid += 1
            except Exception:
                continue
        perm_p_implant = count_implant / n_valid if n_valid else np.nan
        perm_p_anneal = count_anneal / n_valid if n_valid else np.nan

        panel_results[name] = {
            "implant_T_sign": int(np.sign(obs_implant_effect)),
            "implant_T_perm_p": perm_p_implant,
            "anneal_T_sign": int(np.sign(obs_anneal_effect)),
            "anneal_T_perm_p": perm_p_anneal,
        }
    except Exception as e:
        panel_results[name] = {"error": str(e)}

n_implant_pos = sum(1 for v in panel_results.values() if v.get("implant_T_sign", 0) > 0)
n_implant_sig = sum(1 for v in panel_results.values() if v.get("implant_T_perm_p", 1) < 0.1)
n_anneal_pos = sum(1 for v in panel_results.values() if v.get("anneal_T_sign", 0) > 0)
n_anneal_sig = sum(1 for v in panel_results.values() if v.get("anneal_T_perm_p", 1) < 0.1)
print(f"Robustness panel: implant-T positive {n_implant_pos}/{len(models)}, sig(p<0.1) {n_implant_sig}/{len(models)}")
print(f"Robustness panel: anneal-T positive {n_anneal_pos}/{len(models)}, sig(p<0.1) {n_anneal_sig}/{len(models)}")

gan = df[df["system"] == "GaN"].copy()
gan["anneal_T"] = pd.to_numeric(gan["anneal_temp_C"], errors="coerce")
gan["act"] = pd.to_numeric(gan["activation_pct"], errors="coerce")
gan_usable = gan.dropna(subset=["anneal_T", "act"])
print(f"GaN usable (anneal-T, activation%) rows: {len(gan_usable)}, sources: {sorted(gan_usable['source_id'].unique())}")

results = {
    "sic_multivariate_n": len(usable),
    "sic_multivariate_sources": sorted(usable["source_id"].unique().tolist()),
    "point_estimate": {"intercept": float(intercept), "anneal_T_std": float(coef_anneal), "implant_T_std": float(coef_implant)},
    "cluster_bootstrap_ci": {
        "anneal_T_std": [float(ci_anneal[0]), float(ci_anneal[1])],
        "implant_T_std": [float(ci_implant[0]), float(ci_implant[1])],
        "n_valid_resamples": len(boot_anneal),
    },
    "anneal_T_excludes_zero": bool(ci_anneal[0] > 0 or ci_anneal[1] < 0),
    "implant_T_excludes_zero": bool(ci_implant[0] > 0 or ci_implant[1] < 0),
    "loso_cv": {"n_sources": len(unique_sources), "mae": float(loso_mae), "rmse": float(loso_rmse)},
    "robustness_panel": {
        "n_models": len(models),
        "implant_T_positive_count": n_implant_pos,
        "implant_T_sig_count": n_implant_sig,
        "anneal_T_positive_count": n_anneal_pos,
        "anneal_T_sig_count": n_anneal_sig,
        "per_model": panel_results,
    },
    "gan_n_usable": len(gan_usable),
    "gan_sources": sorted(gan_usable["source_id"].unique().tolist()),
}

with open(DATA / "sic_final_model_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Saved data/sic_final_model_results.json")
