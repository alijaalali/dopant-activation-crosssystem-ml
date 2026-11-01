import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


rng = np.random.default_rng(42)
OUT = {}

df = pd.read_csv(DATA / "activation_dataset.csv", na_values=["NA"])
df["anneal_num"] = pd.to_numeric(df["anneal_temp_C"], errors="coerce")
df["implant_num"] = pd.to_numeric(df["implant_temp_C"], errors="coerce")
df["dose_num"] = pd.to_numeric(df["dose_cm-2"], errors="coerce")

s1s2_mask = df["source_id"] == "S1_S2"
df.loc[s1s2_mask & df["implant_num"].isna(), "implant_num"] = 500.0

print("=" * 78)
print("DATASET SUMMARY")
print("=" * 78)
for sys_ in ["SiC", "GaN", "Diamond"]:
    sub = df[df.system == sys_]
    usable = sub[sub.activation_pct.notna() & sub.anneal_num.notna()]
    print(f"{sys_}: total rows={len(sub)}, usable (anneal-T+act%)={len(usable)}, "
          f"independent source_ids={sub.source_id.nunique()}")

OUT["dataset_summary"] = {
    sys_: {
        "n_total_rows": int((df.system == sys_).sum()),
        "n_usable": int(((df.system == sys_) & df.activation_pct.notna() & df.anneal_num.notna()).sum()),
        "n_sources": int(df[df.system == sys_].source_id.nunique()),
    }
    for sys_ in ["SiC", "GaN", "Diamond"]
}

sic_mv = df[(df.system == "SiC") & df.activation_pct.notna() & df.anneal_num.notna()
            & df.implant_num.notna()].copy()
print(f"\nSiC multivariate-eligible n={len(sic_mv)}, sources={sorted(sic_mv.source_id.unique())}")

X = sic_mv[["anneal_num", "implant_num"]].to_numpy(dtype=float)
y = sic_mv["activation_pct"].to_numpy(dtype=float)
Xstd = (X - X.mean(axis=0)) / X.std(axis=0)
Xd = np.column_stack([np.ones(len(Xstd)), Xstd])

def ols_fit(Xd_, y_):
    beta, *_ = np.linalg.lstsq(Xd_, y_, rcond=None)
    resid = y_ - Xd_ @ beta
    ss_res = (resid ** 2).sum()
    ss_tot = ((y_ - y_.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return beta, r2

beta, r2 = ols_fit(Xd, y)
print(f"OLS (2-var, standardized): intercept={beta[0]:.2f}, "
      f"anneal_T_std={beta[1]:.2f}, implant_T_std={beta[2]:.2f}, R2={r2:.3f}")

sources = sic_mv.source_id.unique()
n_boot = 5000
boot_coefs = []
for _ in range(n_boot):
    chosen = rng.choice(sources, size=len(sources), replace=True)
    idx = np.concatenate([np.where(sic_mv.source_id.values == s)[0] for s in chosen])
    if len(np.unique(sic_mv.source_id.values[idx])) < 2:
        continue
    Xb, yb = Xd[idx], y[idx]
    if np.linalg.matrix_rank(Xb) < Xb.shape[1]:
        continue
    b, _ = ols_fit(Xb, yb)
    boot_coefs.append(b)
boot_coefs = np.array(boot_coefs)
ci_lo, ci_hi = np.percentile(boot_coefs, [2.5, 97.5], axis=0)

OUT["sic_multivariate_expanded"] = {
    "n": len(sic_mv),
    "n_sources": int(sic_mv.source_id.nunique()),
    "sources": sorted(sic_mv.source_id.unique().tolist()),
    "r2": float(r2),
    "coef_anneal_T_std": float(beta[1]),
    "coef_implant_T_std": float(beta[2]),
    "cluster_ci_anneal_T": [float(ci_lo[1]), float(ci_hi[1])],
    "cluster_ci_implant_T": [float(ci_lo[2]), float(ci_hi[2])],
    "anneal_T_excludes_zero": bool(ci_lo[1] > 0 or ci_hi[1] < 0),
    "implant_T_excludes_zero": bool(ci_lo[2] > 0 or ci_hi[2] < 0),
    "n_valid_boot": len(boot_coefs),
}
print(f"Cluster-bootstrap 95% CI: anneal_T=[{ci_lo[1]:.2f},{ci_hi[1]:.2f}], "
      f"implant_T=[{ci_lo[2]:.2f},{ci_hi[2]:.2f}]")

sic_mv["conc_num"] = pd.to_numeric(sic_mv["hole_conc_cm-3"], errors="coerce")

sic_dose = sic_mv[sic_mv.dose_num.notna()].copy()
print(f"\nSiC rows with numeric dose: n={len(sic_dose)}, sources={sorted(sic_dose.source_id.unique())}")
OUT["dose_covariate_check"] = {
    "n_with_numeric_dose": len(sic_dose),
    "sources": sorted(sic_dose.source_id.unique().tolist()),
    "note": "Dose is numeric only for part of the SiC rows; most sources report target concentration (cm-3), which is not convertible to dose without implant depth. Dose is therefore not used as a covariate.",
}

from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, BayesianRidge, ARDRegression, HuberRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.svm import SVR

models = {
    "OLS": LinearRegression(),
    "Ridge(a=1)": Ridge(alpha=1.0),
    "Lasso(a=0.5)": Lasso(alpha=0.5),
    "ElasticNet": ElasticNet(alpha=0.5, l1_ratio=0.5),
    "BayesianRidge": BayesianRidge(),
    "ARDRegression": ARDRegression(),
    "HuberRegressor": HuberRegressor(),
    "PLSRegression(1comp)": PLSRegression(n_components=1),
    "DecisionTree(d=2)": DecisionTreeRegressor(max_depth=2, random_state=0),
    "DecisionTree(d=3)": DecisionTreeRegressor(max_depth=3, random_state=0),
    "RandomForest(n=50,d=3)": RandomForestRegressor(n_estimators=50, max_depth=3, random_state=0),
    "ExtraTrees(n=50,d=3)": ExtraTreesRegressor(n_estimators=50, max_depth=3, random_state=0),
    "AdaBoost(n=30)": AdaBoostRegressor(n_estimators=30, random_state=0),
    "KNN(k=3)": KNeighborsRegressor(n_neighbors=3),
    "KNN(k=5)": KNeighborsRegressor(n_neighbors=5),
    "GaussianProcess": GaussianProcessRegressor(random_state=0, normalize_y=True),
    "SVR(rbf)": SVR(kernel="rbf"),
}

panel = {}
sd = Xstd.std(axis=0)
for name, model in models.items():
    try:
        model.fit(Xstd, y)
        n_perm = 300
        implant_effects, anneal_effects = [], []
        rng2 = np.random.default_rng(0)
        base_pred_hi_implant = model.predict(np.array([[0, 1]]))
        base_pred_lo_implant = model.predict(np.array([[0, -1]]))
        implant_effect_obs = float(base_pred_hi_implant[0] - base_pred_lo_implant[0])
        base_pred_hi_anneal = model.predict(np.array([[1, 0]]))
        base_pred_lo_anneal = model.predict(np.array([[-1, 0]]))
        anneal_effect_obs = float(base_pred_hi_anneal[0] - base_pred_lo_anneal[0])
        n_valid = 0
        count_ge_implant, count_ge_anneal = 0, 0
        for _ in range(n_perm):
            yp = rng2.permutation(y)
            try:
                m2 = model.__class__(**model.get_params())
            except Exception:
                m2 = model
            try:
                m2.fit(Xstd, yp)
                e_imp = float(m2.predict(np.array([[0, 1]]))[0] - m2.predict(np.array([[0, -1]]))[0])
                e_ann = float(m2.predict(np.array([[1, 0]]))[0] - m2.predict(np.array([[-1, 0]]))[0])
                if abs(e_imp) >= abs(implant_effect_obs):
                    count_ge_implant += 1
                if abs(e_ann) >= abs(anneal_effect_obs):
                    count_ge_anneal += 1
                n_valid += 1
            except Exception:
                continue
        panel[name] = {
            "implant_T": {"effect": implant_effect_obs, "perm_p": count_ge_implant / max(n_valid, 1), "n_perm_valid": n_valid},
            "anneal_T": {"effect": anneal_effect_obs, "perm_p": count_ge_anneal / max(n_valid, 1), "n_perm_valid": n_valid},
        }
    except Exception as e:
        panel[name] = {"error": str(e)}

n_ok = sum(1 for v in panel.values() if "error" not in v)
implant_pos = sum(1 for v in panel.values() if "error" not in v and v["implant_T"]["effect"] > 0)
implant_sig = sum(1 for v in panel.values() if "error" not in v and v["implant_T"]["perm_p"] < 0.1)
anneal_pos = sum(1 for v in panel.values() if "error" not in v and v["anneal_T"]["effect"] > 0)
anneal_sig = sum(1 for v in panel.values() if "error" not in v and v["anneal_T"]["perm_p"] < 0.1)

OUT["robustness_panel_expanded"] = {
    "n_models": len(models), "n_models_ok": n_ok,
    "implant_T_positive_sign_count": implant_pos, "implant_T_perm_p_lt_0.1_count": implant_sig,
    "anneal_T_positive_sign_count": anneal_pos, "anneal_T_perm_p_lt_0.1_count": anneal_sig,
    "per_model_detail": panel,
}
print(f"\nRobustness panel (expanded, n={len(sic_mv)}): implant-T positive-sign "
      f"{implant_pos}/{n_ok}, sig(p<0.1) {implant_sig}/{n_ok}; anneal-T positive-sign "
      f"{anneal_pos}/{n_ok}, sig(p<0.1) {anneal_sig}/{n_ok}")

dia = df[(df.system == "Diamond")]
dia_fitted = dia[dia.activation_pct.notna() & dia.anneal_num.notna()]
OUT["diamond_descriptive"] = {
    "n_total_rows": len(dia), "n_fitted_eligible": len(dia_fitted),
    "sources": sorted(dia.source_id.unique().tolist()),
    "note": "Only two diamond rows have both activation% and anneal-T; no model is fitted.",
}

with open(DATA / "sic_expanded_results.json", "w") as f:
    json.dump(OUT, f, indent=2)
print("\nSaved data/sic_expanded_results.json")

plt.rcParams.update({"font.size": 9})

fig, ax = plt.subplots(figsize=(5.5, 4.5))
sc = ax.scatter(sic_mv.anneal_num, sic_mv.implant_num, c=sic_mv.activation_pct,
                 s=90, cmap="viridis", edgecolor="k", linewidth=0.6, vmin=0, vmax=100)
for _, r in sic_mv.iterrows():
    ax.annotate(f"{r.activation_pct:.0f}%", (r.anneal_num, r.implant_num),
                fontsize=6.5, xytext=(3, 3), textcoords="offset points")
ax.set_xlabel("Anneal temperature (°C)")
ax.set_ylabel("Implant temperature (°C)")
ax.set_title(f"SiC: implant-T vs. anneal-T vs. activation% (n={len(sic_mv)})")
cb = plt.colorbar(sc, ax=ax)
cb.set_label("Activation (%)")
plt.tight_layout()
plt.savefig(FIG / "sic_implantT_annealT_activation.png", dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(6.5, 6))
names = list(panel.keys())
implant_vals = [panel[n]["implant_T"]["effect"] if "error" not in panel[n] else np.nan for n in names]
anneal_vals = [panel[n]["anneal_T"]["effect"] if "error" not in panel[n] else np.nan for n in names]
implant_p = [panel[n]["implant_T"]["perm_p"] if "error" not in panel[n] else np.nan for n in names]
y_pos = np.arange(len(names))
colors = ["#d62728" if (p is not None and not np.isnan(p) and p < 0.1) else "#7f7f7f" for p in implant_p]
ax.barh(y_pos - 0.2, implant_vals, height=0.35, color=colors, label="Implant-T effect (perm. sig. in red)")
ax.barh(y_pos + 0.2, anneal_vals, height=0.35, color="#1f77b4", alpha=0.6, label="Anneal-T effect")
ax.axvline(0, color="k", linewidth=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=7.5)
ax.set_xlabel("Effect on activation% (high vs. low, standardized ±1 SD)")
ax.set_title(f"Robustness panel: {len(models)} model families, n={len(sic_mv)}")
ax.legend(fontsize=7, loc="lower right")
plt.tight_layout()
plt.savefig(FIG / "robustness_panel.png", dpi=150)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
axes[0].hist(boot_coefs[:, 1], bins=40, color="#1f77b4", alpha=0.8)
axes[0].axvline(0, color="k", linestyle="--")
axes[0].set_title("Anneal-T coefficient\n(source-level cluster bootstrap)")
axes[0].set_xlabel("Standardized coefficient")
axes[1].hist(boot_coefs[:, 2], bins=40, color="#d62728", alpha=0.8)
axes[1].axvline(0, color="k", linestyle="--")
axes[1].set_title("Implant-T coefficient\n(source-level cluster bootstrap)")
axes[1].set_xlabel("Standardized coefficient")
plt.tight_layout()
plt.savefig(FIG / "sic_cluster_bootstrap_coefficients.png", dpi=150)
plt.close()

print("Figures updated: sic_implantT_annealT_activation.png, robustness_panel.png, "
      "sic_cluster_bootstrap_coefficients.png")
print("\nDONE")
