import json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.linear_model import (LinearRegression, Ridge, Lasso, ElasticNet,
                                   BayesianRidge, ARDRegression, HuberRegressor)
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, AdaBoostRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVR

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


rng = np.random.default_rng(1234)
OUT = {}

df = pd.read_csv(DATA / "activation_dataset.csv", na_values=["NA"])
df["anneal_num"] = pd.to_numeric(df["anneal_temp_C"], errors="coerce")
df["implant_num"] = pd.to_numeric(df["implant_temp_C"], errors="coerce")
df["time_num"] = pd.to_numeric(df["anneal_time_min"], errors="coerce")

s1s2_mask = df["source_id"] == "S1_S2"
df.loc[s1s2_mask & df["implant_num"].isna(), "implant_num"] = 500.0

usable = df[df["activation_pct"].notna() & df["anneal_num"].notna()].copy()
sic_all = usable[usable.system == "SiC"].copy()
sic_mv = sic_all[sic_all.implant_num.notna()].copy()
n = len(sic_mv)
print(f"SiC multivariate-eligible n={n}, sources={sorted(sic_mv.source_id.unique())}")

X_raw = sic_mv[["anneal_num", "implant_num"]].to_numpy(dtype=float)
y = sic_mv["activation_pct"].to_numpy(dtype=float)
sources = sic_mv["source_id"].to_numpy()
X_mean, X_std = X_raw.mean(axis=0), X_raw.std(axis=0, ddof=1)
Xz = (X_raw - X_mean) / X_std

print("\n" + "="*78)
print("1. MODEL-FAMILY PANEL: implant-T and anneal-T effects across 17 model families")
print("="*78)

models = {
    "OLS": LinearRegression(),
    "Ridge(a=1)": Ridge(alpha=1.0),
    "Lasso(a=0.5)": Lasso(alpha=0.5, max_iter=5000),
    "ElasticNet": ElasticNet(alpha=0.5, l1_ratio=0.5, max_iter=5000),
    "BayesianRidge": BayesianRidge(),
    "ARDRegression": ARDRegression(),
    "HuberRegressor": HuberRegressor(max_iter=500),
    "PLSRegression(1comp)": PLSRegression(n_components=1),
    "DecisionTree(d=2)": DecisionTreeRegressor(max_depth=2, random_state=0),
    "DecisionTree(d=3)": DecisionTreeRegressor(max_depth=3, random_state=0),
    "RandomForest(n=50,d=3)": RandomForestRegressor(n_estimators=50, max_depth=3, random_state=0),
    "ExtraTrees(n=50,d=3)": ExtraTreesRegressor(n_estimators=50, max_depth=3, random_state=0),
    "AdaBoost(n=30)": AdaBoostRegressor(n_estimators=30, random_state=0),
    "KNN(k=3)": KNeighborsRegressor(n_neighbors=3),
    "KNN(k=5)": KNeighborsRegressor(n_neighbors=5),
    "GaussianProcess": GaussianProcessRegressor(kernel=RBF(length_scale=1.0)+WhiteKernel(noise_level=1.0),
                                                 normalize_y=True, random_state=0),
    "SVR(rbf)": SVR(kernel="rbf", C=10),
}

def permutation_importance_1d(model, X, y, col, n_perm=500, is_tree_rf=False):
    Xc = X.copy()
    other = 1 - col
    x_plus = np.zeros((1, 2)); x_plus[0, other] = 0.0; x_plus[0, col] = 1.0
    x_minus = np.zeros((1, 2)); x_minus[0, other] = 0.0; x_minus[0, col] = -1.0
    try:
        model.fit(Xc, y)
        eff_real = float(model.predict(x_plus)[0] - model.predict(x_minus)[0])
    except Exception:
        return None
    null_effs = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        try:
            model.fit(Xc, yp)
            e = float(model.predict(x_plus)[0] - model.predict(x_minus)[0])
            null_effs.append(e)
        except Exception:
            continue
    if len(null_effs) < 50:
        return {"effect_at_plus1_minus1SD": eff_real, "perm_p": None, "n_perm_valid": len(null_effs)}
    null_effs = np.array(null_effs)
    p = float(np.mean(np.abs(null_effs) >= np.abs(eff_real)))
    return {"effect_at_plus1_minus1SD": eff_real, "perm_p": p, "n_perm_valid": len(null_effs)}

panel_results = {}
for name, model in models.items():
    try:
        eff_anneal = permutation_importance_1d(model, Xz, y, col=0, n_perm=300)
        eff_implant = permutation_importance_1d(model, Xz, y, col=1, n_perm=300)
        panel_results[name] = {"anneal_T": eff_anneal, "implant_T": eff_implant}
        a_sign = "+" if eff_anneal and eff_anneal["effect_at_plus1_minus1SD"] > 0 else "-"
        i_sign = "+" if eff_implant and eff_implant["effect_at_plus1_minus1SD"] > 0 else "-"
        a_p = eff_anneal["perm_p"] if eff_anneal else None
        i_p = eff_implant["perm_p"] if eff_implant else None
        print(f"  {name:26s} anneal_T: sign={a_sign} eff={eff_anneal['effect_at_plus1_minus1SD']:7.2f} p={a_p}   "
              f"implant_T: sign={i_sign} eff={eff_implant['effect_at_plus1_minus1SD']:7.2f} p={i_p}")
    except Exception as e:
        panel_results[name] = {"error": str(e)}
        print(f"  {name:26s} FAILED: {e}")

n_models_run = sum(1 for v in panel_results.values() if "error" not in v)
implant_pos = sum(1 for v in panel_results.values()
                   if "error" not in v and v["implant_T"] and v["implant_T"]["effect_at_plus1_minus1SD"] > 0)
implant_sig = sum(1 for v in panel_results.values()
                   if "error" not in v and v["implant_T"] and v["implant_T"]["perm_p"] is not None and v["implant_T"]["perm_p"] < 0.1)
anneal_pos = sum(1 for v in panel_results.values()
                  if "error" not in v and v["anneal_T"] and v["anneal_T"]["effect_at_plus1_minus1SD"] > 0)
anneal_sig = sum(1 for v in panel_results.values()
                  if "error" not in v and v["anneal_T"] and v["anneal_T"]["perm_p"] is not None and v["anneal_T"]["perm_p"] < 0.1)

print(f"\nSUMMARY: {n_models_run} model families fit successfully.")
print(f"  implant_T effect positive-sign in {implant_pos}/{n_models_run} models; "
      f"permutation p<0.1 in {implant_sig}/{n_models_run} models.")
print(f"  anneal_T effect positive-sign in {anneal_pos}/{n_models_run} models; "
      f"permutation p<0.1 in {anneal_sig}/{n_models_run} models.")

OUT["robustness_panel"] = {
    "n_models_total": len(models), "n_models_run_successfully": n_models_run,
    "implant_T_positive_sign_count": implant_pos, "implant_T_perm_p_lt_0.1_count": implant_sig,
    "anneal_T_positive_sign_count": anneal_pos, "anneal_T_perm_p_lt_0.1_count": anneal_sig,
    "per_model_detail": panel_results,
}

print("\n" + "="*78)
print("2. SOURCE-LEVEL (CLUSTER) BOOTSTRAP")
print("="*78)
print("Resampling whole sources with replacement and refitting on all rows of the resampled sources.")

Xd_full = np.column_stack([np.ones(n), Xz])
def ols_fit(Xd, yv):
    beta, *_ = np.linalg.lstsq(Xd, yv, rcond=None)
    return beta

beta_full = ols_fit(Xd_full, y)
unique_sources = sorted(set(sources))
n_src = len(unique_sources)
src_to_rows = {s: np.where(sources == s)[0] for s in unique_sources}

n_boot = 5000
boot_betas_cluster = []
for _ in range(n_boot):
    samp_srcs = rng.choice(unique_sources, size=n_src, replace=True)
    row_idx = np.concatenate([src_to_rows[s] for s in samp_srcs])
    Xb = np.column_stack([np.ones(len(row_idx)), Xz[row_idx]])
    yb = y[row_idx]
    if np.linalg.matrix_rank(Xb) < 3 or len(set(sources[row_idx])) < 2:
        continue
    b = ols_fit(Xb, yb)
    boot_betas_cluster.append(b)
boot_betas_cluster = np.array(boot_betas_cluster)
ci_lo_c, ci_hi_c = np.percentile(boot_betas_cluster, [2.5, 97.5], axis=0)

print(f"\nCluster (source-level) bootstrap, {len(boot_betas_cluster)}/{n_boot} valid resamples "
      f"(resampling which of the {n_src} sources appear, not individual rows):")
names = ["intercept", "anneal_T (standardized)", "implant_T (standardized)"]
for i, nm in enumerate(names):
    excl = ci_lo_c[i] * ci_hi_c[i] > 0
    print(f"  {nm}: point={beta_full[i]:.2f}  CLUSTER CI [{ci_lo_c[i]:.2f}, {ci_hi_c[i]:.2f}]  "
          f"{'(excludes 0)' if excl else '(includes 0)'}")

anneal_excl_c = ci_lo_c[1]*ci_hi_c[1] > 0
implant_excl_c = ci_lo_c[2]*ci_hi_c[2] > 0

OUT["cluster_bootstrap"] = {
    "n_valid_cluster_resamples": len(boot_betas_cluster),
    "point_estimates": {"intercept": float(beta_full[0]), "anneal_T_std": float(beta_full[1]), "implant_T_std": float(beta_full[2])},
    "cluster_ci": {"intercept": [float(ci_lo_c[0]), float(ci_hi_c[0])],
                   "anneal_T_std": [float(ci_lo_c[1]), float(ci_hi_c[1])],
                   "implant_T_std": [float(ci_lo_c[2]), float(ci_hi_c[2])]},
    "anneal_T_excludes_zero_under_cluster_boot": bool(anneal_excl_c),
    "implant_T_excludes_zero_under_cluster_boot": bool(implant_excl_c),
    }

print("\n" + "="*78)
print("3. THERMAL-BUDGET COVARIATE (anneal T and time combined)")
print("="*78)
sic_time = sic_mv[sic_mv.time_num.notna()].copy()
print(f"SiC multivariate-eligible rows with anneal time reported: n={len(sic_time)} "
      f"(sources: {sorted(sic_time.source_id.unique())})")

if len(sic_time) >= 6 and sic_time.source_id.nunique() >= 2:
    T_K = sic_time.anneal_num.to_numpy(float) + 273.15
    t_min = sic_time.time_num.to_numpy(float)

    thermal_budget = T_K * np.log1p(t_min)
    y_tb = sic_time.activation_pct.to_numpy(float)
    implant_tb = sic_time.implant_num.to_numpy(float)

    r_raw_anneal, p_raw_anneal = stats.pearsonr(sic_time.anneal_num, y_tb)
    r_tb, p_tb = stats.pearsonr(thermal_budget, y_tb)
    print(f"\n  n={len(sic_time)}: raw anneal-T vs activation%: r={r_raw_anneal:.3f} (p={p_raw_anneal:.3f})")
    print(f"  thermal-budget (T_K * ln(1+t_min)) vs activation%: r={r_tb:.3f} (p={p_tb:.3f})")

    Xtb = np.column_stack([np.ones(len(sic_time)),
                            (thermal_budget - thermal_budget.mean())/thermal_budget.std(ddof=1),
                            (implant_tb - implant_tb.mean())/implant_tb.std(ddof=1)])
    beta_tb = ols_fit(Xtb, y_tb)
    resid_tb = y_tb - Xtb @ beta_tb
    r2_tb_model = 1 - np.sum(resid_tb**2)/np.sum((y_tb-y_tb.mean())**2)

    Xraw2 = np.column_stack([np.ones(len(sic_time)),
                              (sic_time.anneal_num.to_numpy(float) - sic_time.anneal_num.mean())/sic_time.anneal_num.std(ddof=1),
                              (implant_tb - implant_tb.mean())/implant_tb.std(ddof=1)])
    beta_raw2 = ols_fit(Xraw2, y_tb)
    resid_raw2 = y_tb - Xraw2 @ beta_raw2
    r2_raw2_model = 1 - np.sum(resid_raw2**2)/np.sum((y_tb-y_tb.mean())**2)

    print(f"\n  Multivariate (thermal-budget + implant-T) in-sample R^2={r2_tb_model:.3f}, "
          f"thermal-budget std coef={beta_tb[1]:.2f}, implant-T std coef={beta_tb[2]:.2f}")
    print(f"  Multivariate (raw anneal-T + implant-T), SAME n={len(sic_time)} subset, "
          f"in-sample R^2={r2_raw2_model:.3f}, anneal-T std coef={beta_raw2[1]:.2f}, "
          f"implant-T std coef={beta_raw2[2]:.2f}")
    improves = r2_tb_model > r2_raw2_model
    print(f"\n  Thermal-budget covariate "
          f"{'improves' if improves else 'does not improve'} in-sample R^2 relative to raw anneal-T "
          f"({r2_tb_model:.3f} vs {r2_raw2_model:.3f}, n={len(sic_time)}).")

    OUT["thermal_budget"] = {
        "n": len(sic_time), "sources": sorted(sic_time.source_id.unique().tolist()),
        "r_raw_anneal_T": float(r_raw_anneal), "r_thermal_budget": float(r_tb),
        "r2_model_thermal_budget": float(r2_tb_model), "r2_model_raw_annealT": float(r2_raw2_model),
        "thermal_budget_improves_fit": bool(improves),
        "caveat": "Anneal-time variation is limited to S3 (1/5/30 min); the 8 Saks2004 rows share a 5 min anneal.",
    }
else:
    print("  Too few rows with both implant-T and anneal time; thermal-budget model skipped.")
    OUT["thermal_budget"] = {"n": len(sic_time), "note": "insufficient for a defensible fit"}

print("\n" + "="*78)
print("4. HALL-FACTOR SENSITIVITY RE-FIT")
print("="*78)
print("Activation% of rows that do not state their r_H convention is divided by the geometric mean "
      "of the 2-5x overstatement range (about 3.16).")

correction_factor = np.sqrt(2.0 * 5.0)
sic_mv_corrected = sic_mv.copy()
not_stated_mask = sic_mv_corrected["hall_factor_convention"].isna()
print(f"  Rows not stating r_H convention (candidates for illustrative correction): "
      f"{int(not_stated_mask.sum())}/{len(sic_mv_corrected)}")
sic_mv_corrected.loc[not_stated_mask, "activation_pct"] = (
    sic_mv_corrected.loc[not_stated_mask, "activation_pct"] / correction_factor
)
sic_mv_corrected["activation_pct"] = sic_mv_corrected["activation_pct"].clip(upper=100)

y_corr = sic_mv_corrected["activation_pct"].to_numpy(float)
Xd_corr = np.column_stack([np.ones(n), Xz])
beta_corr = ols_fit(Xd_corr, y_corr)

boot_betas_corr = []
for _ in range(n_boot):
    samp_srcs = rng.choice(unique_sources, size=n_src, replace=True)
    row_idx = np.concatenate([src_to_rows[s] for s in samp_srcs])
    Xb = np.column_stack([np.ones(len(row_idx)), Xz[row_idx]])
    yb = y_corr[row_idx]
    if np.linalg.matrix_rank(Xb) < 3 or len(set(sources[row_idx])) < 2:
        continue
    b = ols_fit(Xb, yb)
    boot_betas_corr.append(b)
boot_betas_corr = np.array(boot_betas_corr)
ci_lo_corr, ci_hi_corr = np.percentile(boot_betas_corr, [2.5, 97.5], axis=0)

implant_excl_corr = ci_lo_corr[2]*ci_hi_corr[2] > 0
anneal_excl_corr = ci_lo_corr[1]*ci_hi_corr[1] > 0
print(f"\n  Under illustrative Hall-factor correction (r_H-unstated rows / {correction_factor:.2f}):")
print(f"  implant_T: point={beta_corr[2]:.2f}  CI [{ci_lo_corr[2]:.2f}, {ci_hi_corr[2]:.2f}]  "
      f"{'(excludes 0)' if implant_excl_corr else '(includes 0)'}")
print(f"  anneal_T: point={beta_corr[1]:.2f}  CI [{ci_lo_corr[1]:.2f}, {ci_hi_corr[1]:.2f}]  "
      f"{'(excludes 0)' if anneal_excl_corr else '(includes 0)'}")

OUT["hall_factor_sensitivity_refit"] = {
    "illustrative_correction_factor": float(correction_factor),
    "n_rows_corrected": int(not_stated_mask.sum()),
    "implant_T_coef_corrected": float(beta_corr[2]), "implant_T_ci_corrected": [float(ci_lo_corr[2]), float(ci_hi_corr[2])],
    "anneal_T_coef_corrected": float(beta_corr[1]), "anneal_T_ci_corrected": [float(ci_lo_corr[1]), float(ci_hi_corr[1])],
    "implant_T_excludes_zero": bool(implant_excl_corr),
    "note": "Uniform illustrative correction factor (geometric mean of the 2-5x range); not a per-row correction.",
}

fig, ax = plt.subplots(figsize=(8, 5))
model_names = list(panel_results.keys())
implant_effs = [panel_results[m]["implant_T"]["effect_at_plus1_minus1SD"] if "error" not in panel_results[m] and panel_results[m]["implant_T"] else np.nan for m in model_names]
anneal_effs = [panel_results[m]["anneal_T"]["effect_at_plus1_minus1SD"] if "error" not in panel_results[m] and panel_results[m]["anneal_T"] else np.nan for m in model_names]
x = np.arange(len(model_names))
width = 0.35
ax.bar(x - width/2, anneal_effs, width, label="anneal_T effect (+1 to -1 SD)", color="tab:blue")
ax.bar(x + width/2, implant_effs, width, label="implant_T effect (+1 to -1 SD)", color="tab:orange")
ax.axhline(0, color="k", linewidth=0.8)
ax.set_xticks(x); ax.set_xticklabels(model_names, rotation=60, ha="right", fontsize=8)
ax.set_ylabel("Effect on activation % (finite-difference, +1SD vs -1SD)")
ax.set_title(f"Robustness panel: implant-T vs anneal-T effect across {n_models_run} model families (n={n})")
ax.legend()
plt.tight_layout()
plt.savefig(FIG / "robustness_panel.png", dpi=150)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].hist(boot_betas_cluster[:, 1], bins=40, color="tab:blue", alpha=0.7, label="cluster boot")
axes[0].axvline(0, color="k", linestyle="--")
axes[0].set_title("Cluster bootstrap: anneal_T coefficient")
axes[1].hist(boot_betas_cluster[:, 2], bins=40, color="tab:orange", alpha=0.7, label="cluster boot")
axes[1].axvline(0, color="k", linestyle="--")
axes[1].set_title("Cluster bootstrap: implant_T coefficient")
plt.tight_layout()
plt.savefig(FIG / "sic_cluster_bootstrap_coefficients.png", dpi=150)
plt.close()

with open(DATA / "sic_bootstrap_robustness_results.json", "w") as f:
    json.dump(OUT, f, indent=2)
print("\nSaved data/sic_bootstrap_robustness_results.json")
print("Figures: figures/robustness_panel.png, figures/sic_cluster_bootstrap_coefficients.png")
print("\n" + "="*78)
print("DONE")
print("="*78)
