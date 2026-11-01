import itertools
import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn.base import clone
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import AdaBoostRegressor, ExtraTreesRegressor, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.linear_model import (ARDRegression, BayesianRidge, ElasticNet, HuberRegressor,
                                  Lasso, LinearRegression, Ridge)
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor

from common import DATA, FIG, cluster_enumeration, load, ols, r_squared, reml_random_intercept, \
    sic_model_frame, summarize_enumeration, zscore

warnings.filterwarnings("ignore")
OUT = {}
NAMES = ["anneal_T", "implant_T"]

df = load()
sic_all = df[(df.system == "SiC")].dropna(subset=["anneal_T", "act"]).reset_index(drop=True)
mv = sic_model_frame(df)
n = len(mv)
src = mv.source_id.to_numpy()
y = mv.act.to_numpy(float)
Xraw = mv[["anneal_T", "implant_T"]].to_numpy(float)
Xz = zscore(Xraw)
uniq = sorted(set(src))
OUT["data"] = {"n_sic_rows_total": int((df.system == "SiC").sum()),
               "n_sic_rows_anneal_and_activation": len(sic_all),
               "n_model_rows": n, "model_sources": uniq,
               "rows_per_source": {s: int((src == s).sum()) for s in uniq},
               "anneal_T_range_model": [float(Xraw[:, 0].min()), float(Xraw[:, 0].max())],
               "implant_T_range_model": [float(Xraw[:, 1].min()), float(Xraw[:, 1].max())]}
print("model rows", n, OUT["data"]["rows_per_source"])

beta = ols(Xz, y)
r2, r2adj = r_squared(Xz, y)
OUT["ols"] = {"intercept": float(beta[0]), "anneal_T": float(beta[1]), "implant_T": float(beta[2]),
              "r2": float(r2), "r2_adj": float(r2adj)}
print("OLS", OUT["ols"])

enum = cluster_enumeration(src, Xz, y)
OUT["cluster_bootstrap_exact"] = summarize_enumeration(enum, NAMES)
print("cluster bootstrap", json.dumps(OUT["cluster_bootstrap_exact"]["terms"]))

rng = np.random.default_rng(42)
row_betas = []
for _ in range(5000):
    s = rng.choice(n, size=n, replace=True)
    A = np.column_stack([np.ones(n), Xz[s]])
    if np.linalg.matrix_rank(A) < 3:
        continue
    row_betas.append(np.linalg.lstsq(A, y[s], rcond=None)[0])
row_betas = np.array(row_betas)
OUT["row_level_bootstrap_for_comparison"] = {
    nm: [float(v) for v in np.percentile(row_betas[:, j + 1], [2.5, 97.5])] for j, nm in enumerate(NAMES)}

loso = {}
for s in uniq:
    m = src != s
    b = ols(Xz[m], y[m])
    A = np.column_stack([np.ones(m.sum()), Xz[m]])
    loso[s] = {"n_rows_remaining": int(m.sum()), "rank": int(np.linalg.matrix_rank(A)),
               "anneal_T": float(b[1]), "implant_T": float(b[2])}
OUT["leave_one_source_out_coefficients"] = loso
print("LOSO coefficients", loso)

def subset_fit(mask, label):
    m = np.asarray(mask)
    A = np.column_stack([np.ones(m.sum()), Xz[m]])
    rank = int(np.linalg.matrix_rank(A))
    b = np.linalg.lstsq(A, y[m], rcond=None)[0]
    entry = {"n_rows": int(m.sum()), "sources": sorted(set(src[m])), "design_rank": rank,
             "anneal_T": float(b[1]), "implant_T": float(b[2]),
             "implant_T_distinct_values": sorted(set(Xraw[m, 1].tolist()))}
    return label, entry


sub = {}
for label, mask in [
        ("all_four_sources", np.ones(n, bool)),
        ("without_Weisse_net_of_compensation_definition", src != "Weisse2019_AIPAdv"),
        ("without_Negoro_derived_ratio_assumed_depth", src != "S3"),
        ("without_Negoro_and_Weisse", np.isin(src, ["S1_S2", "Saks2004"])),
        ("Saks_only", src == "Saks2004")]:
    k, v = subset_fit(mask, label)
    sub[k] = v
sk = src == "Saks2004"
sub["Saks_only"]["mean_activation_by_implant_T"] = {
    str(int(t)): float(y[sk & (Xraw[:, 1] == t)].mean()) for t in sorted(set(Xraw[sk, 1]))}
hall_direct = src == "S1_S2"
sub["stated_Hall_rows_only_S1_S2"] = {"n_rows": int(hall_direct.sum()),
                                     "implant_T_distinct_values": sorted(set(Xraw[hall_direct, 1].tolist())),
                                     "note": "implant temperature is constant, so the coefficient is not estimable"}
D = np.column_stack([(src == s).astype(float) for s in uniq] + [Xz])
bfe = np.linalg.lstsq(D, y, rcond=None)[0]
sub["source_fixed_effects"] = {"anneal_T": float(bfe[len(uniq)]), "implant_T": float(bfe[len(uniq) + 1]),
                               "design_rank": int(np.linalg.matrix_rank(D))}
OUT["definition_and_source_subsets"] = sub
print("subsets", json.dumps(sub))

fit_mixed = reml_random_intercept(Xz, y, src)
tcrit = stats.t.ppf(0.975, df=len(uniq) - 1)
OUT["random_intercept_model"] = {
    "var_ratio_source_to_residual": fit_mixed["var_ratio_source_to_residual"],
    "df_used_for_ci": len(uniq) - 1,
    **{nm: {"coef": float(fit_mixed["beta"][j + 1]), "se": float(fit_mixed["se"][j + 1]),
            "ci95_approx": [float(fit_mixed["beta"][j + 1] - tcrit * fit_mixed["se"][j + 1]),
                            float(fit_mixed["beta"][j + 1] + tcrit * fit_mixed["se"][j + 1])]}
       for j, nm in enumerate(NAMES)}}
print("random intercept", OUT["random_intercept_model"])


def loso_mae(Xcols, yy, groups):
    preds, truth = [], []
    for s in sorted(set(groups)):
        tr, te = groups != s, groups == s
        mu, sd = Xcols[tr].mean(0), Xcols[tr].std(0, ddof=1)
        b = ols((Xcols[tr] - mu) / sd, yy[tr])
        preds.extend(np.column_stack([np.ones(te.sum()), (Xcols[te] - mu) / sd]) @ b)
        truth.extend(yy[te])
    return float(np.mean(np.abs(np.array(preds) - np.array(truth))))


ua = sic_all.anneal_T.to_numpy(float)[:, None]
OUT["cross_validation"] = {
    "multivariate_mae_4_sources": loso_mae(Xraw, y, src),
    "univariate_anneal_only_mae": loso_mae(ua, sic_all.act.to_numpy(float), sic_all.source_id.to_numpy()),
    "univariate_n_rows": len(sic_all), "univariate_n_sources": int(sic_all.source_id.nunique()),
    "mean_baseline_mae_4_sources": float(np.mean(np.concatenate(
        [np.abs(y[src == s] - y[src != s].mean()) for s in uniq])))}
r_uv, p_uv = stats.pearsonr(sic_all.anneal_T, sic_all.act)
OUT["univariate_anneal_correlation"] = {"n": len(sic_all), "r": float(r_uv), "p": float(p_uv)}
print("CV", OUT["cross_validation"], OUT["univariate_anneal_correlation"])

tb = mv.dropna(subset=["time"]).reset_index(drop=True)
Bud = (tb.anneal_T.to_numpy(float) + 273.15) * np.log1p(tb.time.to_numpy(float))
imp = tb.implant_T.to_numpy(float)
ytb = tb.act.to_numpy(float)
r2_b, adj_b = r_squared(zscore(np.column_stack([Bud, imp])), ytb)
r2_a, adj_a = r_squared(zscore(np.column_stack([tb.anneal_T.to_numpy(float), imp])), ytb)
OUT["thermal_budget"] = {"n": len(tb), "sources": sorted(tb.source_id.unique()),
                         "r2_budget_model": float(r2_b), "r2_adj_budget_model": float(adj_b),
                         "r2_raw_anneal_model": float(r2_a), "r2_adj_raw_anneal_model": float(adj_a)}
print("thermal budget", OUT["thermal_budget"])

sc = np.array([1.0, 5.0])
grid = []
for f in itertools.product(sc, repeat=len(uniq)):
    fac = np.array([f[uniq.index(s)] for s in src])
    b = ols(Xz, y / fac)
    grid.append({"factors": dict(zip(uniq, map(float, f))), "implant_T": float(b[2]), "anneal_T": float(b[1])})
rng = np.random.default_rng(7)
mc = []
for _ in range(5000):
    f = np.exp(rng.uniform(0, np.log(5.0), size=len(uniq)))
    fac = np.array([f[uniq.index(s)] for s in src])
    mc.append(ols(Xz, y / fac)[2])
mc = np.array(mc)
imp_grid = np.array([g["implant_T"] for g in grid])
OUT["hall_factor_sensitivity"] = {
    "design": "per-source overstatement factor f_s applied as activation/f_s; grid f_s in {1,5}; Monte Carlo f_s log-uniform on [1,5]",
    "uniform_scaling_ci_note": "a single factor applied to all rows rescales coefficients and intervals by 1/f and cannot change their sign or significance",
    "grid_n": len(grid), "grid_implant_T_min": float(imp_grid.min()), "grid_implant_T_max": float(imp_grid.max()),
    "grid_n_positive": int((imp_grid > 0).sum()),
    "grid_worst": min(grid, key=lambda g: g["implant_T"]),
    "mc_implant_T_percentiles_2.5_50_97.5": [float(v) for v in np.percentile(mc, [2.5, 50, 97.5])],
    "mc_share_positive": float(np.mean(mc > 0))}
print("hall sensitivity", {k: v for k, v in OUT["hall_factor_sensitivity"].items() if k != "design"})

models = {
    "OLS": LinearRegression(), "Ridge(a=1)": Ridge(alpha=1.0), "Lasso(a=0.5)": Lasso(alpha=0.5),
    "ElasticNet": ElasticNet(alpha=0.5, l1_ratio=0.5), "BayesianRidge": BayesianRidge(),
    "ARDRegression": ARDRegression(), "HuberRegressor": HuberRegressor(max_iter=500),
    "PLSRegression(1comp)": PLSRegression(n_components=1),
    "DecisionTree(d=2)": DecisionTreeRegressor(max_depth=2, random_state=0),
    "DecisionTree(d=3)": DecisionTreeRegressor(max_depth=3, random_state=0),
    "RandomForest(n=50,d=3)": RandomForestRegressor(n_estimators=50, max_depth=3, random_state=0),
    "ExtraTrees(n=50,d=3)": ExtraTreesRegressor(n_estimators=50, max_depth=3, random_state=0),
    "AdaBoost(n=30)": AdaBoostRegressor(n_estimators=30, random_state=0),
    "KNN(k=3)": KNeighborsRegressor(n_neighbors=3), "KNN(k=5)": KNeighborsRegressor(n_neighbors=5),
    "GaussianProcess": GaussianProcessRegressor(kernel=RBF(1.0) + WhiteKernel(1.0), normalize_y=True, random_state=0),
    "SVR(rbf)": SVR(kernel="rbf", C=10),
}
PLUS = {"anneal_T": (np.array([[1.0, 0.0]]), np.array([[-1.0, 0.0]])),
        "implant_T": (np.array([[0.0, 1.0]]), np.array([[0.0, -1.0]]))}


def effects(model, Xf, yf):
    m = clone(model).fit(Xf, yf)
    return {k: float(np.ravel(m.predict(hi))[0] - np.ravel(m.predict(lo))[0]) for k, (hi, lo) in PLUS.items()}


rng = np.random.default_rng(0)
panel = {}
N_PERM = 500
for name, model in models.items():
    full = effects(model, Xz, y)
    null = {k: [] for k in PLUS}
    for _ in range(N_PERM):
        e = effects(model, Xz, rng.permutation(y))
        for k in PLUS:
            null[k].append(abs(e[k]))
    entry = {"effect": full}
    entry["row_permutation_p"] = {k: float(np.mean(np.array(null[k]) >= abs(full[k]))) for k in PLUS}
    drop = {}
    for s in uniq:
        m = src != s
        drop[s] = effects(model, Xz[m], y[m])
    entry["drop_source_effect"] = drop
    entry["implant_T_positive_in_all_drop_source_fits"] = bool(all(d["implant_T"] > 0 for d in drop.values()))
    panel[name] = entry
    print(f"{name:24s} implant={full['implant_T']:7.2f} p={entry['row_permutation_p']['implant_T']:.3f} "
          f"anneal={full['anneal_T']:7.2f} p={entry['row_permutation_p']['anneal_T']:.3f} "
          f"stable={entry['implant_T_positive_in_all_drop_source_fits']}")

summ = {
    "n_models": len(panel),
    "implant_T_positive": sum(v["effect"]["implant_T"] > 0 for v in panel.values()),
    "implant_T_row_perm_p_lt_0.1": sum(v["row_permutation_p"]["implant_T"] < 0.1 for v in panel.values()),
    "anneal_T_positive": sum(v["effect"]["anneal_T"] > 0 for v in panel.values()),
    "anneal_T_row_perm_p_lt_0.1": sum(v["row_permutation_p"]["anneal_T"] < 0.1 for v in panel.values()),
    "implant_T_positive_in_all_drop_source_fits": sum(v["implant_T_positive_in_all_drop_source_fits"] for v in panel.values()),
    "implant_T_positive_without_Saks2004": sum(v["drop_source_effect"]["Saks2004"]["implant_T"] > 0 for v in panel.values()),
    "n_permutations": N_PERM}
OUT["model_panel"] = {"summary": summ, "per_model": panel}
print("panel summary", summ)

with open(DATA / "sic_results.json", "w") as f:
    json.dump(OUT, f, indent=2)

plt.rcParams.update({"font.size": 9})
fig, ax = plt.subplots(figsize=(5.5, 4.5))
sc_ = ax.scatter(mv.anneal_T, mv.implant_T, c=mv.act, s=90, cmap="viridis", edgecolor="k",
                 linewidth=0.6, vmin=0, vmax=100)
for _, r in mv.iterrows():
    ax.annotate(f"{r.act:.0f}%", (r.anneal_T, r.implant_T), fontsize=6.5, xytext=(3, 3),
                textcoords="offset points")
ax.set_xlabel("Anneal temperature (°C)")
ax.set_ylabel("Implant temperature (°C)")
ax.set_title(f"SiC: implant-T vs. anneal-T vs. activation% (n={n})")
plt.colorbar(sc_, ax=ax).set_label("Activation (%)")
plt.tight_layout()
plt.savefig(FIG / "sic_implantT_annealT_activation.png", dpi=150)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
b = enum["betas"]
for ax, j, col, ttl in [(axes[0], 1, "#1f77b4", "Anneal-T coefficient"), (axes[1], 2, "#d62728", "Implant-T coefficient")]:
    ax.hist(b[:, j], bins=30, color=col, alpha=0.8)
    ax.axvline(0, color="k", linestyle="--")
    ax.set_title(f"{ttl}\n(exact source-level bootstrap, {len(b)} resamples)")
    ax.set_xlabel("Standardized coefficient")
plt.tight_layout()
plt.savefig(FIG / "sic_cluster_bootstrap_coefficients.png", dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(6.5, 6))
names = list(panel)
yp = np.arange(len(names))
full_imp = [panel[k]["effect"]["implant_T"] for k in names]
no_saks = [panel[k]["drop_source_effect"]["Saks2004"]["implant_T"] for k in names]
ann = [panel[k]["effect"]["anneal_T"] for k in names]
ax.barh(yp - 0.2, full_imp, height=0.35, color="#d62728", label="Implant-T effect, all sources")
ax.barh(yp + 0.2, ann, height=0.35, color="#1f77b4", alpha=0.6, label="Anneal-T effect, all sources")
ax.scatter(no_saks, yp - 0.2, color="k", marker="D", s=14, zorder=3, label="Implant-T effect, Saks et al. removed")
ax.axvline(0, color="k", linewidth=0.8)
ax.set_yticks(yp)
ax.set_yticklabels(names, fontsize=7.5)
ax.set_xlabel("Change in activation% from -1 SD to +1 SD")
ax.set_title(f"17 model specifications fitted to the same {n} rows")
ax.legend(fontsize=7, loc="lower right")
plt.tight_layout()
plt.savefig(FIG / "robustness_panel.png", dpi=150)
plt.close()
print("DONE")
