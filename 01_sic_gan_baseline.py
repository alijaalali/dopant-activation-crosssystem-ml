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
cov = pd.read_csv(DATA / "system_covariates.csv", comment="#")
df["anneal_num"] = pd.to_numeric(df["anneal_temp_C"], errors="coerce")
df["implant_num"] = pd.to_numeric(df["implant_temp_C"], errors="coerce")

s1s2_mask = df["source_id"] == "S1_S2"
FILLED_IMPLANT_T_FLAG = s1s2_mask & df["implant_num"].isna()
df.loc[FILLED_IMPLANT_T_FLAG, "implant_num"] = 500.0
n_filled = int(FILLED_IMPLANT_T_FLAG.sum())

print("="*78)
print("DATA ASSEMBLY")
print("="*78)
print(f"S1_S2 rows with implant_temp_C set to 500 C (stated in the source methods): {n_filled}")

usable = df[df["activation_pct"].notna() & df["anneal_num"].notna()].copy()
sic_all = usable[usable.system == "SiC"].copy()
gan_all = usable[usable.system == "GaN"].copy()
print(f"\nUsable rows (anneal-T and activation%): n={len(usable)} "
      f"(SiC={len(sic_all)}, GaN={len(gan_all)})")

sic_mv = sic_all[sic_all.implant_num.notna()].copy()
missing_implant_sources = sorted(sic_all[sic_all.implant_num.isna()].source_id.unique().tolist())
print(f"\nSiC rows with BOTH anneal-T and implant-T known: n={len(sic_mv)} "
      f"across sources {sorted(sic_mv.source_id.unique().tolist())}")
print(f"SiC rows excluded from the multivariate model (implant-T not reported): "
      f"n={len(sic_all)-len(sic_mv)} from sources {missing_implant_sources}")

OUT["n_total_usable"] = {"SiC": len(sic_all), "GaN": len(gan_all)}
OUT["sic_multivariate_n"] = len(sic_mv)
OUT["sic_multivariate_sources"] = sorted(sic_mv.source_id.unique().tolist())
OUT["sic_excluded_sources_no_implantT"] = missing_implant_sources
OUT["implant_T_filled_rows"] = n_filled

row6_flag = sic_mv[(sic_mv.source_id=="S3") & (sic_mv.orientation=="11-20")]
if len(row6_flag):
    print(f"\nS3 (11-20) orientation row (row_id={int(row6_flag.row_id.iloc[0])}): orientation, "
          f"implant temperature and carbon co-implantation all differ from the other S3 rows")

print("\n" + "="*78)
print("1. PRIMARY SiC MULTIVARIATE MODEL: activation_pct ~ anneal_T + implant_T")
print("="*78)

X_raw = sic_mv[["anneal_num", "implant_num"]].to_numpy(dtype=float)
y = sic_mv["activation_pct"].to_numpy(dtype=float)
n = len(y)

X_mean, X_std = X_raw.mean(axis=0), X_raw.std(axis=0, ddof=1)
Xz = (X_raw - X_mean) / X_std
Xd = np.column_stack([np.ones(n), Xz])

def ols_fit(Xd, y):
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    resid = y - Xd @ beta
    return beta, resid

beta, resid = ols_fit(Xd, y)
ss_res = np.sum(resid**2)
ss_tot = np.sum((y - y.mean())**2)
r2_insample = 1 - ss_res/ss_tot
print(f"n={n} (4 independent sources: {OUT['sic_multivariate_sources']})")
print(f"In-sample R^2 = {r2_insample:.3f}")
print(f"Standardized coefficients: intercept={beta[0]:.2f}, "
      f"anneal_T_z={beta[1]:.2f}, implant_T_z={beta[2]:.2f}")

n_boot = 5000
boot_betas = []
idx = np.arange(n)
for _ in range(n_boot):
    samp = rng.choice(idx, size=n, replace=True)
    Xb, yb = Xd[samp], y[samp]
    if np.linalg.matrix_rank(Xb) < 3:
        continue
    b, _ = ols_fit(Xb, yb)
    boot_betas.append(b)
boot_betas = np.array(boot_betas)
ci_lo, ci_hi = np.percentile(boot_betas, [2.5, 97.5], axis=0)
print(f"\nBootstrap 95% CI ({len(boot_betas)}/{n_boot} valid resamples):")
for i, name in enumerate(["intercept", "anneal_T (standardized)", "implant_T (standardized)"]):
    print(f"  {name}: {beta[i]:.2f}  CI [{ci_lo[i]:.2f}, {ci_hi[i]:.2f}]"
          f"  {'(excludes 0)' if ci_lo[i]*ci_hi[i]>0 else '(includes 0)'}")

OUT["sic_multivariate_model"] = {
    "n": n, "r2_insample": float(r2_insample),
    "coef_intercept": float(beta[0]), "coef_anneal_T_std": float(beta[1]),
    "coef_implant_T_std": float(beta[2]),
    "ci_intercept": [float(ci_lo[0]), float(ci_hi[0])],
    "ci_anneal_T_std": [float(ci_lo[1]), float(ci_hi[1])],
    "ci_implant_T_std": [float(ci_lo[2]), float(ci_hi[2])],
    "X_mean": X_mean.tolist(), "X_std": X_std.tolist(),
}

print("\n" + "="*78)
print("2. LEAVE-ONE-SOURCE-OUT CV")
print("="*78)
print("Leave-one-source-out CV: 4 folds for the multivariate model; 5 folds for the "
      "anneal-T-only model on the full SiC set.")

def loso_multivariate(data):
    srcs = sorted(data.source_id.unique())
    preds, actuals, fold_srcs = [], [], []
    for held in srcs:
        train = data[data.source_id != held]
        test = data[data.source_id == held]
        Xtr = train[["anneal_num","implant_num"]].to_numpy(float)
        ytr = train["activation_pct"].to_numpy(float)
        Xte = test[["anneal_num","implant_num"]].to_numpy(float)
        yte = test["activation_pct"].to_numpy(float)
        mu, sd = Xtr.mean(0), Xtr.std(0, ddof=1)
        Xtrz = np.column_stack([np.ones(len(Xtr)), (Xtr-mu)/sd])
        Xtez = np.column_stack([np.ones(len(Xte)), (Xte-mu)/sd])
        b, _ = ols_fit(Xtrz, ytr)
        yhat = Xtez @ b
        preds.extend(yhat.tolist()); actuals.extend(yte.tolist())
        fold_srcs.extend([held]*len(yte))
        mae_fold = np.mean(np.abs(yhat-yte))
        print(f"  held-out source={held:12s} n_test={len(yte)}  fold MAE={mae_fold:6.1f} pct-points")
    preds, actuals = np.array(preds), np.array(actuals)
    mae = np.mean(np.abs(preds-actuals))
    rmse = np.sqrt(np.mean((preds-actuals)**2))
    r_cv, p_cv = stats.pearsonr(preds, actuals) if np.std(preds)>0 else (np.nan,np.nan)
    return mae, rmse, r_cv, p_cv, preds, actuals, fold_srcs

print("\n-- 4-fold LOSO, multivariate model (n=%d) --" % n)
mae_mv, rmse_mv, r_mv, p_mv, preds_mv, acts_mv, folds_mv = loso_multivariate(sic_mv)
print(f"Overall LOSO: MAE={mae_mv:.1f} pct-points, RMSE={rmse_mv:.1f}, "
      f"pred-vs-actual r={r_mv:.3f} (p={p_mv:.3f})")
print("Only 4 source-level folds are available; the estimate is noisy.")

def loso_univariate(data):
    srcs = sorted(data.source_id.unique())
    preds, actuals = [], []
    for held in srcs:
        train = data[data.source_id != held]
        test = data[data.source_id == held]
        xtr, ytr = train["anneal_num"].to_numpy(float), train["activation_pct"].to_numpy(float)
        xte, yte = test["anneal_num"].to_numpy(float), test["activation_pct"].to_numpy(float)
        slope, intercept, *_ = stats.linregress(xtr, ytr)
        yhat = slope*xte + intercept
        preds.extend(yhat.tolist()); actuals.extend(yte.tolist())
    preds, actuals = np.array(preds), np.array(actuals)
    mae = np.mean(np.abs(preds-actuals))
    return mae, preds, actuals

print("\n-- 5-fold LOSO, univariate anneal-T-only model, full SiC set --")
mae_uv, preds_uv, acts_uv = loso_univariate(sic_all)
print(f"Overall LOSO MAE (anneal-T-only, n={len(sic_all)}, 5 sources): {mae_uv:.1f} pct-points")
print(f"Comparison: multivariate (anneal-T+implant-T) LOSO MAE={mae_mv:.1f} vs "
      f"univariate (anneal-T-only) LOSO MAE={mae_uv:.1f} -- "
      f"{'multivariate model has lower error' if mae_mv < mae_uv else 'univariate model has lower or comparable error'}")

OUT["loso_cv"] = {
    "multivariate_4fold": {"n_sources": 4, "mae": float(mae_mv), "rmse": float(rmse_mv),
                            "pred_vs_actual_r": float(r_mv) if not np.isnan(r_mv) else None},
    "univariate_5fold_anneal_only": {"n_sources": 5, "n_rows": len(sic_all), "mae": float(mae_uv)},
}

print("\n" + "="*78)
print("3. GaN rows and pooled cross-system correlation")
print("="*78)
print(gan_all[["source_id","anneal_num","activation_pct"]].to_string(index=False))
gan_context = df[(df.system=="GaN") & df.activation_pct.notna() & df.anneal_num.isna()]
print(f"\nGaN rows with activation% but no numeric anneal-T (n={len(gan_context)}): "
      f"{gan_context.source_id.tolist()}")
print(gan_context[["source_id","activation_pct","activation_pct_note"]].to_string(index=False))

usable = usable[usable.system.isin(["SiC", "GaN"])].copy()

ceilings = {}
for sysname, ceiling_cov in [("SiC","practical_anneal_ceiling"), ("GaN","practical_anneal_ceiling_UHPA")]:
    row = cov[(cov.system==sysname)&(cov.covariate==ceiling_cov)]
    ceilings[sysname] = float(row["value"].iloc[0])
usable["anneal_frac_of_ceiling"] = usable.apply(lambda r: r.anneal_num/ceilings[r.system], axis=1)

r_raw, p_raw = stats.pearsonr(usable.anneal_num, usable.activation_pct)
r_norm, p_norm = stats.pearsonr(usable.anneal_frac_of_ceiling, usable.activation_pct)

def boot_r(x, y, n_boot=5000):
    idx = np.arange(len(x)); boots=[]
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        xs, ys = x.iloc[s] if hasattr(x,'iloc') else x[s], y.iloc[s] if hasattr(y,'iloc') else y[s]
        if np.std(xs)>0:
            boots.append(stats.pearsonr(xs,ys)[0])
    return np.percentile(boots, [2.5,97.5])

ci_raw = boot_r(usable.anneal_num, usable.activation_pct)
ci_norm = boot_r(usable.anneal_frac_of_ceiling, usable.activation_pct)
print(f"\nPooled cross-system (n={len(usable)}): raw-T r={r_raw:.3f} (p={p_raw:.3f}) "
      f"bootstrap CI [{ci_raw[0]:.3f},{ci_raw[1]:.3f}]")
print(f"Pooled cross-system normalized-T r={r_norm:.3f} (p={p_norm:.3f}) "
      f"bootstrap CI [{ci_norm[0]:.3f},{ci_norm[1]:.3f}]")

OUT["gan_contrast"] = {
    "n_fitted": len(gan_all),
    "rows": gan_all[["source_id","anneal_num","activation_pct"]].to_dict("records"),
    "n_qualitative_context_rows": len(gan_context),
    "pooled_raw_r": float(r_raw), "pooled_raw_p": float(p_raw), "pooled_raw_ci": [float(ci_raw[0]),float(ci_raw[1])],
    "pooled_norm_r": float(r_norm), "pooled_norm_p": float(p_norm), "pooled_norm_ci": [float(ci_norm[0]),float(ci_norm[1])],
}

print("\n" + "="*78)
print("4. HALL-FACTOR-CONVENTION TRANSPARENCY GAP + BOUNDING CALCULATION")
print("="*78)
stated = df[df["hall_factor_convention"].notna() & ~df["hall_factor_convention"].astype(str).str.contains("not_Hall", na=False)]
hall_rows = df[~df["measurement_technique"].astype(str).str.contains("not_Hall", na=False)]
print(f"Rows using Hall-effect measurement: {len(hall_rows)} / {len(df)} total rows.")
print(f"Of those, rows EXPLICITLY stating their Hall-scattering-factor (r_H) "
      f"convention: {len(stated)} / {len(hall_rows)} "
      f"({100*len(stated)/len(hall_rows):.0f}%).")

r_H_naive = 1.0
r_H_lit_range = (0.2, 0.5)
overstatement_factor_range = (r_H_naive/r_H_lit_range[1], r_H_naive/r_H_lit_range[0])
print(f"\nSensitivity bound: r_H in {r_H_lit_range[0]}-{r_H_lit_range[1]} versus the r_H=1 assumption.")
print("  R_H = r_H/(p*e), so p_apparent = 1/(e*R_H) and p_true = r_H * p_apparent.")
print(f"  Activation computed under r_H=1 could exceed the true value by a factor of "
      f"{overstatement_factor_range[0]:.1f} to {overstatement_factor_range[1]:.1f}.")
print("  The bound is illustrative and is not applied to any individual row.")

OUT["hall_factor_gap"] = {
    "n_hall_rows": len(hall_rows), "n_stating_convention": len(stated),
    "fraction_stating": len(stated)/len(hall_rows),
    "r_H_naive_assumption": r_H_naive, "r_H_literature_range": list(r_H_lit_range),
    "overstatement_factor_range": [float(overstatement_factor_range[0]), float(overstatement_factor_range[1])],
    "note": "Illustrative sensitivity bound from literature r_H values for p-type 4H-SiC; not applied as a per-row correction."
}

print("\n" + "="*78)
print("5. EFFECT SIGNS AND CONFIDENCE INTERVALS")
print("="*78)
sign_anneal = "positive" if beta[1] > 0 else "negative"
sign_implant = "positive" if beta[2] > 0 else "negative"
ci_anneal_real = ci_lo[1]*ci_hi[1] > 0
ci_implant_real = ci_lo[2]*ci_hi[2] > 0
print(f"Standardized effect of anneal temperature: {sign_anneal}, bootstrap CI [{ci_lo[1]:.2f},{ci_hi[1]:.2f}] "
      f"({'excludes' if ci_anneal_real else 'includes'} 0)")
print(f"Standardized effect of implant temperature: {sign_implant}, bootstrap CI [{ci_lo[2]:.2f},{ci_hi[2]:.2f}] "
      f"({'excludes' if ci_implant_real else 'includes'} 0)")
OUT["effect_signs"] = {"anneal_T": sign_anneal, "implant_T": sign_implant,
                        "anneal_T_CI_excludes_zero": bool(ci_anneal_real),
                        "implant_T_CI_excludes_zero": bool(ci_implant_real)}

fig, ax = plt.subplots(figsize=(6,5))
sc = ax.scatter(sic_mv.anneal_num, sic_mv.implant_num, c=sic_mv.activation_pct,
                 s=120, cmap="viridis", edgecolor="k")
for _, r in sic_mv.iterrows():
    ax.annotate(f"{r.activation_pct:.0f}%", (r.anneal_num, r.implant_num),
                textcoords="offset points", xytext=(6,4), fontsize=8)
ax.set_xlabel("Anneal temperature (C)")
ax.set_ylabel("Implant temperature (C)")
ax.set_title(f"SiC: activation pct vs anneal-T and implant-T (n={n}, 4 sources)")
plt.colorbar(sc, label="Activation %")
plt.tight_layout()
plt.savefig(FIG / "sic_implantT_annealT_activation.png", dpi=150)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(10,4))
axes[0].hist(boot_betas[:,1], bins=40, color="tab:blue", alpha=0.7)
axes[0].axvline(0, color="k", linestyle="--")
axes[0].set_title("Bootstrap: anneal_T coefficient (std.)")
axes[1].hist(boot_betas[:,2], bins=40, color="tab:orange", alpha=0.7)
axes[1].axvline(0, color="k", linestyle="--")
axes[1].set_title("Bootstrap: implant_T coefficient (std.)")
plt.tight_layout()
plt.savefig(FIG / "sic_bootstrap_coefficients.png", dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(6,5))
ax.scatter(sic_all.anneal_frac_of_ceiling if "anneal_frac_of_ceiling" in sic_all else sic_all.anneal_num/ceilings["SiC"],
           sic_all.activation_pct, label=f"SiC (n={len(sic_all)})", color="tab:blue", s=70)
ax.scatter(gan_all.anneal_num/ceilings["GaN"], gan_all.activation_pct,
           label=f"GaN (n={len(gan_all)})", color="tab:red", s=70, marker="^")
ax.set_xlabel("Anneal temperature / system practical ceiling")
ax.set_ylabel("Activation %")
ax.set_title("Cross-system comparison, normalized anneal-T")
ax.legend()
plt.tight_layout()
plt.savefig(FIG / "cross_system_normalized_comparison.png", dpi=150)
plt.close()

print("\nFigures saved.")

with open(DATA / "sic_gan_baseline_results.json", "w") as f:
    json.dump(OUT, f, indent=2)
print("\nSaved data/sic_gan_baseline_results.json")
print("\n" + "="*78)
print("DONE")
print("="*78)
