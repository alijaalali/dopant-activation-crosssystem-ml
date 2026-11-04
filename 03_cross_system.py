import itertools
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import DATA, FIG, load

df = load()
OUT = {}
cov = pd.read_csv(DATA / "system_covariates.csv")
ceil = {"SiC": float(cov[(cov.system == "SiC") & (cov.covariate == "practical_anneal_ceiling")].value.iloc[0]),
        "GaN": float(cov[(cov.system == "GaN") & (cov.covariate == "practical_anneal_ceiling_UHPA")].value.iloc[0])}
pool = df[df.system.isin(["SiC", "GaN"])].dropna(subset=["anneal_T", "act"]).reset_index(drop=True)
basis = df.groupby(["system", "value_basis"]).size().unstack(fill_value=0)
OUT["value_basis_counts"] = {s: {k: int(v) for k, v in basis.loc[s].items()} for s in basis.index}
OUT["text_basis_counts"] = {k: int(v) for k, v in df.text_basis.value_counts().items()}
hall = df[df.measurement_technique.str.contains("Hall", na=False) & ~df.measurement_technique.str.contains("not_Hall", na=False)]
OUT["hall_convention_by_text_basis"] = {
    tb: {"hall_rows": int(len(g)), "state_convention": int(g.hall_factor_convention.notna().sum())}
    for tb, g in hall.groupby("text_basis")}
print(OUT["value_basis_counts"], OUT["text_basis_counts"], OUT["hall_convention_by_text_basis"])
OUT["pooled_n"] = {s: int((pool.system == s).sum()) for s in ["SiC", "GaN"]}
OUT["pooled_sources"] = {s: int(pool[pool.system == s].source_id.nunique()) for s in ["SiC", "GaN"]}


def normalize(frame, kind, cs=None, cg=None):
    T = frame.anneal_T.to_numpy(float)
    out = np.empty(len(frame))
    for s in ["SiC", "GaN"]:
        m = (frame.system == s).to_numpy()
        if kind == "raw":
            out[m] = T[m]
        elif kind == "ceiling":
            out[m] = T[m] / {"SiC": cs or ceil["SiC"], "GaN": cg or ceil["GaN"]}[s]
        elif kind == "rank":
            out[m] = (stats.rankdata(T[m]) - 0.5) / m.sum()
        elif kind == "zscore":
            out[m] = (T[m] - T[m].mean()) / T[m].std(ddof=1)
    return out


def pooled_r(frame, kind, **kw):
    return float(stats.pearsonr(normalize(frame, kind, **kw), frame.act)[0])


rng = np.random.default_rng(3)
by_sys = {s: pool[pool.system == s] for s in ["SiC", "GaN"]}
srcs = {s: sorted(by_sys[s].source_id.unique()) for s in by_sys}


def boot_ci(kind):
    vals = []
    for _ in range(5000):
        parts = []
        for s in ["SiC", "GaN"]:
            pick = rng.choice(srcs[s], size=len(srcs[s]), replace=True)
            parts.append(pd.concat([by_sys[s][by_sys[s].source_id == p] for p in pick]))
        fr = pd.concat(parts).reset_index(drop=True)
        if fr.anneal_T.std() == 0 or fr.act.std() == 0 or any(fr[fr.system == s].anneal_T.nunique() < 2 for s in ["SiC", "GaN"]):
            continue
        vals.append(pooled_r(fr, kind))
    return [float(v) for v in np.percentile(vals, [2.5, 97.5])], len(vals)


OUT["pooled_correlation"] = {}
for kind in ["raw", "ceiling", "rank", "zscore"]:
    ci, nv = boot_ci(kind)
    r = pooled_r(pool, kind)
    OUT["pooled_correlation"][kind] = {"r": r, "source_level_bootstrap_ci95": ci, "n_valid_resamples": nv,
                                       "row_level_p": float(stats.pearsonr(normalize(pool, kind), pool.act)[1])}
    print(kind, OUT["pooled_correlation"][kind])

pool_noapprox = pool[pool.value_basis != "approximate"].reset_index(drop=True)
OUT["pooled_correlation_excluding_approximate_rows"] = {
    "n": {s: int((pool_noapprox.system == s).sum()) for s in ["SiC", "GaN"]},
    **{kind: pooled_r(pool_noapprox, kind) for kind in ["raw", "ceiling", "rank", "zscore"]}}
print("excluding approximate rows", OUT["pooled_correlation_excluding_approximate_rows"])

grid = []
for cs, cg in itertools.product([2100, 2300, 2500, 2700], [1480, 1600, 1800, 2000]):
    grid.append({"ceiling_SiC": cs, "ceiling_GaN": cg, "r": pooled_r(pool, "ceiling", cs=cs, cg=cg)})
rs = [g["r"] for g in grid]
OUT["ceiling_sweep"] = {"r_min": min(rs), "r_max": max(rs), "grid": grid}
print("ceiling sweep", min(rs), max(rs))

OUT["within_system_correlation"] = {
    s: {"r": float(stats.pearsonr(by_sys[s].anneal_T, by_sys[s].act)[0]), "n": len(by_sys[s])} for s in by_sys}
print(OUT["within_system_correlation"])

rows = []
for s in ["SiC", "GaN", "Diamond", "GaAs", "InP", "Si"]:
    sub = df[df.system == s]
    T = sub.anneal_T.dropna()
    use = sub.dropna(subset=["anneal_T", "act"])
    rows.append({"system": s, "n_rows": int(len(sub)), "n_sources": int(sub.source_id.nunique()),
                 "n_rows_anneal_and_activation": int(len(use)), "n_sources_anneal_and_activation": int(use.source_id.nunique()),
                 "anneal_T_min": float(T.min()), "anneal_T_max": float(T.max()),
                 "anneal_T_min_excluding_unverified": float(sub[sub.row_id != 2].anneal_T.min())
                 if s == "SiC" else float(T.min())})
OUT["six_system_summary"] = rows
for r in rows:
    print(r)
OUT["dataset_totals"] = {"rows": int(len(df)), "sources": int(df.source_id.nunique())}

with open(DATA / "cross_system_results.json", "w") as f:
    json.dump(OUT, f, indent=2)

order = ["Si", "GaAs", "InP", "Diamond", "GaN", "SiC"]
colors = {"SiC": "#b03a2e", "GaN": "#2471a3", "Diamond": "#7d3c98", "GaAs": "#229954", "InP": "#d68910", "Si": "#566573"}
fig, ax = plt.subplots(figsize=(6.2, 3.5))
for i, s in enumerate(order):
    r = next(x for x in rows if x["system"] == s)
    lo, hi = r["anneal_T_min_excluding_unverified"], r["anneal_T_max"]
    if lo == hi:
        lo, hi = lo - 15, hi + 15
    ax.barh(i, hi - lo, left=lo, height=0.55, color=colors[s], alpha=0.85 if s == "SiC" else 0.55,
            edgecolor="black", linewidth=1.2 if s == "SiC" else 0.6)
    lab = f"{s}  (n={r['n_rows']} rows, {r['n_sources']} sources)"
    ax.text(hi + 40, i, lab, va="center", fontsize=8.5, fontweight="bold" if s == "SiC" else "normal")
ax.set_yticks([])
ax.set_xlabel("Numeric post-implantation anneal temperatures in the compiled literature (°C)")
ax.set_xlim(0, 3300)
ax.set_title("Anneal-temperature range across six compiled dopant-implantation systems\n(bar = reported range, not activation level)", fontsize=9)
ax.spines[["top", "right", "left"]].set_visible(False)
plt.tight_layout()
plt.savefig(FIG / "sixsystem_comparison.png", dpi=200)
plt.savefig(FIG / "sixsystem_comparison.pdf")
plt.close()

fig, ax = plt.subplots(figsize=(4.7, 3.0))
labs = {"raw": "Raw", "ceiling": "Divided by\nsystem ceiling", "rank": "Within-system\nrank", "zscore": "Within-system\nz-score"}
for i, k in enumerate(labs):
    c = OUT["pooled_correlation"][k]
    ax.errorbar(i, c["r"], yerr=[[c["r"] - c["source_level_bootstrap_ci95"][0]], [c["source_level_bootstrap_ci95"][1] - c["r"]]],
                fmt="o", color="k", capsize=4)
ax.axhline(0, color="gray", linestyle="--")
ax.set_xticks(range(len(labs)))
ax.set_xticklabels(labs.values())
ax.set_ylabel("Pooled Pearson r (anneal-T vs. activation)")
ax.set_title("SiC + GaN, error bars: 95% source-level bootstrap")
plt.tight_layout()
plt.savefig(FIG / "cross_system_normalized_comparison.png", dpi=150)
plt.savefig(FIG / "cross_system_normalized_comparison.pdf")
plt.close()
print("DONE")
