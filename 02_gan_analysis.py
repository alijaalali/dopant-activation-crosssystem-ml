import json

import numpy as np

from common import DATA, cluster_enumeration, load, ols, r_squared, summarize_enumeration, zscore

df = load()
gan = df[df.system == "GaN"]
OUT = {}


def run(frame, cols, label):
    d = frame.dropna(subset=cols + ["act"]).reset_index(drop=True)
    X = d[cols].to_numpy(float)
    y = d.act.to_numpy(float)
    src = d.source_id.to_numpy()
    Xz = zscore(X)
    b = ols(Xz, y)
    r2, r2adj = r_squared(Xz, y)
    enum = cluster_enumeration(src, Xz, y)
    res = {"n": len(d), "n_sources": int(len(set(src))),
           "rows_per_source": {s: int((src == s).sum()) for s in sorted(set(src))},
           "coef": {c: float(b[j + 1]) for j, c in enumerate(cols)}, "intercept": float(b[0]),
           "r2": float(r2), "r2_adj": float(r2adj),
           "bootstrap": summarize_enumeration(enum, cols)}
    drop = {}
    for s in sorted(set(src)):
        m = src != s
        A = np.column_stack([np.ones(m.sum()), Xz[m]])
        bb = np.linalg.lstsq(A, y[m], rcond=None)[0]
        drop[s] = {"n_rows_remaining": int(m.sum()), "rank": int(np.linalg.matrix_rank(A)),
                   **{c: float(bb[j + 1]) for j, c in enumerate(cols)}}
    res["drop_source"] = drop
    print(label, json.dumps({k: v for k, v in res.items() if k != "bootstrap"}))
    print("   bootstrap", json.dumps(res["bootstrap"]))
    return res


OUT["anneal_T_only"] = run(gan.rename(columns={}), ["anneal_T"], "GaN anneal-T only")
OUT["anneal_T_only_without_Rocco2020"] = run(gan[gan.source_id != "Rocco2020_JAP"], ["anneal_T"],
                                              "GaN anneal-T only, Rocco removed")
OUT["anneal_T_and_time"] = run(gan, ["anneal_T", "time"], "GaN anneal-T + time")

br = gan[gan.source_id == "Breckenridge2021"].dropna(subset=["anneal_T", "time", "act"])
Xb = zscore(br[["anneal_T", "time"]].to_numpy(float))
bb = ols(Xb, br.act.to_numpy(float))
OUT["breckenridge_only_anneal_T_and_time"] = {"n": len(br), "anneal_T": float(bb[1]), "time": float(bb[2])}
print("Breckenridge only", OUT["breckenridge_only_anneal_T_and_time"])

d7 = gan.dropna(subset=["anneal_T", "time", "act"])
OUT["time_variation_by_source"] = {s: sorted(set(g.time)) for s, g in d7.groupby("source_id")}
print(OUT["time_variation_by_source"])

with open(DATA / "gan_results.json", "w") as f:
    json.dump(OUT, f, indent=2)
