import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


d = json.load(open(DATA / "six_system_summary.json"))
rows = d["six_system_summary"]

order = ["Si", "GaAs", "InP", "Diamond", "GaN", "SiC"]
rows = sorted(rows, key=lambda r: order.index(r["system"]))

fig, ax = plt.subplots(figsize=(9.0, 5.0))
colors = {"SiC": "#b03a2e", "GaN": "#2471a3", "Diamond": "#7d3c98",
          "GaAs": "#229954", "InP": "#d68910", "Si": "#566573"}

for i, r in enumerate(rows):
    sysname = r["system"]
    lo, hi = r["anneal_T_min"], r["anneal_T_max"]
    if lo == hi:
        lo -= 15
        hi += 15
    ax.barh(i, hi - lo, left=lo, height=0.55, color=colors[sysname],
            alpha=0.85 if sysname == "SiC" else 0.55,
            edgecolor="black", linewidth=1.2 if sysname == "SiC" else 0.6)
    label = f"{sysname}  (n={r['n_rows']} rows, {r['n_sources']} sources"
    label += ", fitted model)" if r["fittable_model"] else ", qualitative/contrast only)"
    ax.text(hi + 40, i, label, va="center", fontsize=9.5,
            fontweight="bold" if sysname == "SiC" else "normal")

ax.set_yticks([])
ax.set_xlabel("Post-implantation anneal temperature range reported in the compiled literature (°C)")
ax.set_xlim(0, 2500)
ax.set_title("Anneal-temperature regime across six compiled dopant-implantation\nsystems (bar = reported range, not activation level)")
ax.spines[["top", "right", "left"]].set_visible(False)
plt.tight_layout()
plt.savefig(FIG / "sixsystem_comparison.png", dpi=200)
plt.savefig(FIG / "sixsystem_comparison.pdf")
print("Saved sixsystem_comparison.png/.pdf")
