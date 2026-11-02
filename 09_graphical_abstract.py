import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import DATA, FIG

h = json.load(open(DATA / "hierarchical_results.json"))
k = json.load(open(DATA / "kinetic_validation_results.json"))
share = h["main"]["share"]["between_groups_total"]
eff = h["main"]["marginal_effect_pp_per_sd"]["implant_temp_C"]
eff_drop = h["sensitivity"]["drop_Saks2004"]["marginal_pp"]
schemes = [("leave_one_source_out", "Held-out\nsources (6)"), ("leave_one_lab_out", "Held-out\nlaboratories (5)"),
           ("temporal", "Later\npublications (2)")]
models = [("training_mean", "Training mean", "#1b9e77"), ("linear_mixed", "Linear mixed", "#7570b3"),
          ("linear_mixed_without_implantT", "Linear, no implant T", "#a6a1d4"),
          ("arrhenius_kinetic_mixed", "Arrhenius time", "#e7298a"), ("random_forest", "Random forest", "#d95f02")]

plt.rcParams.update({"font.size": 11})
fig, axes = plt.subplots(1, 3, figsize=(13.5, 5), gridspec_kw={"width_ratios": [0.75, 1.0, 1.7]})

ax = axes[0]
s = 100 * share["mean"]
ax.bar([0], [s], color="#d95f02", width=0.75)
ax.bar([0], [100 - s], bottom=[s], color="#bdbdbd", width=0.75)
ax.set_xlim(-0.7, 0.7)
ax.set_ylim(0, 100)
ax.set_xticks([])
ax.set_ylabel("Share of variance in logit activation (%)")
ax.text(0, s / 2, f"Between sources\nand laboratories\n{s:.0f}%\n(95% interval\n{100 * share['ci95'][0]:.0f}–{100 * share['ci95'][1]:.0f}%)", ha="center", va="center", color="w", fontweight="bold", fontsize=10)
ax.text(0, s + (100 - s) / 2, f"Within\nsources\n{100 - s:.0f}%", ha="center", va="center")
ax.set_title("Between-study\nvariation", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
rows = [("Without\nSaks 2004", eff_drop), ("All 6 sources", eff)]
for i, (lab, e) in enumerate(rows):
    ax.errorbar(e["mean"], i, xerr=[[e["mean"] - e["ci95"][0]], [e["ci95"][1] - e["mean"]]], fmt="o", color="#1b9e77", capsize=4, markersize=6)
    ax.text(e["ci95"][1] + 2, i, f"{e['mean']:+.0f}", va="center")
ax.axvline(0, color="k", linestyle="--", linewidth=0.8)
ax.set_yticks([0, 1])
ax.set_yticklabels([r[0] for r in rows])
ax.set_ylim(-0.6, 1.6)
ax.set_xlim(-20, 70)
ax.set_xlabel("Association with implant temperature\n(points per 281 °C, 95% interval)")
ax.set_title("Association rests on\none study", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[2]
w = 0.16
for j, (key, label, col) in enumerate(models):
    vals = [k[sch][key]["mae_pp"] for sch, _ in schemes]
    xs = np.arange(len(schemes)) + (j - 2) * w
    ax.bar(xs, vals, width=w * 0.92, color=col, label=label)
    for x, v, (sch, _) in zip(xs, vals, schemes):
        ax.text(x, v + 0.8, f"{v:.1f}", ha="center", fontsize=7.5, rotation=90, va="bottom")
ax.set_xticks(range(len(schemes)))
ax.set_xticklabels([lab for _, lab in schemes], fontsize=10)
ax.set_ylabel("Error on held-out data (points of activation)")
ax.set_ylim(0, 78)
ax.legend(fontsize=8, ncol=2, loc="upper left", frameon=False)
ax.set_title("No consistent gain over\nthe training mean", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(FIG / "graphical_abstract.png", dpi=150)
plt.savefig(FIG / "graphical_abstract.pdf")
print("saved")
