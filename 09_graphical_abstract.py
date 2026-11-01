import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import DATA, FIG

h = json.load(open(DATA / "hierarchical_results.json"))
k = json.load(open(DATA / "kinetic_validation_results.json"))
m = h["main"]
share_source = m["sd_source"] ** 2 / (m["sd_source"] ** 2 + m["sd_residual"] ** 2)
eff = h["main"]["marginal_effect_pp_per_sd"]["implant_temp_C"]
eff_drop = h["sensitivity"]["drop_Saks2004"]["marginal_pp"]
loso = k["leave_one_source_out"]
names = {"training_mean": "Mean", "linear_mixed": "Linear", "linear_mixed_without_implantT": "Linear\nno impl. T",
         "arrhenius_kinetic_mixed": "Kinetic", "random_forest": "Forest"}

plt.rcParams.update({"font.size": 8.5})
fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.0), gridspec_kw={"width_ratios": [0.8, 1.1, 1.3]})

ax = axes[0]
ax.bar([0], [100 * share_source], color="#d95f02", width=0.75, label="Between sources")
ax.bar([0], [100 * (1 - share_source)], bottom=[100 * share_source], color="#bdbdbd", width=0.75, label="Within sources")
ax.set_xlim(-0.7, 0.7)
ax.set_ylim(0, 100)
ax.set_xticks([])
ax.set_ylabel("Share of variance in\nlogit activation (%)")
ax.text(0, 100 * share_source / 2, f"Between\nsources\n{100 * share_source:.0f}%", ha="center", va="center", color="w", fontweight="bold")
ax.text(0, 100 * share_source + 100 * (1 - share_source) / 2, f"Within\nsources\n{100 * (1 - share_source):.0f}%", ha="center", va="center")
ax.set_title("Between-study\nvariation", fontsize=9)
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
ax.set_xlabel("Implant-temperature effect\n(points per 281 °C, 95% interval)")
ax.set_title("Effect rests on\none study", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[2]
vals = [loso[key]["mae_pp"] for key in names]
cols = ["#1b9e77"] + ["#7570b3"] * (len(vals) - 1)
ax.bar(range(len(vals)), vals, color=cols, width=0.65)
for i, v in enumerate(vals):
    ax.text(i, v + 0.6, f"{v:.0f}", ha="center")
ax.set_xticks(range(len(vals)))
ax.set_xticklabels(list(names.values()), fontsize=7.5)
ax.set_ylabel("Error on held-out sources\n(points of activation)")
ax.set_ylim(0, 38)
ax.set_title("No model beats the mean", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(FIG / "graphical_abstract.png", dpi=200)
plt.savefig(FIG / "graphical_abstract.pdf")
print("saved")
