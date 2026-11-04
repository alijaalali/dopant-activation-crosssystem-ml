import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from common import DATA, FIG

res = json.load(open(DATA / "sensitivity_results.json"))
flow = res["dataset_flow"]
n_labs = json.load(open(DATA / "hierarchical_results.json"))["n_labs"]
origin = {r["origin"]: r for r in flow["by_origin"]}
system = {r["system"]: r for r in flow["by_system"]}
dopant = {r["dopant"]: r for r in flow["sic_by_dopant"]}
steps = flow["sic_al_steps"]
init = flow["initial_sic_steps"]
other_rows = flow["database_rows"] - system["SiC"]["rows"]
other_src = sum(r["sources"] for s, r in system.items() if s != "SiC")


def plural(n, word):
    return f"{n} {word}" + ("" if n == 1 else "s")


fig, ax = plt.subplots(figsize=(6.5, 7.8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 12)
ax.axis("off")


def box(x, y, w, h, text, fc="#f2f2f2", fs=8.6, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fc, ec="#444444", lw=1.0))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", linespacing=1.3)


def arrow(x0, y0, x1, y1, style="-|>", ls="-"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle=style, lw=1.0, color="#444444", linestyle=ls, shrinkA=0, shrinkB=0))


MX, MW = 4.7, 5.0
RX, RW = 8.7, 2.5
box(3.6, 11.2, 4.6, 1.0, f"Initial compilation\n{origin['initial_compilation']['rows']} rows, {origin['initial_compilation']['sources']} sources")
box(7.9, 11.2, 3.6, 1.0, "13 full-text papers read\n2 gave no usable rows; 72 rows\nextracted, 2 removed, 5 corrected", fs=8.2)
arrow(3.6, 10.7, 4.2, 9.95)
arrow(7.9, 10.7, 5.8, 9.95)
box(MX, 9.5, MW, 0.9, f"Database: {flow['database_rows']} rows, {flow['database_sources']} sources\n"
    f"({origin['expansion']['rows']} rows added from full texts)", fc="#dce9f5")
arrow(MX, 9.05, MX, 8.3)
arrow(MX, 8.7, RX - RW / 2, 8.7)
box(MX, 7.85, MW, 0.9, f"4H-SiC: {system['SiC']['rows']} rows, {system['SiC']['sources']} sources")
box(RX, 8.7, RW, 1.3, f"Other systems\n{other_rows} rows, {other_src} sources\nqualitative only\n(Supplementary Material)", fs=8.0)
arrow(MX, 7.4, MX, 6.65)
box(MX, 6.2, MW, 0.9, f"Al-implanted: {steps[0][1]} rows, {steps[0][2]} sources")
arrow(MX, 7.0, RX - RW / 2, 7.0)
box(RX, 7.0, RW, 1.0, f"B-implanted\n{plural(dopant['B']['rows'], 'row')}, {plural(dopant['B']['sources'], 'source')}\nnot modelled", fs=8.0)
arrow(MX, 5.75, MX, 5.1)
removed = [st for st in steps[1:] if st[3] > 0]
lines = ["Rows removed before modelling"] + [
    f"{st[3]}: {st[0].replace('after removing rows with ', '').replace('after removing ', '')}" for st in removed]
box(MX, 4.3, MW, 1.6, "\n".join(lines), fc="#fbe9e0", fs=8.2)
arrow(MX, 3.5, MX, 2.8)
box(MX, 2.35, MW, 0.9, f"Hierarchical model and validation\n{steps[-1][1]} rows, {plural(steps[-1][2], 'source')}, {n_labs} laboratory groups", fc="#dce9f5", bold=True)

box(MX, 0.7, 7.0, 1.0, "Initial multivariate SiC model (Section 3.2)\n" + "\n".join(
    f"{st[0]}: {st[1]} rows, {st[2]} sources" for st in init[-2:]), fc="#e5f2e0", fs=8)
ax.plot([0.35, 0.35, 1.0], [11.2, 0.7, 0.7], color="#444444", lw=1.0, ls="--")
ax.plot([0.35, 1.3], [11.2, 11.2], color="#444444", lw=1.0, ls="--")
ax.text(0.5, 6.0, f"SiC rows of the initial compilation: {init[0][1]} rows, {init[0][2]} sources", rotation=90, ha="left", va="center", fontsize=8, style="italic")
arrow(1.0, 0.7, 1.2, 0.7)
plt.tight_layout()
plt.savefig(FIG / "dataset_flowchart.png", dpi=200)
plt.savefig(FIG / "dataset_flowchart.pdf")
print("saved")
