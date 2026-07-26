import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 2:
    print("Usage: plot_func_comparison.py <results.csv> [output_dir]")
    print("  CSV columns: dataset,rate_type,func,AIC")
    sys.exit(1)

csv_path   = sys.argv[1]
output_dir = sys.argv[2] if len(sys.argv) > 2 else "figures"

FUNCS = ["CONST", "LINEAR_BD", "LINEAR", "EXP", "POLYNOMIAL"]

# ── Load CSV ─────────────────────────────────────────────────────────────────
# aic[rate_type][dataset][func] = AIC value
aic = {}
with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        aic.setdefault(row["rate_type"], {}).setdefault(row["dataset"], {})[row["func"]] = float(row["AIC"])

COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2"]


def make_plot(rt_data, rate_type, outpath):
    datasets = list(rt_data.keys())
    n        = len(datasets)
    funcs    = FUNCS

    values  = np.array([[rt_data[ds].get(f, float("nan")) for f in funcs] for ds in datasets])
    best    = np.nanmin(values, axis=1, keepdims=True)
    delta   = values - best

    fig, ax = plt.subplots(figsize=(9, 5))
    x       = np.arange(len(funcs))
    offsets = np.linspace(-0.3, 0.3, n)

    for i, (dataset, color) in enumerate(zip(datasets, COLORS)):
        xpos   = x + offsets[i]
        yvals  = delta[i]
        winners = yvals == 0

        for xp, yv in zip(xpos, yvals):
            ax.plot([xp, xp], [0, yv], color=color, lw=1.5, alpha=0.7)

        ax.scatter(xpos[~winners], yvals[~winners], color=color, s=60, zorder=3)
        ax.scatter(xpos[winners],  yvals[winners] - max(delta.max() * 0.03, 1),
                   color=color, s=140, zorder=4, marker="*")

    for dataset, color in zip(datasets, COLORS):
        ax.scatter([], [], color=color, s=60, label=dataset)
    ax.scatter([], [], color="gray", s=140, marker="*", label="Winner")

    ax.set_xticks(x)
    ax.set_xticklabels(funcs, fontsize=11)
    ax.set_ylabel("ΔAIC from best", fontsize=11)
    ax.set_title(f"{rate_type.capitalize()} function comparison", fontsize=12)
    ax.legend(fontsize=10, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved {outpath}")


for rate_type, rt_data in aic.items():
    make_plot(rt_data, rate_type, f"{output_dir}/{rate_type}_func_comparison.png")
