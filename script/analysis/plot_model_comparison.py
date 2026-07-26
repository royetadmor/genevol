import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

if len(sys.argv) < 2:
    print("Usage: plot_model_comparison.py <results.csv> [output_dir]")
    print("  CSV columns: dataset,model,AIC")
    sys.exit(1)

csv_path   = sys.argv[1]
output_dir = sys.argv[2] if len(sys.argv) > 2 else "figures"

# ── Load CSV ─────────────────────────────────────────────────────────────────
# aic[dataset][model] = AIC value
aic = {}
with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        aic.setdefault(row["dataset"], {})[row["model"]] = float(row["AIC"])

datasets = list(aic.keys())
models   = sorted({m for ds in aic.values() for m in ds})

COLORS = ["#8172B2", "#4C72B0", "#DD8452", "#55A868", "#C44E52"]
model_colors = {m: COLORS[i % len(COLORS)] for i, m in enumerate(models)}

SMALL_THRESHOLD = 10_000

def make_group(small):
    rows = [(ds, aic[ds]) for ds in datasets
            if (min(aic[ds].values()) < SMALL_THRESHOLD) == small]
    rows.sort(key=lambda r: r[1].get(models[-1], float("inf")), reverse=True)
    return rows

def draw_panel(ax, rows, title):
    y     = np.arange(len(rows))
    bar_h = 0.8 / len(models)

    for i, model in enumerate(models):
        offset = (i - len(models) / 2 + 0.5) * bar_h
        values = [row[1].get(model, float("nan")) for row in rows]
        ax.barh(y + offset, values, height=bar_h, color=model_colors[model])

    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=11)
    ax.set_xlabel("AIC", fontsize=11)
    ax.set_title(title, fontsize=12, pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)

small_rows = make_group(small=True)
large_rows = make_group(small=False)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
if small_rows:
    draw_panel(axes[0], small_rows, "Smaller datasets")
if large_rows:
    draw_panel(axes[1], large_rows, "Larger datasets")

legend_handles = [mpatches.Patch(color=model_colors[m], label=m) for m in models]
fig.legend(handles=legend_handles, fontsize=10, loc="upper center",
           ncol=len(models), bbox_to_anchor=(0.5, 1.06))

plt.tight_layout()
outpath = f"{output_dir}/model_comparison.png"
plt.savefig(outpath, dpi=150, bbox_inches="tight")
print(f"Saved {outpath}")
