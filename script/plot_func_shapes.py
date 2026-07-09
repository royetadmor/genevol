import sys
import csv
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 2:
    print("Usage: plot_func_shapes.py <results.csv> [output_dir]")
    print("  CSV columns: dataset,rate_type,func,params,is_winner")
    print("  params field uses semicolons to separate values (e.g. 2.82621;0.102448)")
    sys.exit(1)

csv_path   = sys.argv[1]
output_dir = sys.argv[2] if len(sys.argv) > 2 else "figures"

# ── Function implementations (mirroring GeneCountDependencyFunction.cpp) ────
def f_const(x, p):       return np.full_like(x, p[0])
def f_linear_bd(x, p):   return p[0] * x
def f_linear(x, p):      return p[0] + p[1] * x
def f_exp(x, p):         return p[0] * np.exp(p[1] * x)
def f_polynomial(x, p):  return p[0] * np.power(x + p[1], p[2])

FUNC_IMPLS = {
    "CONST":      f_const,
    "LINEAR_BD":  f_linear_bd,
    "LINEAR":     f_linear,
    "EXP":        f_exp,
    "POLYNOMIAL": f_polynomial,
}

FUNC_COLORS = {
    "CONST":      "#4C72B0",
    "LINEAR_BD":  "#DD8452",
    "LINEAR":     "#55A868",
    "EXP":        "#C44E52",
    "POLYNOMIAL": "#8172B2",
}

# ── Load CSV ─────────────────────────────────────────────────────────────────
# data[rate_type][dataset][func] = params list
# winners[rate_type][dataset] = func name
data    = {}
winners = {}

with open(csv_path, newline="") as f:
    for row in csv.DictReader(f):
        ds       = row["dataset"]
        rt       = row["rate_type"]
        func     = row["func"]
        params   = [float(v) for v in row["params"].split(";")]
        is_win   = row["is_winner"].strip().lower() == "true"

        data.setdefault(rt, {}).setdefault(ds, {})[func] = params
        if is_win:
            winners.setdefault(rt, {})[ds] = func

X = np.linspace(1, 20, 300)


def plot_panel(ax, dataset, funcs, winner):
    winner_vals = FUNC_IMPLS[winner](X, funcs[winner])
    winner_vals = np.nan_to_num(winner_vals, nan=0, posinf=0, neginf=0)
    y_cap = max(np.max(winner_vals) * 1.5, 1.0)

    for fname, params in funcs.items():
        vals = FUNC_IMPLS[fname](X, params)
        vals = np.nan_to_num(vals, nan=0, posinf=y_cap, neginf=0)
        is_winner = fname == winner
        ax.plot(X, vals,
                color=FUNC_COLORS[fname],
                lw=2.5 if is_winner else 1.0,
                alpha=1.0 if is_winner else 0.35,
                zorder=3 if is_winner else 2,
                label=fname)

    ax.set_ylim(0, y_cap)
    ax.set_title(dataset, fontsize=11)
    ax.set_xlabel("Copy number", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_figure(rt_data, rt_winners, ylabel, outpath):
    datasets = list(rt_data.keys())
    n = len(datasets)
    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, dataset in zip(axes, datasets):
        plot_panel(ax, dataset, rt_data[dataset], rt_winners[dataset])

    axes[0].set_ylabel(ylabel, fontsize=11)

    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(FUNC_IMPLS),
               fontsize=10, bbox_to_anchor=(0.5, 1.05))

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Saved {outpath}")


for rate_type, rt_data in data.items():
    outpath = f"{output_dir}/{rate_type}_func_shapes.png"
    make_figure(rt_data, winners[rate_type], f"{rate_type.capitalize()} rate", outpath)
