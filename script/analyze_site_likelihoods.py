"""
Parse log.out and analyse per-site log-likelihood changes across WGD iterations.

For each accepted WGD the script shows which sites improved, which worsened,
and by how much — helping verify that each WGD is driven by genuine signal
rather than noise on a handful of outlier sites.
"""

import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

LOG_FILE = "log.out"


# ---------------------------------------------------------------------------
# 1. Parse
# ---------------------------------------------------------------------------

def parse_log(path):
    """Return a dict {label: np.array of per-site log-likelihoods}."""
    blocks = {}
    current_label = None
    current_sites = []

    site_re = re.compile(r"^\s+site\s+(\d+):\s+(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
    header_re = re.compile(r"^(Baseline|After WGD #\d+) \(per-site log-likelihoods\):")

    with open(path) as fh:
        for line in fh:
            m = header_re.match(line)
            if m:
                if current_label is not None:
                    blocks[current_label] = np.array(current_sites)
                current_label = m.group(1)
                current_sites = []
                continue
            if current_label is not None:
                m = site_re.match(line)
                if m:
                    current_sites.append(float(m.group(2)))

    if current_label is not None and current_sites:
        blocks[current_label] = np.array(current_sites)

    return blocks


# ---------------------------------------------------------------------------
# 2. Order labels
# ---------------------------------------------------------------------------

def ordered_labels(blocks):
    def key(lbl):
        if lbl == "Baseline":
            return 0
        return int(re.search(r"\d+", lbl).group())
    return sorted(blocks.keys(), key=key)


# ---------------------------------------------------------------------------
# 3. Analysis helpers
# ---------------------------------------------------------------------------

def delta_matrix(blocks, labels):
    """Return (n_sites, n_deltas) matrix of Δ log-L between consecutive labels."""
    arrays = [blocks[l] for l in labels]
    n_sites = arrays[0].shape[0]
    deltas = np.column_stack([arrays[i+1] - arrays[i] for i in range(len(arrays)-1)])
    return deltas  # positive = improvement (less negative log-L)


def top_movers(delta_col, n=20):
    """Indices of the n sites with the largest absolute Δ."""
    return np.argsort(np.abs(delta_col))[::-1][:n]


# ---------------------------------------------------------------------------
# 4. Plot
# ---------------------------------------------------------------------------

def plot(blocks, labels, out_path="figures/site_likelihood_analysis.png"):
    deltas = delta_matrix(blocks, labels)
    n_sites, n_wgd = deltas.shape
    sites = np.arange(n_sites)

    wgd_labels = [f"WGD #{i+1}" for i in range(n_wgd)]

    fig = plt.figure(figsize=(16, 4 * (n_wgd + 2)))
    n_rows = n_wgd + 2
    gs = gridspec.GridSpec(n_rows, 1, figure=fig, hspace=0.55)

    # ---- Row 0: baseline log-likelihoods ----
    ax0 = fig.add_subplot(gs[0])
    ax0.bar(sites, blocks["Baseline"], width=1, color="steelblue", alpha=0.7)
    ax0.set_title("Baseline per-site log-likelihood", fontsize=11)
    ax0.set_xlabel("Site"); ax0.set_ylabel("log L")
    ax0.grid(axis="y", alpha=0.3)

    # ---- Rows 1..n_wgd: Δ log-L for each accepted WGD ----
    vmax = np.abs(deltas).max()
    for k in range(n_wgd):
        ax = fig.add_subplot(gs[k + 1])
        col = deltas[:, k]
        colors = np.where(col >= 0, "forestgreen", "tomato")
        ax.bar(sites, col, width=1, color=colors, alpha=0.8)
        ax.axhline(0, color="black", lw=0.7)
        total = col.sum()
        ax.set_title(
            f"Δ log-L after {wgd_labels[k]}  (total Δ = {total:.2f})",
            fontsize=11
        )
        ax.set_xlabel("Site"); ax.set_ylabel("Δ log L")
        ax.grid(axis="y", alpha=0.3)

        # annotate the top 5 movers
        for idx in top_movers(col, n=5):
            ax.annotate(
                str(idx),
                xy=(idx, col[idx]),
                xytext=(0, 6 if col[idx] >= 0 else -12),
                textcoords="offset points",
                ha="center", fontsize=7, color="black"
            )

    # ---- Last row: cumulative Δ log-L (sum across all WGDs) ----
    ax_cum = fig.add_subplot(gs[-1])
    cumulative = deltas.sum(axis=1)
    colors = np.where(cumulative >= 0, "forestgreen", "tomato")
    ax_cum.bar(sites, cumulative, width=1, color=colors, alpha=0.8)
    ax_cum.axhline(0, color="black", lw=0.7)
    ax_cum.set_title(
        f"Cumulative Δ log-L across all WGDs  (total = {cumulative.sum():.2f})",
        fontsize=11
    )
    ax_cum.set_xlabel("Site"); ax_cum.set_ylabel("Δ log L")
    ax_cum.grid(axis="y", alpha=0.3)

    fig.suptitle("Per-site log-likelihood changes at each WGD acceptance", fontsize=13, fontweight="bold")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {out_path}")


# ---------------------------------------------------------------------------
# 5. Text summary
# ---------------------------------------------------------------------------

def summarize(blocks, labels):
    deltas = delta_matrix(blocks, labels)
    n_sites, n_wgd = deltas.shape
    wgd_labels = [f"WGD #{i+1}" for i in range(n_wgd)]

    print(f"\n{'='*60}")
    print(f"Sites: {n_sites}    WGDs accepted: {n_wgd}")
    print(f"{'='*60}\n")

    for k in range(n_wgd):
        col = deltas[:, k]
        improved = (col > 0).sum()
        worsened = (col < 0).sum()
        print(f"{wgd_labels[k]}:  total Δ = {col.sum():.3f}  |  "
              f"improved: {improved}  worsened: {worsened}")
        top = top_movers(col, n=10)
        print(f"  Top 10 movers (by |Δ|):")
        for idx in top:
            direction = "↑" if col[idx] >= 0 else "↓"
            print(f"    site {idx:5d}:  {col[idx]:+.3f} {direction}")
        print()

    print("Cumulative Δ log-L per site (top 15 movers overall):")
    cumulative = deltas.sum(axis=1)
    for idx in top_movers(cumulative, n=15):
        print(f"  site {idx:5d}:  {cumulative[idx]:+.3f}")
    print()


# ---------------------------------------------------------------------------
# 6. Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    log_path = sys.argv[1] if len(sys.argv) > 1 else LOG_FILE
    blocks = parse_log(log_path)
    labels = ordered_labels(blocks)

    print(f"Parsed {len(labels)} iterations: {labels}")
    summarize(blocks, labels)
    plot(blocks, labels)
