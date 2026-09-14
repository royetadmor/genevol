from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print("Usage: python3 plot_tree.py '<newick string>' [output.png]")
    sys.exit(1)

NEWICK = sys.argv[1]
OUTPATH = sys.argv[2] if len(sys.argv) >= 3 else "figures/tree.png"

tree = Phylo.read(StringIO(NEWICK), "newick")


def compute_positions(tree):
    """Return {clade: (x, y)} mirroring Bio.Phylo's phylogram layout."""
    # Assign y: leaf order in pre-order traversal, 1-indexed from bottom
    terminals = tree.get_terminals()
    y_of = {t: i + 1 for i, t in enumerate(terminals)}

    def _y(clade):
        if clade.is_terminal():
            return y_of[clade]
        child_ys = [_y(c) for c in clade.clades]
        y = sum(child_ys) / len(child_ys)
        y_of[clade] = y
        return y

    _y(tree.root)

    # Assign x: cumulative branch length from root
    x_of = {}

    def _x(clade, parent_x):
        x_of[clade] = parent_x + (clade.branch_length or 0)
        for child in clade.clades:
            _x(child, x_of[clade])

    _x(tree.root, 0)

    return {c: (x_of[c], y_of[c]) for c in x_of}


positions = compute_positions(tree)

# Collect WGD nodes: non-root clades with branch_length == 0
wgd_clades = [
    c for c in tree.find_clades()
    if c is not tree.root and (c.branch_length or 0) == 0
]

fig, ax = plt.subplots(figsize=(12, 4))
Phylo.draw(tree, axes=ax, do_show=False,
           label_func=lambda c: c.name if c.is_terminal() else "")

# Group markers by (x, y) and offset overlapping ones vertically
from collections import defaultdict
groups = defaultdict(list)
for clade in wgd_clades:
    x, y = positions[clade]
    groups[(round(x, 6), round(y, 4))].append((x, y))

JITTER = 0.18  # vertical offset between stacked markers
for (_, _), pts in groups.items():
    n = len(pts)
    offsets = [JITTER * (i - (n - 1) / 2) for i in range(n)]
    for (x, y), dy in zip(pts, offsets):
        ax.plot(x, y + dy, marker="D", ms=6, color="#27ae60", zorder=5, clip_on=False)

ax.set_title("")
ax.set_xlabel("Branch length", fontsize=11)
ax.set_ylabel("")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_yticks([])

plt.tight_layout()
plt.savefig(OUTPATH, dpi=150, bbox_inches="tight")
print(f"Saved {OUTPATH}")
