"""Conceptual schematic: the borrowing dial (sigma = between-basket SD).
Not a data figure -- a diagram of how one hyperprior sets the degree of shrinkage,
from full pooling (sigma -> 0) to no borrowing (sigma -> infinity), with the paper's
two simulation settings, Half-N(0.3) and Half-N(3), marked."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

INK="#20201E"; RUST="#C8613A"; TEAL="#3E7C7B"; GOLD="#B08A3E"; GREY="#9A948C"
LGREY="#D9D5CE"
plt.rcParams.update({'font.family':'DejaVu Sans','axes.edgecolor':INK,'figure.dpi':150})

fig, ax = plt.subplots(figsize=(10.6, 5.0))
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')

mean_y   = 0.62      # common (pooled) mean line
bar_y    = 0.30      # the dial bar
max_amp  = 0.110     # max vertical spread of basket dots

# four settings along the dial: x position, spread factor, label lines, highlight colour
settings = [
    (0.19, 0.00, ["σ → 0", "full pooling"],               None),
    (0.42, 0.34, ["Half-N(0.3)", "strong borrowing"],     TEAL),
    (0.66, 0.82, ["Half-N(3)", "moderate borrowing"],     GOLD),
    (0.90, 1.21, ["σ → ∞", "no borrowing", "(separate)"], None),
]
dot_dx = np.array([-0.028, -0.010, 0.010, 0.028])
base   = np.array([-1.5, -0.5, 0.5, 1.5])

# common-mean reference line
ax.plot([0.05, 0.95], [mean_y, mean_y], ls=(0,(4,3)), lw=1.1, color=GREY, zorder=1)
ax.text(0.05, mean_y+0.022, "common mean", fontsize=8.5, color=GREY, va='bottom')

# basket clusters
for x, spread, labels, hl in settings:
    ys = mean_y + base*spread*max_amp
    for k in range(4):
        ax.plot([x+dot_dx[k], x+dot_dx[k]], [mean_y, ys[k]], lw=0.9, color=LGREY, zorder=2)
    ax.scatter(x+dot_dx, ys, s=42, color=INK, zorder=3, edgecolor='white', linewidth=0.6)

# the dial bar (double-headed continuum)
ax.annotate("", xy=(0.95, bar_y), xytext=(0.05, bar_y),
            arrowprops=dict(arrowstyle='<|-|>', color=GREY, lw=3))
ax.text(0.05, bar_y-0.075, "more shrinkage", fontsize=8.5, color=GREY, ha='left', style='italic')
ax.text(0.95, bar_y-0.075, "more separation", fontsize=8.5, color=GREY, ha='right', style='italic')

# ticks + labels; highlight the two simulation settings
for x, spread, labels, hl in settings:
    ax.plot([x, x], [bar_y-0.018, bar_y+0.018], lw=2, color=INK, zorder=4)
    col = hl if hl else INK
    weight = 'bold' if hl else 'normal'
    for i, ln in enumerate(labels):
        yy = bar_y - 0.11 - i*0.052
        fs = 10.5 if i == 0 else 9.2
        ax.text(x, yy, ln, ha='center', va='top', fontsize=fs, color=col,
                fontweight=(weight if i < 2 else 'normal'))

# highlight band tying the two paper settings to "the borrowing dial"
ax.annotate("", xy=(0.66, bar_y+0.06), xytext=(0.42, bar_y+0.06),
            arrowprops=dict(arrowstyle='<->', color=RUST, lw=1.3))
ax.text(0.54, bar_y+0.085, "the paper's two settings differ only here",
        ha='center', fontsize=8.6, color=RUST, style='italic')

# titles
ax.text(0.02, 0.965, "The borrowing dial", fontsize=15, color=INK, fontweight='bold', va='top')
ax.text(0.02, 0.905,
        "One hyperprior, the prior on σ, the between-basket SD, sets how far "
        "each basket is pulled toward the common mean.",
        fontsize=9.6, color=INK, va='top')

plt.tight_layout()
output_path = Path(__file__).resolve().with_name("figure_borrowing_dial.png")
plt.savefig(output_path, bbox_inches='tight')
print(f"wrote {output_path}")
