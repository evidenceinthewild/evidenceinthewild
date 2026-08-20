"""
Create title-free figures for the Ghost CMS blog post.
Titles and details go in Ghost's image caption fields instead.

Usage:
    Invoke this script by path from any working directory. Inputs and outputs
    are resolved relative to the script location. From the project root:
        python code/create_blog_figures.py

    Requires: numpy, matplotlib
    Input:  results/*.npy files (produced by eteplirsen_permutation_analysis.py
            and sensitivity_analysis.py)
    Output: figures/blog/*.png
"""
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
BLOG_DIR = os.path.join(PROJECT_ROOT, 'figures', 'blog')
os.makedirs(BLOG_DIR, exist_ok=True)

perm_dist_50 = np.load(os.path.join(RESULTS_DIR, 'perm_dist_50v_placebo.npy'))
perm_dist_mitt = np.load(os.path.join(RESULTS_DIR, 'perm_dist_mitt.npy'))
sensitivity_50 = np.load(os.path.join(RESULTS_DIR, 'sensitivity_50v_pvals.npy'))
sensitivity_mitt = np.load(os.path.join(RESULTS_DIR, 'sensitivity_mitt_pvals.npy'))

with open(os.path.join(RESULTS_DIR, 'permutation_results.json')) as f:
    results_data = json.load(f)
displayed_50 = results_data['results']['50mg_vs_placebo']['observed']
displayed_mitt_val = results_data['results']['mitt_vs_placebo']['observed']

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
})


# Figure 1: 50 mg/kg vs placebo permutation distribution.
fig, ax = plt.subplots(figsize=(10, 5))

displayed_diff = displayed_50
n_extreme = int(np.sum(np.abs(perm_dist_50) >= abs(displayed_diff)))
p_two = n_extreme / len(perm_dist_50)

bins = np.arange(np.floor(perm_dist_50.min() / 10) * 10 - 5,
                 np.ceil(perm_dist_50.max() / 10) * 10 + 15, 10)
_, bin_edges, patches = ax.hist(
    perm_dist_50, bins=bins, color='#bdc3c7', edgecolor='white', linewidth=0.5
)
for patch, left, right in zip(patches, bin_edges[:-1], bin_edges[1:]):
    in_bin = perm_dist_50[(perm_dist_50 >= left) & (perm_dist_50 < right)]
    if np.any(np.abs(in_bin) >= abs(displayed_diff)):
        patch.set_facecolor('#2f6b8a')

ax.axvline(displayed_diff, color='#2f6b8a', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-displayed_diff, color='#2f6b8a', linewidth=1, linestyle=':', alpha=0.6, zorder=5)
ax.annotate(
    f'Displayed construction\nΔ = +{displayed_diff:.1f}m; p = {p_two:.3f}',
    xy=(displayed_diff, ax.get_ylim()[1] * 0.75),
    xytext=(displayed_diff + 18, ax.get_ylim()[1] * 0.85),
    fontsize=11,
    fontweight='bold',
    color='#2f6b8a',
    arrowprops=dict(arrowstyle='->', color='#2f6b8a', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#2f6b8a', alpha=0.9),
)
ax.annotate(
    'MMRM: p ≈ 0.56\n(derived; different analysis)',
    xy=(-75, ax.get_ylim()[1] * 0.92),
    fontsize=10,
    fontstyle='italic',
    color='#7f8c8d',
    bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', edgecolor='#bdc3c7'),
)
ax.set_xlabel('Difference in mean 6MWT change (50mg − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 70)', fontsize=12)

focus_patch = mpatches.Patch(
    color='#2f6b8a', label=f'|Δ| ≥ displayed split ({n_extreme}/70; p={p_two:.3f})'
)
grey_patch = mpatches.Patch(
    color='#bdc3c7', label=f'|Δ| < displayed split ({len(perm_dist_50) - n_extreme}/70)'
)
ax.legend(
    handles=[focus_patch, grey_patch],
    loc='upper center',
    framealpha=0.9,
    bbox_to_anchor=(0.5, -0.15),
    ncol=2,
)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
plt.savefig(os.path.join(BLOG_DIR, 'fig1_permutation_distribution_50mg.png'), dpi=200, bbox_inches='tight')
print('Blog figure 1 saved.')
plt.close()


# Figure 2: exploratory mITT permutation distribution.
fig, ax = plt.subplots(figsize=(10, 5))

displayed_mitt = displayed_mitt_val
n_extreme_mitt = int(np.sum(np.abs(perm_dist_mitt) >= abs(displayed_mitt)))
p_two_mitt = n_extreme_mitt / len(perm_dist_mitt)

bins2 = np.arange(np.floor(perm_dist_mitt.min() / 10) * 10 - 5,
                  np.ceil(perm_dist_mitt.max() / 10) * 10 + 15, 10)
_, bin_edges2, patches2 = ax.hist(
    perm_dist_mitt, bins=bins2, color='#bdc3c7', edgecolor='white', linewidth=0.5
)
for patch, left, right in zip(patches2, bin_edges2[:-1], bin_edges2[1:]):
    in_bin = perm_dist_mitt[(perm_dist_mitt >= left) & (perm_dist_mitt < right)]
    if np.any(np.abs(in_bin) >= abs(displayed_mitt)):
        patch.set_facecolor('#c6922f')

ax.axvline(displayed_mitt, color='#c6922f', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-displayed_mitt, color='#c6922f', linewidth=1, linestyle=':', alpha=0.6, zorder=5)
ax.annotate(
    f'Displayed mITT construction\nΔ = +{displayed_mitt:.1f}m; p = {p_two_mitt:.3f}',
    xy=(displayed_mitt, ax.get_ylim()[1] * 0.75),
    xytext=(displayed_mitt + 20, ax.get_ylim()[1] * 0.85),
    fontsize=11,
    fontweight='bold',
    color='#8a621f',
    arrowprops=dict(arrowstyle='->', color='#c6922f', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#c6922f', alpha=0.9),
)
ax.set_xlabel('Difference in mean 6MWT change (exploratory mITT − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 210)', fontsize=12)

focus_patch2 = mpatches.Patch(
    color='#c6922f', label=f'|Δ| ≥ displayed split ({n_extreme_mitt}/210; p={p_two_mitt:.3f})'
)
grey_patch2 = mpatches.Patch(
    color='#bdc3c7', label=f'|Δ| < displayed split ({len(perm_dist_mitt) - n_extreme_mitt}/210)'
)
ax.legend(
    handles=[focus_patch2, grey_patch2],
    loc='upper center',
    framealpha=0.9,
    bbox_to_anchor=(0.5, -0.15),
    ncol=2,
)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
plt.savefig(os.path.join(BLOG_DIR, 'fig2_permutation_distribution_mitt.png'), dpi=200, bbox_inches='tight')
print('Blog figure 2 saved.')
plt.close()


# Figure 3: generator-based sensitivity analysis.
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

pct_50 = 100 * np.mean(sensitivity_50 < 0.05)
pct_mitt = 100 * np.mean(sensitivity_mitt < 0.05)

ax1 = axes[0]
ax1.hist(sensitivity_50, bins=30, color='#4c78a8', edgecolor='white', alpha=0.85)
ax1.axvline(0.05, color='#b65f24', linewidth=2, linestyle='--', label='0.05 reference')
ax1.axvline(p_two, color='#263238', linewidth=2, linestyle='-', label=f'Displayed construction ({p_two:.3f})')
ax1.set_xlabel('Permutation p-value', fontsize=11)
ax1.set_ylabel('Generated configurations (n=1000)', fontsize=11)
ax1.set_title('50 mg/kg vs placebo', fontsize=12, fontweight='bold')
ax1.text(
    0.98,
    0.80,
    f'{pct_50:.0f}% below 0.05\nmedian = {np.median(sensitivity_50):.3f}',
    transform=ax1.transAxes,
    ha='right',
    fontsize=10,
    bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#bdc3c7'),
)
ax1.legend(fontsize=9, loc='center right')

ax2 = axes[1]
ax2.hist(sensitivity_mitt, bins=30, color='#d4a72c', edgecolor='white', alpha=0.85)
ax2.axvline(0.05, color='#b65f24', linewidth=2, linestyle='--', label='0.05 reference')
ax2.axvline(p_two_mitt, color='#263238', linewidth=2, linestyle='-', label=f'Displayed construction ({p_two_mitt:.3f})')
ax2.set_xlabel('Permutation p-value', fontsize=11)
ax2.set_ylabel('Generated configurations (n=1000)', fontsize=11)
ax2.set_title('Exploratory mITT eteplirsen vs placebo', fontsize=12, fontweight='bold')
ax2.text(
    0.98,
    0.80,
    f'{pct_mitt:.0f}% below 0.05\nmedian = {np.median(sensitivity_mitt):.3f}',
    transform=ax2.transAxes,
    ha='right',
    fontsize=10,
    bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#bdc3c7'),
)
ax2.legend(fontsize=9, loc='center right')

plt.tight_layout()
fig.subplots_adjust(bottom=0.18)
fig.text(
    0.5,
    0.02,
    'Selected group-mean targets held fixed; ranges are ad hoc. '
    'Left-panel p-values occur in steps of 2/70.',
    ha='center',
    fontsize=9,
    color='#5f6b70',
)
plt.savefig(os.path.join(BLOG_DIR, 'fig3_sensitivity_analysis.png'), dpi=200, bbox_inches='tight')
print('Blog figure 3 saved.')
plt.close()

print(f'\nBlog figures saved to {BLOG_DIR}/')
