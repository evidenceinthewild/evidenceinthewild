"""
Create publication-quality visualizations for the permutation test analysis.

Usage:
    Invoke this script by path from any working directory. Inputs and outputs
    are resolved relative to the script location. From the project root:
        python code/create_visualizations.py

    Requires: numpy, matplotlib
    Input:  results/*.npy files (produced by eteplirsen_permutation_analysis.py
            and sensitivity_analysis.py)
    Output: figures/png/*.png and figures/pdf/*.pdf
"""
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# PATH SETUP
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
PNG_DIR = os.path.join(PROJECT_ROOT, 'figures', 'png')
PDF_DIR = os.path.join(PROJECT_ROOT, 'figures', 'pdf')
os.makedirs(PNG_DIR, exist_ok=True)
os.makedirs(PDF_DIR, exist_ok=True)

# Load saved data
import json

perm_dist_50 = np.load(os.path.join(RESULTS_DIR, 'perm_dist_50v_placebo.npy'))
perm_dist_mitt = np.load(os.path.join(RESULTS_DIR, 'perm_dist_mitt.npy'))
sensitivity_50 = np.load(os.path.join(RESULTS_DIR, 'sensitivity_50v_pvals.npy'))
sensitivity_mitt = np.load(os.path.join(RESULTS_DIR, 'sensitivity_mitt_pvals.npy'))

# Read values for the displayed construction from results JSON (avoid hardcoding)
with open(os.path.join(RESULTS_DIR, 'permutation_results.json')) as f:
    results_data = json.load(f)
displayed_50 = results_data['results']['50mg_vs_placebo']['observed']
displayed_mitt_val = results_data['results']['mitt_vs_placebo']['observed']

# Style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.facecolor': 'white',
})

# =========================================================================
# FIGURE 1: Permutation Distribution — 50mg vs Placebo
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 5.5))

displayed_diff = displayed_50
# Two-sided: count assignments with |diff| >= the displayed split
n_extreme = int(np.sum(np.abs(perm_dist_50) >= abs(displayed_diff)))
p_two = n_extreme / len(perm_dist_50)

# Bin edges: 10m-wide bins spanning the range of permutation diffs
bins = np.arange(np.floor(perm_dist_50.min() / 10) * 10 - 5,
                 np.ceil(perm_dist_50.max() / 10) * 10 + 15, 10)
# Color each bar if it contains a value at least as extreme as the displayed split.
n_vals, bin_edges, patches = ax.hist(perm_dist_50, bins=bins, color='#bdc3c7',
                                      edgecolor='white', linewidth=0.5)
for patch, left, right in zip(patches, bin_edges[:-1], bin_edges[1:]):
    centre = (left + right) / 2
    if abs(centre) >= abs(displayed_diff) - 5:  # bin contains extreme values
        # Check if any actual permutation value in this bin is extreme
        in_bin = perm_dist_50[(perm_dist_50 >= left) & (perm_dist_50 < right)]
        if np.any(np.abs(in_bin) >= abs(displayed_diff)):
            patch.set_facecolor('#2f6b8a')

ax.axvline(displayed_diff, color='#2f6b8a', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-displayed_diff, color='#2f6b8a', linewidth=1, linestyle=':', alpha=0.6, zorder=5)
ax.annotate(f'Displayed construction\nΔ = +{displayed_diff:.1f}m; p = {p_two:.3f}',
            xy=(displayed_diff, ax.get_ylim()[1]*0.75),
            xytext=(displayed_diff + 18, ax.get_ylim()[1]*0.85),
            fontsize=11, fontweight='bold', color='#2f6b8a',
            arrowprops=dict(arrowstyle='->', color='#2f6b8a', lw=1.5),
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#2f6b8a', alpha=0.9))

ax.annotate('MMRM: p ≈ 0.56\n(derived; different analysis)',
            xy=(-75, ax.get_ylim()[1]*0.92), fontsize=10, fontstyle='italic', color='#7f8c8d',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', edgecolor='#bdc3c7'))

ax.set_xlabel('Difference in mean 6MWT change (50mg − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 70)', fontsize=12)
ax.set_title('Permutation distribution: 50 mg/kg eteplirsen vs placebo\n'
             'Constructed Week 24 outcomes; conditional on 70 assumed allocations',
             fontsize=13, fontweight='bold')

focus_patch = mpatches.Patch(color='#2f6b8a', label=f'|Δ| ≥ displayed split ({n_extreme}/70; p={p_two:.3f})')
grey_patch = mpatches.Patch(color='#bdc3c7', label=f'|Δ| < displayed split ({len(perm_dist_50) - n_extreme}/70)')
ax.legend(handles=[focus_patch, grey_patch], loc='upper center', framealpha=0.9,
          bbox_to_anchor=(0.5, -0.12), ncol=2)

plt.tight_layout()
fig.subplots_adjust(bottom=0.2)
plt.savefig(os.path.join(PNG_DIR, 'fig1_permutation_distribution_50mg.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig1_permutation_distribution_50mg.pdf'), bbox_inches='tight')
print("Figure 1 saved.")
plt.close()

# =========================================================================
# FIGURE 2: Permutation Distribution — mITT
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 5.5))

displayed_mitt = displayed_mitt_val
n_extreme_mitt = int(np.sum(np.abs(perm_dist_mitt) >= abs(displayed_mitt)))
p_two_mitt = n_extreme_mitt / len(perm_dist_mitt)

# Bin edges: 10m-wide bins
bins2 = np.arange(np.floor(perm_dist_mitt.min() / 10) * 10 - 5,
                  np.ceil(perm_dist_mitt.max() / 10) * 10 + 15, 10)
n_vals2, bin_edges2, patches2 = ax.hist(perm_dist_mitt, bins=bins2, color='#bdc3c7',
                                         edgecolor='white', linewidth=0.5)
for patch, left, right in zip(patches2, bin_edges2[:-1], bin_edges2[1:]):
    in_bin = perm_dist_mitt[(perm_dist_mitt >= left) & (perm_dist_mitt < right)]
    if np.any(np.abs(in_bin) >= abs(displayed_mitt)):
        patch.set_facecolor('#c6922f')

ax.axvline(displayed_mitt, color='#c6922f', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-displayed_mitt, color='#c6922f', linewidth=1, linestyle=':', alpha=0.6, zorder=5)
ax.annotate(f'Displayed mITT construction\nΔ = +{displayed_mitt:.1f}m; p = {p_two_mitt:.3f}',
            xy=(displayed_mitt, ax.get_ylim()[1]*0.75),
            xytext=(displayed_mitt + 20, ax.get_ylim()[1]*0.85),
            fontsize=11, fontweight='bold', color='#8a621f',
            arrowprops=dict(arrowstyle='->', color='#c6922f', lw=1.5),
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#c6922f', alpha=0.9))

ax.set_xlabel('Difference in mean 6MWT change (mITT eteplirsen − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 210)', fontsize=12)
ax.set_title('Exploratory permutation distribution: mITT eteplirsen vs placebo\n'
             'Constructed outcomes after post-randomisation exclusion',
             fontsize=13, fontweight='bold')

focus_patch2 = mpatches.Patch(color='#c6922f', label=f'|Δ| ≥ displayed split ({n_extreme_mitt}/210; p={p_two_mitt:.3f})')
grey_patch2 = mpatches.Patch(color='#bdc3c7', label=f'|Δ| < displayed split ({len(perm_dist_mitt) - n_extreme_mitt}/210)')
ax.legend(handles=[focus_patch2, grey_patch2], loc='upper center', framealpha=0.9,
          bbox_to_anchor=(0.5, -0.12), ncol=2)

plt.tight_layout()
fig.subplots_adjust(bottom=0.2)
plt.savefig(os.path.join(PNG_DIR, 'fig2_permutation_distribution_mitt.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig2_permutation_distribution_mitt.pdf'), bbox_inches='tight')
print("Figure 2 saved.")
plt.close()

# =========================================================================
# FIGURE 3: Sensitivity Analysis
# =========================================================================
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
ax1.text(0.98, 0.80, f'{pct_50:.0f}% below 0.05\nmedian = {np.median(sensitivity_50):.3f}',
         transform=ax1.transAxes, ha='right', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#bdc3c7'))
ax1.legend(fontsize=9, loc='center right')

ax2 = axes[1]
ax2.hist(sensitivity_mitt, bins=30, color='#d4a72c', edgecolor='white', alpha=0.85)
ax2.axvline(0.05, color='#b65f24', linewidth=2, linestyle='--', label='0.05 reference')
ax2.axvline(p_two_mitt, color='#263238', linewidth=2, linestyle='-', label=f'Displayed construction ({p_two_mitt:.3f})')
ax2.set_xlabel('Permutation p-value', fontsize=11)
ax2.set_ylabel('Generated configurations (n=1000)', fontsize=11)
ax2.set_title('Exploratory mITT eteplirsen vs placebo', fontsize=12, fontweight='bold')
ax2.text(0.98, 0.80, f'{pct_mitt:.0f}% below 0.05\nmedian = {np.median(sensitivity_mitt):.3f}',
         transform=ax2.transAxes, ha='right', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#bdc3c7'))
ax2.legend(fontsize=9, loc='center right')

fig.suptitle('Generator-based sensitivity analysis\n'
             'Selected group-mean targets held fixed; individual ranges are ad hoc',
             fontsize=13, fontweight='bold', y=1.02)

plt.tight_layout()
fig.subplots_adjust(bottom=0.18)
fig.text(0.5, 0.02,
         'Selected group-mean targets held fixed; ranges are ad hoc. '
         'Left-panel p-values occur in steps of 2/70.',
         ha='center', fontsize=9, color='#5f6b70')
plt.savefig(os.path.join(PNG_DIR, 'fig3_sensitivity_analysis.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig3_sensitivity_analysis.pdf'), bbox_inches='tight')
print("Figure 3 saved.")
plt.close()

# =========================================================================
# FIGURE 4: Supplementary p-value comparison (not used in the current post)
# =========================================================================
fig, ax = plt.subplots(figsize=(8, 5))

comparisons = ['50mg vs Placebo\n(displayed construction)',
               'mITT Eteplirsen\n(exploratory subset)',
               'All Eteplirsen\n(displayed construction)']
p_two_itt = results_data['results']['itt_all_vs_placebo']['p_two']
perm_pvals = [p_two, p_two_mitt, p_two_itt]

x = np.arange(len(comparisons))
width = 0.3

bars1 = ax.bar(x - width/2, perm_pvals, width,
               label='Enumeration using constructed outcomes', color='#4c78a8', edgecolor='white')
bars2 = ax.bar([x[0] + width/2], [0.563], width,
               label='MMRM (derived; different analysis)', color='#d4a72c', edgecolor='white')

ax.set_yscale('log')
ax.axhline(0.05, color='black', linewidth=1, linestyle=':', alpha=0.5, label='α = 0.05')

for bar, pv in zip(bars1, perm_pvals):
    ax.text(bar.get_x() + bar.get_width()/2, pv * 1.5,
            f'p={pv:.3f}', ha='center', fontsize=10, fontweight='bold', color='#2f5f85')

ax.text(bars2[0].get_x() + bars2[0].get_width()/2, 0.563 * 1.2,
        'p≈0.56', ha='center', fontsize=10, fontweight='bold', color='#8a621f')

ax.set_xticks(x)
ax.set_xticklabels(comparisons, fontsize=10)
ax.set_ylabel('p-value (two-sided, log scale)', fontsize=12)
ax.set_title('P-values from different analysis specifications\n'
             'Constructed Week 24 outcomes; values are not directly comparable',
             fontsize=14, fontweight='bold')
ax.set_ylim(0.002, 2.0)
ax.set_yticks([0.005, 0.01, 0.05, 0.1, 0.5, 1.0])
ax.set_yticklabels(['0.005', '0.01', '0.05', '0.10', '0.50', '1.00'])
ax.legend(loc='upper center', fontsize=10, bbox_to_anchor=(0.5, -0.12), ncol=3)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
plt.savefig(os.path.join(PNG_DIR, 'fig4_pvalue_comparison.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig4_pvalue_comparison.pdf'), bbox_inches='tight')
print("Figure 4 saved.")
plt.close()

print(f"\nAll visualizations saved to {PNG_DIR}/ and {PDF_DIR}/")
