"""
Create publication-quality visualizations for the permutation test analysis.

Usage:
    Run from the project root (the parent of code/, results/, figures/):
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

# Read observed values from results JSON (avoid hardcoding)
with open(os.path.join(RESULTS_DIR, 'permutation_results.json')) as f:
    results_data = json.load(f)
observed_50 = results_data['results']['50mg_vs_placebo']['observed']
observed_mitt_val = results_data['results']['mitt_vs_placebo']['observed']

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

observed_diff = observed_50
# Two-sided: count permutations with |diff| >= |observed|
n_extreme = int(np.sum(np.abs(perm_dist_50) >= abs(observed_diff)))
p_two = n_extreme / len(perm_dist_50)

unique_vals, counts = np.unique(np.round(perm_dist_50, 1), return_counts=True)
colors = ['#c0392b' if abs(v) >= abs(observed_diff) else '#bdc3c7' for v in unique_vals]
ax.bar(unique_vals, counts, width=6, color=colors, edgecolor='white', linewidth=0.5)

ax.axvline(observed_diff, color='#c0392b', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-observed_diff, color='#c0392b', linewidth=1, linestyle=':', alpha=0.5, zorder=5)
ax.annotate(f'Observed Δ = +{observed_diff:.1f}m\ntwo-sided p = {p_two:.3f}',
            xy=(observed_diff, max(counts)*0.85), xytext=(observed_diff + 20, max(counts)*0.9),
            fontsize=12, fontweight='bold', color='#c0392b',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5),
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#c0392b', alpha=0.9))

ax.annotate('MMRM: p ≈ 0.56\n(derived; not significant)',
            xy=(-20, max(counts)*0.65), fontsize=10, fontstyle='italic', color='#7f8c8d',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', edgecolor='#bdc3c7'))

ax.set_xlabel('Difference in mean 6MWT change (50mg − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 70)', fontsize=12)
ax.set_title('Exact Permutation Distribution: 50 mg/kg Eteplirsen vs Placebo\n'
             'Study 201 (Mendell et al. 2013), Week 24 6MWT Change', fontsize=13, fontweight='bold')

red_patch = mpatches.Patch(color='#c0392b', label=f'|Δ| ≥ observed ({n_extreme}/70 = two-sided p={p_two:.3f})')
grey_patch = mpatches.Patch(color='#bdc3c7', label=f'|Δ| < observed ({len(perm_dist_50) - n_extreme}/70)')
ax.legend(handles=[red_patch, grey_patch], loc='upper left', framealpha=0.9)

plt.tight_layout()
plt.savefig(os.path.join(PNG_DIR, 'fig1_permutation_distribution_50mg.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig1_permutation_distribution_50mg.pdf'), bbox_inches='tight')
print("Figure 1 saved.")
plt.close()

# =========================================================================
# FIGURE 2: Permutation Distribution — mITT
# =========================================================================
fig, ax = plt.subplots(figsize=(10, 5.5))

observed_mitt = observed_mitt_val
n_extreme_mitt = int(np.sum(np.abs(perm_dist_mitt) >= abs(observed_mitt)))
p_two_mitt = n_extreme_mitt / len(perm_dist_mitt)

unique_vals2, counts2 = np.unique(np.round(perm_dist_mitt, 1), return_counts=True)
colors2 = ['#c0392b' if abs(v) >= abs(observed_mitt) else '#bdc3c7' for v in unique_vals2]
ax.bar(unique_vals2, counts2, width=8, color=colors2, edgecolor='white', linewidth=0.5)

ax.axvline(observed_mitt, color='#c0392b', linewidth=2, linestyle='--', zorder=5)
ax.axvline(-observed_mitt, color='#c0392b', linewidth=1, linestyle=':', alpha=0.5, zorder=5)
ax.annotate(f'Observed Δ = +{observed_mitt:.1f}m\ntwo-sided p = {p_two_mitt:.3f}',
            xy=(observed_mitt, max(counts2)*0.8), xytext=(observed_mitt + 25, max(counts2)*0.85),
            fontsize=12, fontweight='bold', color='#c0392b',
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5),
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#c0392b', alpha=0.9))

ax.set_xlabel('Difference in mean 6MWT change (mITT eteplirsen − placebo), metres', fontsize=12)
ax.set_ylabel('Number of permutations (out of 210)', fontsize=12)
ax.set_title('Exact Permutation Distribution: mITT Eteplirsen (n=6) vs Placebo (n=4)\n'
             'Excluding 2 patients who lost ambulation', fontsize=13, fontweight='bold')

red_patch2 = mpatches.Patch(color='#c0392b', label=f'|Δ| ≥ observed ({n_extreme_mitt}/210 = two-sided p={p_two_mitt:.3f})')
grey_patch2 = mpatches.Patch(color='#bdc3c7', label=f'|Δ| < observed ({len(perm_dist_mitt) - n_extreme_mitt}/210)')
ax.legend(handles=[red_patch2, grey_patch2], loc='upper left', framealpha=0.9)

plt.tight_layout()
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
ax1.hist(sensitivity_50, bins=30, color='#3498db', edgecolor='white', alpha=0.8)
ax1.axvline(0.05, color='#c0392b', linewidth=2, linestyle='--', label='α = 0.05')
ax1.axvline(p_two, color='#2c3e50', linewidth=2, linestyle='-', label=f'Our estimate ({p_two:.3f})')
ax1.set_xlabel('Permutation p-value', fontsize=11)
ax1.set_ylabel('Count (of 1000 simulations)', fontsize=11)
ax1.set_title(f'50 mg/kg vs Placebo\n{pct_50:.0f}% below α=0.05', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)

ax2 = axes[1]
ax2.hist(sensitivity_mitt, bins=30, color='#2ecc71', edgecolor='white', alpha=0.8)
ax2.axvline(0.05, color='#c0392b', linewidth=2, linestyle='--', label='α = 0.05')
ax2.axvline(p_two_mitt, color='#2c3e50', linewidth=2, linestyle='-', label=f'Our estimate ({p_two_mitt:.3f})')
ax2.set_xlabel('Permutation p-value', fontsize=11)
ax2.set_ylabel('Count (of 1000 simulations)', fontsize=11)
ax2.set_title(f'mITT Eteplirsen vs Placebo\n{pct_mitt:.0f}% below α=0.05', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)

fig.suptitle('Sensitivity Analysis: Permutation p-values across 1000 plausible\n'
             'individual patient allocations (group means held fixed)',
             fontsize=13, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig(os.path.join(PNG_DIR, 'fig3_sensitivity_analysis.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig3_sensitivity_analysis.pdf'), bbox_inches='tight')
print("Figure 3 saved.")
plt.close()

# =========================================================================
# FIGURE 4: The p-value comparison — the headline result
# =========================================================================
fig, ax = plt.subplots(figsize=(8, 5))

comparisons = ['50mg vs Placebo\n(dose arm)', 'mITT Eteplirsen\nvs Placebo', 'All Eteplirsen\n(ITT) vs Placebo']
p_two_itt = results_data['results']['itt_all_vs_placebo']['p_two']
perm_pvals = [p_two, p_two_mitt, p_two_itt]
mmrm_pvals = [0.563, None, None]

x = np.arange(len(comparisons))
width = 0.3

bars1 = ax.bar(x - width/2, perm_pvals, width, label='Permutation test (two-sided)', color='#2ecc71', edgecolor='white')
mmrm_x = [0]
mmrm_y = [0.563]
bars2 = ax.bar([x[0] + width/2], mmrm_y, width, label='MMRM (derived, two-sided)', color='#e74c3c', edgecolor='white')

ax.axhline(0.05, color='black', linewidth=1, linestyle=':', alpha=0.5, label='α = 0.05')

for i, (bar, pv) in enumerate(zip(bars1, perm_pvals)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
            f'p={pv:.3f}', ha='center', fontsize=10, fontweight='bold', color='#27ae60')

ax.text(bars2[0].get_x() + bars2[0].get_width()/2, bars2[0].get_height() + 0.015,
        'p≈0.56', ha='center', fontsize=10, fontweight='bold', color='#c0392b')

ax.set_xticks(x)
ax.set_xticklabels(comparisons, fontsize=10)
ax.set_ylabel('p-value (two-sided)', fontsize=12)
ax.set_title('The Permutation Test Nobody Ran\nEteplirsen Study 201 — Week 24 6MWT', fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.1)
ax.legend(loc='upper right', fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(PNG_DIR, 'fig4_pvalue_comparison.png'), dpi=200, bbox_inches='tight')
plt.savefig(os.path.join(PDF_DIR, 'fig4_pvalue_comparison.pdf'), bbox_inches='tight')
print("Figure 4 saved.")
plt.close()

print(f"\nAll visualizations saved to {PNG_DIR}/ and {PDF_DIR}/")
