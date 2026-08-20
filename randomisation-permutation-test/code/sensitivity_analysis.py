"""
Sensitivity Analysis for Eteplirsen Permutation Tests
=====================================================
Examines sensitivity of the permutation p-value across generated individual
outcome configurations that preserve selected group-mean targets.

The generation ranges are ad hoc and do not define a probability model for
the unknown patient data. Reported percentages are properties of this
generator, not probabilities about the unobserved original-data p-value.

Usage:
    Invoke this script by path from any working directory. Inputs and outputs
    are resolved relative to the script location. From the project root:
        python code/sensitivity_analysis.py
"""

import numpy as np
from itertools import combinations
import os
import random

random.seed(42)
np.random.seed(42)

# =============================================================================
# PATH SETUP
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# Selected targets and constraints; see ../PROVENANCE.md.
# Patient 010: Week 24 change = -213 (reported in the paper).
# Patient 009: Week 24 change = -294 (reconstructed; not a reported value).
# 50mg group mean change = 12.5m (analysis input requiring source verification;
#                                not separately reported by Figure 6).
# Placebo group mean change = -55.0m (approximate graphical target).
# mITT eteplirsen mean change = +15.0m (approximate combined-cohort target).

def two_arm_perm_pvalue(treatment, control):
    """Exact two-sided permutation p-value."""
    all_vals = np.array(list(treatment) + list(control))
    n_t = len(treatment)
    observed = np.mean(treatment) - np.mean(control)
    count_extreme = 0
    total = 0
    for idx in combinations(range(len(all_vals)), n_t):
        trt = all_vals[list(idx)]
        ctrl_vals = [all_vals[i] for i in range(len(all_vals)) if i not in idx]
        diff = np.mean(trt) - np.mean(ctrl_vals)
        if abs(diff) >= abs(observed):
            count_extreme += 1
        total += 1
    return count_extreme / total

def generate_configurations(group_mean, n, constraints=None, n_samples=5000):
    """
    Generate individual change sets that sum to n * group_mean.

    The ranges are deliberately broad but ad hoc. They support a stress test;
    they are not a validated model for the missing patient outcomes.
    """
    target_sum = n * group_mean
    samples = []

    for _ in range(n_samples * 10):  # oversample then filter
        if constraints == 'placebo':
            # Placebo DMD patients: expect decline, range roughly -150 to +20
            vals = np.random.uniform(-150, 20, n)
        elif constraints == '50mg':
            # 50mg patients (remained ambulant): range roughly -50 to +60
            vals = np.random.uniform(-50, 60, n)
        elif constraints == '30mg_mitt':
            # 30mg mITT patients (remained ambulant): range roughly -30 to +60
            vals = np.random.uniform(-30, 60, n)
        else:
            vals = np.random.uniform(-100, 100, n)

        # Adjust to hit target sum
        current_sum = np.sum(vals)
        adjustment = (target_sum - current_sum) / n
        vals = vals + adjustment

        # Check all values remain within the configured ad hoc range.
        if constraints == 'placebo' and np.all(vals >= -200) and np.all(vals <= 50):
            samples.append(vals.tolist())
        elif constraints == '50mg' and np.all(vals >= -80) and np.all(vals <= 80):
            samples.append(vals.tolist())
        elif constraints == '30mg_mitt' and np.all(vals >= -60) and np.all(vals <= 80):
            samples.append(vals.tolist())

        if len(samples) >= n_samples:
            break

    return samples[:n_samples]

print("=" * 70)
print("GENERATOR-BASED SENSITIVITY ANALYSIS (two-sided p-values)")
print("Varying individual patient changes while preserving selected group-mean targets")
print("=" * 70)

# Generate many outcome configurations under the ad hoc ranges
n_sims = 5000
placebo_samples = generate_configurations(-55.0, 4, 'placebo', n_sims)
arm50_samples = generate_configurations(12.5, 4, '50mg', n_sims)
arm30_mitt_samples = generate_configurations(20.0, 2, '30mg_mitt', n_sims)
# mITT 30mg ambulant mean = (25 + 15) / 2 = 20.0 in the displayed construction;
# combined with the selected 50mg target, this yields the approximate +15m target.

n_available = min(len(placebo_samples), len(arm50_samples), len(arm30_mitt_samples))
print(f"\nGenerated {n_available} outcome configurations per arm")

# ---- Test A: 50mg vs placebo (key comparison) ----
print(f"\n--- 50 mg/kg vs Placebo (each n=4), two-sided ---")
pvals_50v = []
for i in range(min(1000, n_available)):
    pv = two_arm_perm_pvalue(arm50_samples[i], placebo_samples[i])
    pvals_50v.append(pv)
    if (i+1) % 200 == 0:
        print(f"  ... computed {i+1} simulations")

pvals_50v = np.array(pvals_50v)
print(f"\n  Simulations run: {len(pvals_50v)}")
print(f"  Mean p-value:    {np.mean(pvals_50v):.4f}")
print(f"  Median p-value:  {np.median(pvals_50v):.4f}")
print(f"  Range:           [{np.min(pvals_50v):.4f}, {np.max(pvals_50v):.4f}]")
print(f"  % below p<0.05:  {100*np.mean(pvals_50v < 0.05):.1f}%")
print(f"  % below p<0.10:  {100*np.mean(pvals_50v < 0.10):.1f}%")
print("  4-v-4 grid:      2/70, 4/70, 6/70, ...")
print("  Therefore p<0.05 means exactly p=2/70 for these configurations.")

# ---- Test B: mITT eteplirsen vs placebo ----
print(f"\n--- mITT Eteplirsen (n=6) vs Placebo (n=4), two-sided ---")
pvals_mitt = []
for i in range(min(1000, n_available)):
    mitt_trt = arm30_mitt_samples[i] + arm50_samples[i]  # 2 + 4 = 6
    pv = two_arm_perm_pvalue(mitt_trt, placebo_samples[i])
    pvals_mitt.append(pv)
    if (i+1) % 200 == 0:
        print(f"  ... computed {i+1} simulations")

pvals_mitt = np.array(pvals_mitt)
print(f"\n  Simulations run: {len(pvals_mitt)}")
print(f"  Mean p-value:    {np.mean(pvals_mitt):.4f}")
print(f"  Median p-value:  {np.median(pvals_mitt):.4f}")
print(f"  Range:           [{np.min(pvals_mitt):.4f}, {np.max(pvals_mitt):.4f}]")
print(f"  % below p<0.05:  {100*np.mean(pvals_mitt < 0.05):.1f}%")
print(f"  % below p<0.10:  {100*np.mean(pvals_mitt < 0.10):.1f}%")

# Save results
np.save(os.path.join(RESULTS_DIR, 'sensitivity_50v_pvals.npy'), pvals_50v)
np.save(os.path.join(RESULTS_DIR, 'sensitivity_mitt_pvals.npy'), pvals_mitt)

print(f"\n{'='*70}")
print("INTERPRETATION")
print(f"{'='*70}")
print(f"""
50mg vs placebo (two-sided):
  Below p<0.05 in {100*np.mean(pvals_50v < 0.05):.0f}% of generated configurations.
  Below p<0.10 in {100*np.mean(pvals_50v < 0.10):.0f}% of generated configurations.

mITT eteplirsen vs placebo (two-sided):
  Below p<0.05 in {100*np.mean(pvals_mitt < 0.05):.0f}% of generated configurations.
  Below p<0.10 in {100*np.mean(pvals_mitt < 0.10):.0f}% of generated configurations.

The difference between the selected group-mean targets is fixed regardless of
how individual values are distributed. The p-value varies because the
permutation distribution depends on the full set of values being permuted.

The individual-level ranges used here are ad hoc but intentionally wide. The
group means are selected reconstruction targets, not a complete identification
of the missing outcome vector. These frequencies describe the generator and
must not be interpreted as probabilities about the original patient data.

The mITT calculation follows post-randomisation exclusion and is exploratory;
this sensitivity analysis does not resolve that selection problem.

Outputs saved to {RESULTS_DIR}/
""")
