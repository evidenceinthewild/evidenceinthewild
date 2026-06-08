"""
Sensitivity Analysis for Eteplirsen Permutation Tests
=====================================================
Tests robustness of the permutation p-value across many plausible
individual patient allocations that still match published group means.

The key question: does the permutation result depend on how we
distributed individual changes within each group, or is it robust?

Usage:
    Run from the project root (the parent of code/, results/, figures/):
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

# Fixed constraints
# Patient 010: change = -213 (exact from paper)
# Patient 009: change = -294 (derived from CE-23; baseline 346m from paper)
# 50mg group mean change = 12.5m (from Figure 6)
# Placebo group mean change = -55.0m (from Figure 6)
# mITT eteplirsen mean change = +15.0m (2 ambulant 30mg + 4 50mg) / 6

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

def generate_plausible_changes(group_mean, n, constraints=None, n_samples=5000):
    """
    Generate plausible individual change sets that sum to n * group_mean.

    For DMD patients in a 24-week study, individual changes are constrained
    to physiologically plausible ranges. Ranges are intentionally wide to
    avoid biasing toward significance; see note at end of output.
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

        # Check all values are in plausible range
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
print("SENSITIVITY ANALYSIS (two-sided p-values)")
print("Varying individual patient changes while preserving group means")
print("=" * 70)

# Generate many plausible individual allocations
n_sims = 5000
placebo_samples = generate_plausible_changes(-55.0, 4, 'placebo', n_sims)
arm50_samples = generate_plausible_changes(12.5, 4, '50mg', n_sims)
arm30_mitt_samples = generate_plausible_changes(20.0, 2, '30mg_mitt', n_sims)
# mITT 30mg ambulant mean = (25 + 15) / 2 = 20.0 (matching published ~+15m for full mITT)

n_available = min(len(placebo_samples), len(arm50_samples), len(arm30_mitt_samples))
print(f"\nGenerated {n_available} plausible patient allocations per arm")

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
  Below p<0.05 in {100*np.mean(pvals_50v < 0.05):.0f}% of plausible allocations.
  Below p<0.10 in {100*np.mean(pvals_50v < 0.10):.0f}% of plausible allocations.

mITT eteplirsen vs placebo (two-sided):
  Below p<0.05 in {100*np.mean(pvals_mitt < 0.05):.0f}% of plausible allocations.
  Below p<0.10 in {100*np.mean(pvals_mitt < 0.10):.0f}% of plausible allocations.

The observed difference in group means is fixed regardless of how individual
values are distributed. The p-value varies because the permutation distribution
depends on the full set of values being permuted.

The individual-level ranges used here are ad hoc but intentionally wide. The
primary constraint is the published group mean; within-group variability is
a secondary influence on the permutation p-value. Narrower physiological
constraints (e.g., using the published within-group SDs) would reduce
the spread of p-values.

Outputs saved to {RESULTS_DIR}/
""")
