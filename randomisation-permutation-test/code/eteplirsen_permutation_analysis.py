"""
Eteplirsen Study 201 (NCT01396239) — Permutation Test Analysis
==============================================================
"The Permutation Test Nobody Ran"

This script constructs one individual-level 6-Minute Walk Test (6MWT) outcome
configuration from selected published constraints for the eteplirsen Phase IIb
trial (Mendell et al. 2013), then exhaustively enumerates permutation reference
distributions. The input records are not original or recovered patient data.

Data sources:
- Mendell et al. (2013) Annals of Neurology 74:637-647 (primary publication)
- Sarepta Therapeutics Advisory Committee Presentation, April 25, 2016
- FDA Clinical Review (NDA 206488)

Methodological note:
- Horn (1983), "Some Easy t Statistics," motivates the discussion of
  alternative small-sample statistics. Its one-sample critical values are not
  used as randomisation critical values in this analysis.

Usage:
    Invoke this script by path from any working directory. Inputs and outputs
    are resolved relative to the script location. From the project root:
        python code/eteplirsen_permutation_analysis.py

Author: Evidence in the Wild (evidenceinthewild.com)
Collaborators: Robert Matthews, Maggie
"""

import numpy as np
from itertools import combinations
import json
import os
from datetime import datetime

# =============================================================================
# PATH SETUP — all outputs go to sibling directories
# =============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
FIGURES_DIR = os.path.join(PROJECT_ROOT, 'figures')
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(os.path.join(FIGURES_DIR, 'png'), exist_ok=True)
os.makedirs(os.path.join(FIGURES_DIR, 'pdf'), exist_ok=True)

# =============================================================================
# RECONSTRUCTED INDIVIDUAL PATIENT DATA
# =============================================================================
# Baseline values: integer selections matching published group means, SDs,
# minima, and maxima from Mendell et al. 2013 Table 1. A reproducible search
# over rounding tolerances is not implemented here, so do not call these
# solutions unique. See ../PROVENANCE.md.
#
# All three arms match their published baseline summary statistics to the
# displayed precision. That does not make the constructed patient records
# observed data.
#
# Patient 010 baseline (261m) is stated in the paper text.
# Patient 009 baseline (346m) is stated in the paper text.
#
# Patient IDs by arm (from Figure 3 biopsy images):
#   30 mg/kg: 02, 06, 09, 10
#   50 mg/kg: 03, 04, 12, 15
#   Placebo:  05, 07, 08, 13

patients = {
    # --- 30 mg/kg arm (n=4) ---
    # Baselines: {261, 346, 372, 442} → mean=355.25, SD=74.78, min=261, max=442
    # Published Table 1: mean=355.2, SD=74.78, min=261, max=442 ✓
    '010': {'arm': '30mg', 'baseline': 261, 'week24': 48,  'change': -213, 'source': 'baseline and Week 24 decline reported in paper; Week 24 value derived'},
    '009': {'arm': '30mg', 'baseline': 346, 'week24': 52,  'change': -294, 'source': 'baseline reported; Week 24 value reconstructed using narrative and combined CE-23 trajectory'},
    '02A': {'arm': '30mg', 'baseline': 372, 'week24': 397, 'change': 25,   'source': 'reconstructed (constrained by mITT group mean change)'},
    '06A': {'arm': '30mg', 'baseline': 442, 'week24': 457, 'change': 15,   'source': 'reconstructed (constrained by mITT group mean change)'},

    # --- 50 mg/kg arm (n=4) ---
    # Baselines: constructed integer values matching published summaries
    '03A': {'arm': '50mg', 'baseline': 365, 'week24': 360, 'change': -5,   'source': 'reconstructed (constrained by group mean change)'},
    '04A': {'arm': '50mg', 'baseline': 389, 'week24': 409, 'change': 20,   'source': 'reconstructed (constrained by group mean change)'},
    '12A': {'arm': '50mg', 'baseline': 401, 'week24': 426, 'change': 25,   'source': 'reconstructed (constrained by group mean change)'},
    '15A': {'arm': '50mg', 'baseline': 429, 'week24': 439, 'change': 10,   'source': 'reconstructed (constrained by group mean change)'},

    # --- Placebo arm (n=4) ---
    # Baselines: constructed integer values matching published summaries
    '05P': {'arm': 'placebo', 'baseline': 364, 'week24': 334, 'change': -30,  'source': 'reconstructed (constrained by group mean change)'},
    '07P': {'arm': 'placebo', 'baseline': 370, 'week24': 330, 'change': -40,  'source': 'reconstructed (constrained by group mean change)'},
    '08P': {'arm': 'placebo', 'baseline': 388, 'week24': 338, 'change': -50,  'source': 'reconstructed (constrained by group mean change)'},
    '13P': {'arm': 'placebo', 'baseline': 456, 'week24': 356, 'change': -100, 'source': 'reconstructed (constrained by group mean change)'},
}

# Exploratory post-randomisation subset used in the mITT calculation. Encode
# the clinical exclusion explicitly rather than inferring it from an outcome.
MITT_EXCLUDED_PATIENT_IDS = {'009', '010'}

# =============================================================================
# VERIFICATION: Check reconstructed data matches published statistics
# =============================================================================
def verify_data():
    """Verify reconstructed data matches published summary statistics."""
    arms = {'30mg': [], '50mg': [], 'placebo': []}
    for p in patients.values():
        arms[p['arm']].append(p)

    print("=" * 70)
    print("VERIFICATION: Reconstructed vs Published Summary Statistics")
    print("=" * 70)

    # Published baseline stats from Mendell 2013 Table 1
    published = {
        '30mg':    {'mean_bl': 355.2, 'sd_bl': 74.8, 'min_bl': 261, 'max_bl': 442, 'n': 4},
        '50mg':    {'mean_bl': 396.0, 'sd_bl': 26.6, 'min_bl': 365, 'max_bl': 429, 'n': 4},
        'placebo': {'mean_bl': 394.5, 'sd_bl': 42.0, 'min_bl': 364, 'max_bl': 456, 'n': 4},
    }

    for arm_name in ['30mg', '50mg', 'placebo']:
        bl = [p['baseline'] for p in arms[arm_name]]
        ch = [p['change'] for p in arms[arm_name]]
        pub = published[arm_name]

        calc_mean = np.mean(bl)
        calc_sd = np.std(bl, ddof=1)

        mean_match = "✓" if abs(calc_mean - pub['mean_bl']) < 0.15 else "✗"
        sd_match = "✓" if abs(calc_sd - pub['sd_bl']) < 0.5 else "≈"
        range_match = "✓" if min(bl) == pub['min_bl'] and max(bl) == pub['max_bl'] else "✗"

        print(f"\n{arm_name} arm (n={len(bl)}):")
        print(f"  Baseline mean:  {calc_mean:.1f}  (published: {pub['mean_bl']}) {mean_match}")
        print(f"  Baseline SD:    {calc_sd:.1f}  (published: {pub['sd_bl']}) {sd_match}")
        print(f"  Baseline range: {min(bl)}-{max(bl)} (published: {pub['min_bl']}-{pub['max_bl']}) {range_match}")
        print(f"  Mean change:    {np.mean(ch):.1f}m")

    # Selected Week 24 group-mean targets. Figure 6 shows a combined six-patient
    # ambulation-evaluable trajectory, not a separate 50 mg/kg raw mean. The
    # 50 mg/kg value used by this reconstruction requires independent source
    # verification and must not be described as published by Figure 6.
    placebo_mean = np.mean([p['change'] for p in arms['placebo']])
    eteplirsen_mitt = [p['change'] for patient_id, p in patients.items()
                       if p['arm'] == '50mg'
                       or (p['arm'] == '30mg'
                           and patient_id not in MITT_EXCLUDED_PATIENT_IDS)]
    mitt_mean = np.mean(eteplirsen_mitt)

    print(f"\nSelected Week 24 group-mean targets used in this reconstruction:")
    print(f"  Placebo mean change:         {placebo_mean:.1f}m (approximate graphical target)")
    print(f"  mITT eteplirsen mean change: {mitt_mean:.1f}m  (approximate combined-cohort target)")
    print("  See PROVENANCE.md; individual outcomes and the 50 mg/kg raw mean are not tabulated in Figure 6.")

# =============================================================================
# PERMUTATION TEST FUNCTIONS
# =============================================================================
def two_arm_permutation_test(treatment_changes, control_changes, n_treatment=None):
    """
    Exact two-arm permutation test.

    Tests H0: no treatment effect by enumerating all possible assignments
    of the combined data into two groups of the original sizes.

    Returns: observed difference, one-sided p (trt > ctrl),
             two-sided p (|diff| >= |observed|), full permutation distribution
    """
    all_changes = np.array(list(treatment_changes) + list(control_changes))
    n_total = len(all_changes)
    if n_treatment is None:
        n_treatment = len(treatment_changes)

    observed_diff = np.mean(treatment_changes) - np.mean(control_changes)

    # Enumerate all C(n_total, n_treatment) permutations
    perm_diffs = []
    for idx in combinations(range(n_total), n_treatment):
        trt = all_changes[list(idx)]
        ctrl = all_changes[[i for i in range(n_total) if i not in idx]]
        perm_diffs.append(np.mean(trt) - np.mean(ctrl))

    perm_diffs = np.array(perm_diffs)

    # One-sided p-value: proportion of permutations >= observed
    p_one = np.mean(perm_diffs >= observed_diff)
    # Two-sided p-value: proportion with |diff| >= |observed|
    p_two = np.mean(np.abs(perm_diffs) >= abs(observed_diff))

    return observed_diff, p_one, p_two, perm_diffs

def two_arm_permutation_welch_t(treatment_changes, control_changes, n_treatment=None):
    """
    Two-arm randomisation test using the Welch t-statistic.

    This uses the same assignment enumeration as the difference-in-means test,
    but studentizes the contrast with the arm-specific sample variances. Under
    the sharp null used here, exactness comes from the assignment mechanism and
    does not require equal population variances. Studentization is included as
    an alternative statistic; it can be useful for weak-null interpretations
    and heterogeneous outcomes, but this finite enumeration is not
    automatically exact for a weak null of zero average treatment effect.

    With n=4 per arm and only 70 permutations, the entire distribution is
    enumerated—no Monte Carlo approximation is involved.

    Returns: observed Welch t, one-sided p, two-sided p, permutation distribution
    """
    all_changes = np.array(list(treatment_changes) + list(control_changes))
    n_total = len(all_changes)
    if n_treatment is None:
        n_treatment = len(treatment_changes)
    n_control = n_total - n_treatment

    def welch_t(trt, ctrl):
        n1, n2 = len(trt), len(ctrl)
        var1 = np.var(trt, ddof=1) if n1 > 1 else 0.0
        var2 = np.var(ctrl, ddof=1) if n2 > 1 else 0.0
        denom = np.sqrt(var1 / n1 + var2 / n2)
        if denom == 0:
            return 0.0
        return (np.mean(trt) - np.mean(ctrl)) / denom

    observed_t = welch_t(np.array(treatment_changes), np.array(control_changes))

    perm_ts = []
    for idx in combinations(range(n_total), n_treatment):
        trt = all_changes[list(idx)]
        ctrl = all_changes[[i for i in range(n_total) if i not in idx]]
        perm_ts.append(welch_t(trt, ctrl))

    perm_ts = np.array(perm_ts)

    p_one = np.mean(perm_ts >= observed_t)
    p_two = np.mean(np.abs(perm_ts) >= abs(observed_t))

    return observed_t, p_one, p_two, perm_ts


def three_arm_permutation_test(arm1, arm2, arm3):
    """
    Exact three-arm permutation test using F-statistic analog.

    Enumerates all C(12,4) x C(8,4) = 34,650 possible allocations.
    The F-test is inherently non-directional, so only one p-value is returned.
    """
    all_changes = np.array(list(arm1) + list(arm2) + list(arm3))
    n1, n2, n3 = len(arm1), len(arm2), len(arm3)
    n_total = n1 + n2 + n3

    def f_stat(a1, a2, a3):
        gm = np.mean(np.concatenate([a1, a2, a3]))
        ss_between = (len(a1) * (np.mean(a1) - gm)**2 +
                      len(a2) * (np.mean(a2) - gm)**2 +
                      len(a3) * (np.mean(a3) - gm)**2)
        ss_within = (np.sum((a1 - np.mean(a1))**2) +
                     np.sum((a2 - np.mean(a2))**2) +
                     np.sum((a3 - np.mean(a3))**2))
        if ss_within == 0:
            return float('inf')
        k = 3
        return (ss_between / (k - 1)) / (ss_within / (n_total - k))

    observed_f = f_stat(np.array(arm1), np.array(arm2), np.array(arm3))

    count_ge = 0
    total = 0
    indices = list(range(n_total))

    for idx1 in combinations(indices, n1):
        remaining = [i for i in indices if i not in idx1]
        for idx2 in combinations(remaining, n2):
            idx3 = [i for i in remaining if i not in idx2]
            a1 = all_changes[list(idx1)]
            a2 = all_changes[list(idx2)]
            a3 = all_changes[list(idx3)]
            f = f_stat(a1, a2, a3)
            if f >= observed_f:
                count_ge += 1
            total += 1

    p_value = count_ge / total
    return observed_f, p_value, total

# =============================================================================
# MAIN ANALYSIS
# =============================================================================
def run_analysis():
    verify_data()

    # Extract changes by arm
    arm_30 = [p['change'] for p in patients.values() if p['arm'] == '30mg']
    arm_50 = [p['change'] for p in patients.values() if p['arm'] == '50mg']
    placebo = [p['change'] for p in patients.values() if p['arm'] == 'placebo']

    # Exploratory mITT: exclude the two participants who lost ambulation.
    arm_30_mitt = [p['change'] for patient_id, p in patients.items()
                   if p['arm'] == '30mg'
                   and patient_id not in MITT_EXCLUDED_PATIENT_IDS]

    print("\n" + "=" * 70)
    print("PERMUTATION TEST RESULTS (one-sided and two-sided)")
    print("=" * 70)

    results = {}

    # --- Test 1: ITT all eteplirsen vs placebo ---
    all_eteplirsen = arm_30 + arm_50
    obs, p1, p2, perm_dist = two_arm_permutation_test(all_eteplirsen, placebo, n_treatment=8)
    results['itt_all_vs_placebo'] = {'observed': obs, 'p_one': p1, 'p_two': p2, 'n_perms': len(perm_dist)}
    print(f"\n1. ITT: All eteplirsen (n=8) vs placebo (n=4)")
    print(f"   Displayed Δ:         {obs:+.1f}m")
    print(f"   Permutation p (1-sided): {p1:.4f}")
    print(f"   Permutation p (2-sided): {p2:.4f}")
    print(f"   Total permutations:  C(12,8) = {len(perm_dist)}")

    # --- Test 2: mITT eteplirsen vs placebo ---
    mitt_eteplirsen = arm_30_mitt + arm_50
    obs2, p1_2, p2_2, perm_dist2 = two_arm_permutation_test(mitt_eteplirsen, placebo, n_treatment=6)
    results['mitt_vs_placebo'] = {'observed': obs2, 'p_one': p1_2, 'p_two': p2_2, 'n_perms': len(perm_dist2)}
    print(f"\n2. EXPLORATORY mITT: Eteplirsen (n=6, post-randomisation exclusion) vs placebo (n=4)")
    print(f"   Displayed Δ:         {obs2:+.1f}m")
    print(f"   Permutation p (1-sided): {p1_2:.4f}")
    print(f"   Permutation p (2-sided): {p2_2:.4f}")
    print(f"   Total permutations:  C(10,6) = {len(perm_dist2)}")

    # --- Test 3: 50 mg/kg vs placebo (displayed construction) ---
    obs3, p1_3, p2_3, perm_dist3 = two_arm_permutation_test(arm_50, placebo)
    results['50mg_vs_placebo'] = {'observed': obs3, 'p_one': p1_3, 'p_two': p2_3, 'n_perms': len(perm_dist3)}
    print(f"\n3. ★ 50 mg/kg (n=4) vs placebo (n=4) — CONSTRUCTION-CONDITIONAL COMPARISON")
    print(f"   Displayed Δ:         {obs3:+.1f}m")
    print(f"   Permutation p (1-sided): {p1_3:.4f}")
    print(f"   Permutation p (2-sided): {p2_3:.4f}")
    print(f"   Approx. MMRM p (2-sided, derived): 0.56")
    print(f"   Total permutations:  C(8,4) = {len(perm_dist3)}")

    # --- Test 4: 30 mg/kg (ITT) vs placebo ---
    obs4, p1_4, p2_4, perm_dist4 = two_arm_permutation_test(arm_30, placebo)
    results['30mg_vs_placebo'] = {'observed': obs4, 'p_one': p1_4, 'p_two': p2_4, 'n_perms': len(perm_dist4)}
    print(f"\n4. 30 mg/kg ITT (n=4) vs placebo (n=4)")
    print(f"   Displayed Δ:         {obs4:+.1f}m")
    print(f"   Permutation p (1-sided): {p1_4:.4f}")
    print(f"   Permutation p (2-sided): {p2_4:.4f}")
    print(f"   Published MMRM p (2-sided): 0.026 (favoring placebo)")
    print(f"   Total permutations:  C(8,4) = {len(perm_dist4)}")

    # --- Test 5: Three-arm F-test (ITT) ---
    print(f"\n5. Three-arm F-test (ITT, all 12 patients)")
    obs_f, pval_f, total_perms = three_arm_permutation_test(arm_30, arm_50, placebo)
    results['three_arm_f'] = {'observed_f': obs_f, 'p_value': pval_f, 'n_perms': total_perms}
    print(f"   Displayed F:         {obs_f:.3f}")
    print(f"   Permutation p:       {pval_f:.4f}")
    print(f"   Total permutations:  C(12,4)×C(8,4) = {total_perms}")

    # --- Sensitivity to the choice of test statistic ---
    print("\n" + "=" * 70)
    print("SENSITIVITY CHECK: Welch-studentized randomisation statistic")
    print("(Alternative statistic; not an equal-variance repair under the sharp null)")
    print("=" * 70)

    obs_wt3, pw1_3, pw2_3, wt_dist3 = two_arm_permutation_welch_t(arm_50, placebo)
    results['50mg_vs_placebo_welch'] = {
        'observed_t': obs_wt3, 'p_one': pw1_3, 'p_two': pw2_3,
        'n_perms': len(wt_dist3),
    }
    print(f"\n6a. ★ 50 mg/kg vs placebo — Welch t")
    print(f"   Displayed Welch t:       {obs_wt3:+.3f}")
    print(f"   Permutation p (1-sided): {pw1_3:.4f}")
    print(f"   Permutation p (2-sided): {pw2_3:.4f}")
    print(f"   [cf. difference-in-means p (2-sided): {p2_3:.4f}]")

    obs_wt2, pw1_2m, pw2_2m, wt_dist2 = two_arm_permutation_welch_t(mitt_eteplirsen, placebo, n_treatment=6)
    results['mitt_vs_placebo_welch'] = {
        'observed_t': obs_wt2, 'p_one': pw1_2m, 'p_two': pw2_2m,
        'n_perms': len(wt_dist2),
    }
    print(f"\n6b. EXPLORATORY mITT eteplirsen vs placebo — Welch t")
    print(f"   Displayed Welch t:       {obs_wt2:+.3f}")
    print(f"   Permutation p (1-sided): {pw1_2m:.4f}")
    print(f"   Permutation p (2-sided): {pw2_2m:.4f}")
    print(f"   [cf. difference-in-means p (2-sided): {p2_2:.4f}]")

    print("\n   Interpretation: for these constructed outcomes and displayed splits,")
    print("   the Welch-studentized and mean-difference statistics give the same")
    print("   tail counts. That agreement does not validate the reconstructed")
    print("   outcomes, mITT selection, or assumed assignment set.")

    # =============================================================================
    # SUMMARY TABLE
    # =============================================================================
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'Comparison':<40} {'Disp Δ':>8} {'p(1-sided)':>11} {'p(2-sided)':>11} {'MMRM p':>10}")
    print("-" * 80)
    for label, key in [
        ('ITT: All eteplirsen vs placebo', 'itt_all_vs_placebo'),
        ('Exploratory mITT: Eteplirsen vs placebo', 'mitt_vs_placebo'),
        ('★ 50 mg/kg vs placebo', '50mg_vs_placebo'),
        ('30 mg/kg (ITT) vs placebo', '30mg_vs_placebo'),
    ]:
        r = results[key]
        mmrm = {'50mg_vs_placebo': '≈0.56†', '30mg_vs_placebo': '0.026*'}.get(key, 'n/a')
        print(f"{label:<40} {r['observed']:>+7.1f}m {r['p_one']:>11.4f} {r['p_two']:>11.4f} {mmrm:>10}")

    f_val = results['three_arm_f']['observed_f']
    print(f"{'3-arm F-test (ITT)':<40} {'F=%.2f' % f_val:>8} {'—':>11} {results['three_arm_f']['p_value']:>11.4f} {'n/a':>10}")
    print("\n* Published MMRM p=0.026 favors placebo over 30mg")
    print("† Approximate value derived from published Week 24 MMRM adjusted means and SEs;")
    print("  Mendell et al. report no significant difference but not this pairwise p-value.")
    print("  The MMRM and permutation calculations test different analysis questions.")
    print("  All permutation p-values are conditional on constructed outcomes; the mITT")
    print("  calculation also follows post-randomisation selection. See PROVENANCE.md.")

    # Save permutation distributions
    np.save(os.path.join(RESULTS_DIR, 'perm_dist_50v_placebo.npy'), perm_dist3)
    np.save(os.path.join(RESULTS_DIR, 'perm_dist_mitt.npy'), perm_dist2)

    # Save results as JSON
    results_json = {
        'analysis_date': datetime.now().isoformat(),
        'study': 'Eteplirsen Study 201 (NCT01396239)',
        'source': 'Mendell et al. 2013, Annals of Neurology 74:637-647',
        'patients': {k: v for k, v in patients.items()},
        'results': {k: {kk: float(vv) if isinstance(vv, (np.floating, float)) else vv
                        for kk, vv in v.items()} for k, v in results.items()},
    }
    with open(os.path.join(RESULTS_DIR, 'permutation_results.json'), 'w') as f:
        json.dump(results_json, f, indent=2, default=str)

    print(f"\nOutputs saved to {RESULTS_DIR}/")
    return results, perm_dist3, perm_dist2

if __name__ == '__main__':
    results, perm_dist_50, perm_dist_mitt = run_analysis()
