#!/usr/bin/env python3
"""
FDA Approvals Shrinkage Analysis
Evidence in the Wild: Pivotal Trials Meet the Shrinkage Estimator

Uses van Zwet's conditional() function with mixture parameters from
Table 1 of van Zwet, Schwab & Senn (2021), Statistics in Medicine, 40:6107-6117.

Applies shrinkage to pivotal trials across therapeutic areas.
"""

import numpy as np
from scipy.stats import norm


# ---- Van Zwet's conditional() function ----
def conditional(z, p, tau):
    tau2 = np.array(tau) ** 2
    q = np.array(p) * norm.pdf(z, 0, np.sqrt(tau2 + 1))
    q = q / q.sum()
    m = z * tau2 / (tau2 + 1)
    v = tau2 / (tau2 + 1)
    sigma = np.sqrt(v)
    return q, m, sigma

# Table 1 mixture parameters
MIX_P = [0.32, 0.31, 0.30, 0.07]
MIX_TAU = [0.61, 1.42, 2.16, 5.64]


def shrinkage_estimate(z, s):
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    return s * np.sum(q * m)

def conditional_power(z, n_mc=2_000_000):
    np.random.seed(42)
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    comp = np.random.choice(len(q), size=n_mc, p=q)
    snr = np.random.normal(m[comp], sigma[comp])
    pw = norm.cdf(-1.96 - snr) + 1 - norm.cdf(1.96 - snr)
    return pw.mean()

def p_direction_correct(z):
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    if z >= 0:
        return np.sum(q * norm.cdf(m / sigma))
    else:
        return np.sum(q * norm.cdf(-m / sigma))

def shrinkage_interval(z, s, level=0.95):
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    alpha = 1 - level
    def mix_cdf(x):
        return np.sum(q * norm.cdf(x, m, sigma))
    def mix_quantile(p_target, lo=-20, hi=20, tol=1e-8):
        for _ in range(200):
            mid = (lo + hi) / 2
            if mix_cdf(mid) < p_target:
                lo = mid
            else:
                hi = mid
            if hi - lo < tol:
                break
        return (lo + hi) / 2
    snr_lo = mix_quantile(alpha / 2)
    snr_hi = mix_quantile(1 - alpha / 2)
    return s * snr_lo, s * snr_hi


# ---- FDA-Approved Pivotal Trials ----
# (name, therapeutic_area, type, estimate, CI_lower, CI_upper, endpoint, year, notes)
#
# Numbers are from published pivotal trial results.
# VERIFY AGAINST PUBLICATIONS BEFORE FINAL USE.

TRIALS = [
    # --- Oncology ---
    ("DESTINY-Breast03\n(trastuzumab deruxtecan)",
     "Oncology", "HR", 0.28, 0.22, 0.37, "PFS", 2022,
     "T-DXd vs T-DM1 in HER2+ breast cancer. NEJM 2022."),

    ("KEYNOTE-024\n(pembrolizumab)",
     "Oncology", "HR", 0.60, 0.41, 0.89, "OS", 2016,
     "Pembro vs chemo in NSCLC PD-L1>=50%. NEJM 2016; Reck et al."),

    ("FLAURA\n(osimertinib)",
     "Oncology", "HR", 0.80, 0.64, 1.00, "OS", 2020,
     "Osimertinib vs comparator EGFR-TKI in EGFR+ NSCLC. NEJM 2020. 95.05% CI."),

    ("ASCENT\n(sacituzumab govitecan)",
     "Oncology", "HR", 0.48, 0.38, 0.59, "OS", 2021,
     "SG vs chemo in mTNBC. NEJM 2021; Bardia et al. Primary efficacy population."),

    ("EV-302/KEYNOTE-A39\n(enfortumab vedotin + pembro)",
     "Oncology", "HR", 0.47, 0.38, 0.58, "OS", 2024,
     "EV+pembro vs chemo in advanced urothelial. NEJM 2024."),

    # --- Cardiology ---
    ("DAPA-HF\n(dapagliflozin)",
     "Cardiology", "HR", 0.74, 0.65, 0.85, "CV death/HF hosp", 2019,
     "Dapagliflozin vs placebo in HFrEF. NEJM 2019."),

    ("EMPEROR-Preserved\n(empagliflozin)",
     "Cardiology", "HR", 0.79, 0.69, 0.90, "CV death/HF hosp", 2021,
     "Empagliflozin vs placebo in HFpEF. NEJM 2021."),

    ("PARADIGM-HF\n(sacubitril/valsartan)",
     "Cardiology", "HR", 0.80, 0.73, 0.87, "CV death/HF hosp", 2014,
     "Sacubitril/valsartan vs enalapril in HFrEF. NEJM 2014."),

    ("VICTORIA\n(vericiguat)",
     "Cardiology", "HR", 0.90, 0.82, 0.98, "CV death/HF hosp", 2020,
     "Vericiguat vs placebo in worsening HF. NEJM 2020."),

    ("REDUCE-IT\n(icosapent ethyl)",
     "Cardiology", "HR", 0.75, 0.68, 0.83, "MACE", 2019,
     "Icosapent ethyl vs mineral oil placebo. NEJM 2019. Controversial control."),

    # --- Endocrine / CV ---
    ("SELECT\n(semaglutide 2.4mg)",
     "Endocrine/CV", "HR", 0.80, 0.72, 0.90, "MACE", 2023,
     "Semaglutide vs placebo in obesity with CVD. NEJM 2023."),

    # --- Nephrology ---
    ("DAPA-CKD\n(dapagliflozin)",
     "Nephrology", "HR", 0.61, 0.51, 0.72, "Kidney composite", 2020,
     "Dapagliflozin vs placebo in CKD. NEJM 2020."),

    # --- Hematology ---
    ("RESONATE-2\n(ibrutinib)",
     "Hematology", "HR", 0.16, 0.09, 0.28, "PFS", 2015,
     "Ibrutinib vs chlorambucil in CLL 1L. NEJM 2015."),

    # --- Infectious Disease ---
    ("EPIC-HR\n(nirmatrelvir/ritonavir)",
     "Infectious Disease", "RR", 0.12, 0.06, 0.25, "Hosp/death", 2022,
     "Paxlovid vs placebo in high-risk COVID. NEJM 2022. Final 5-day mITT: 8/1039 vs 66/1046."),

    # --- Callback to GBM AGILE essay ---
    ("REGOMA\n(regorafenib)",
     "Oncology", "HR", 0.50, 0.33, 0.75, "OS", 2019,
     "Known false positive. GBM AGILE showed HR 1.07. Lancet Oncol 2019."),
]


# ---- Compute ----

if __name__ == "__main__":

    print()
    print("=" * 95)
    print("What Would the Shrinkage Estimator Say About FDA Approvals?")
    print("Evidence in the Wild, by Lu Qian")
    print("=" * 95)
    print()

    # Summary table
    print(f"{'Trial':<32} {'Area':<16} {'Orig':>6} {'Shrunk':>6} "
          f"{'Shrink':>7} {'Power':>7} {'P(dir)':>7} {'z':>6}")
    print("-" * 95)

    results = []

    for entry in TRIALS:
        name, area, typ, est, lo, hi, endpoint, year, notes = entry

        b = np.log(est)
        s = (np.log(hi) - np.log(lo)) / 3.92
        z = b / s
        p_val = 2 * (1 - norm.cdf(abs(z)))

        beta_hat = shrinkage_estimate(z, s)
        est_shrunk = np.exp(beta_hat)
        shrink_pct = (1 - abs(beta_hat) / abs(b)) * 100
        pwr = conditional_power(z) * 100
        p_dir = p_direction_correct(z) * 100
        ci_lo, ci_hi = shrinkage_interval(z, s)

        # Clean up multiline name for table
        label = name.replace('\n', ' ')
        if len(label) > 30:
            label = label[:28] + ".."

        p_str = ">99%" if p_dir > 99 else f"{p_dir:.0f}%"

        print(f"{label:<32} {area:<16} {est:6.2f} {est_shrunk:6.2f} "
              f"{shrink_pct:6.0f}% {pwr:6.0f}% {p_str:>7} {z:6.2f}")

        results.append({
            'name': name, 'area': area, 'type': typ,
            'est': est, 'lo': lo, 'hi': hi, 'endpoint': endpoint,
            'year': year, 'notes': notes,
            'b': b, 's': s, 'z': z, 'p_val': p_val,
            'beta_hat': beta_hat, 'est_shrunk': est_shrunk,
            'shrink_pct': shrink_pct, 'pwr': pwr, 'p_dir': p_dir,
            'ci_lo_shrunk': np.exp(ci_lo), 'ci_hi_shrunk': np.exp(ci_hi)
        })

    # ---- By therapeutic area ----
    print()
    print("=" * 95)
    print("Summary by therapeutic area")
    print("=" * 95)

    areas = {}
    for r in results:
        a = r['area']
        if a not in areas:
            areas[a] = []
        areas[a].append(r)

    for area, trials in areas.items():
        shrinks = [t['shrink_pct'] for t in trials]
        powers = [t['pwr'] for t in trials]
        print(f"\n{area} (n={len(trials)})")
        print(f"  Median shrinkage: {np.median(shrinks):.0f}%")
        print(f"  Median conditional power: {np.median(powers):.0f}%")
        print(f"  Range of z-values: {min(t['z'] for t in trials):.2f} to {max(t['z'] for t in trials):.2f}")

    # ---- Detailed output ----
    print()
    print("=" * 95)
    print("Detailed calculations")
    print("=" * 95)

    for r in results:
        name = r['name'].replace('\n', ' ')
        print(f"\n{name} ({r['type']}) — {r['area']}")
        print(f"  {r['endpoint']}, {r['year']}")
        print(f"  Original: {r['type']} = {r['est']:.2f} (95% CI: {r['lo']:.2f} to {r['hi']:.2f})")
        print(f"  b = {r['b']:.4f}, s = {r['s']:.4f}, z = {r['z']:.2f}, p = {r['p_val']:.6f}")
        print(f"  Shrinkage estimate: {r['type']} = {r['est_shrunk']:.2f} (beta_hat = {r['beta_hat']:.4f})")
        print(f"  Shrinkage: {r['shrink_pct']:.0f}%")
        print(f"  Shrinkage 95% interval: ({r['ci_lo_shrunk']:.2f}, {r['ci_hi_shrunk']:.2f})")
        print(f"  Conditional power: {r['pwr']:.0f}%")
        print(f"  P(direction correct): {r['p_dir']:.0f}%")
        print(f"  Notes: {r['notes']}")

    print()
