#!/usr/bin/env python3
"""
Shrinkage table supplement for:
"The P-Value Got a Lecture. It Needed a Lookup Table."
Evidence in the Wild, by Lu Qian

Implements van Zwet's conditional() function with mixture
parameters from Table 1 of van Zwet, Schwab & Senn (2021),
Statistics in Medicine, 40:6107-6117.

Reproduces the 7-trial shrinkage table from the essay.

Requirements: numpy, scipy
"""

import numpy as np
from scipy.stats import norm


# ---- Van Zwet's conditional() function ----
# Appendix of van Zwet, Schwab & Senn (2021), Stat Med 40:6107-6117.
#
# Computes the conditional distribution of SNR given the observed
# z-value, when the marginal SNR distribution is a mixture of
# zero-mean normals with proportions p and SDs tau.
#
# Returns (q, m, sigma):
#   q     - conditional mixing proportions
#   m     - conditional component means
#   sigma - conditional component standard deviations

def conditional(z, p, tau):
    tau2 = np.array(tau) ** 2
    q = np.array(p) * norm.pdf(z, 0, np.sqrt(tau2 + 1))
    q = q / q.sum()
    m = z * tau2 / (tau2 + 1)
    v = tau2 / (tau2 + 1)
    sigma = np.sqrt(v)
    return q, m, sigma


# Table 1 mixture parameters (van Zwet, Schwab & Senn 2021)
# Estimated from 23,551 z-statistics in the Cochrane Database.
MIX_P = [0.32, 0.31, 0.30, 0.07]  # proportions
MIX_TAU = [0.61, 1.42, 2.16, 5.64]  # SNR component SDs


# ---- Derived quantities ----

def shrinkage_estimate(z, s):
    """Shrinkage estimator: beta_hat = s * E(SNR | z)."""
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    return s * np.sum(q * m)


def conditional_power(z, n_mc=2_000_000):
    """E[power(SNR) | z] via Monte Carlo.

    Power at alpha = 0.05 (two-sided):
      power(snr) = Phi(-1.96 - snr) + 1 - Phi(1.96 - snr)
    """
    np.random.seed(42)
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    comp = np.random.choice(len(q), size=n_mc, p=q)
    snr = np.random.normal(m[comp], sigma[comp])
    pw = norm.cdf(-1.96 - snr) + 1 - norm.cdf(1.96 - snr)
    return pw.mean()


def p_direction_correct(z):
    """P(sign(SNR) = sign(z) | z)."""
    q, m, sigma = conditional(z, MIX_P, MIX_TAU)
    if z >= 0:
        return np.sum(q * norm.cdf(m / sigma))
    else:
        return np.sum(q * norm.cdf(-m / sigma))


def shrinkage_interval(z, s, level=0.95):
    """Conditional interval for beta, scaled from SNR quantiles."""
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


# ---- Trial data ----
# (name, type, estimate, CI_lower, CI_upper)
#
# b = log(estimate)
# s = (log(CI_upper) - log(CI_lower)) / 3.92
# z = b / s
#
# Sources:
#   ANDROMEDA-SHOCK:    Hernandez et al. 2019, JAMA 321(7):654-664
#   RECOVERY overall:   RECOVERY Collaborative Group 2021, NEJM 384:693-704
#   RECOVERY ventilated: ibid., invasive mechanical ventilation subgroup
#   RECOVERY no O2:     ibid., no supplemental oxygen subgroup
#   Ronco:              Ronco et al. 2000, Lancet 356:26-30
#   TOPCAT Americas:    Pfeffer et al. 2015, Circulation 131(1):34-42
#   CITRIS-ALI:         Fowler et al. 2019, JAMA 322(13):1261-1270
#
# NOTE: The RECOVERY subgroup CIs (ventilated: 0.48-0.88; no O2: 0.86-1.74)
# and point estimates (0.65; 1.22) differ slightly from the published
# age-adjusted rate ratios in the NEJM paper (0.64, 0.51-0.81; 1.19, 0.91-1.55).
# The values here reproduce the essay table; verify against your source extraction.

TRIALS = [
    ("ANDROMEDA-SHOCK",     "HR", 0.75, 0.55, 1.02),
    ("RECOVERY overall",    "RR", 0.83, 0.75, 0.93),
    ("RECOVERY ventilated", "RR", 0.65, 0.48, 0.88),
    ("RECOVERY no O2",      "RR", 1.22, 0.86, 1.74),
    ("Ronco hemofiltration", "OR", 1.91, 1.13, 3.24),
    ("TOPCAT Americas",     "HR", 0.82, 0.69, 0.98),
    ("CITRIS-ALI mortality", "OR", 0.49, 0.26, 0.93),
]


# ---- Main ----

if __name__ == "__main__":

    print()
    print("=" * 80)
    print("Shrinkage Table")
    print("Evidence in the Wild: The P-Value Got a Lecture. It Needed a Lookup Table.")
    print("=" * 80)
    print()
    print(f"{'Trial':<28} {'Original':>8} {'Shrunk':>8} {'Shrink':>7} "
          f"{'Power':>7} {'P(dir)':>7}")
    print("-" * 72)

    for name, typ, est, lo, hi in TRIALS:
        b = np.log(est)
        s = (np.log(hi) - np.log(lo)) / 3.92
        z = b / s

        beta_hat = shrinkage_estimate(z, s)
        est_shrunk = np.exp(beta_hat)
        shrink_pct = (1 - abs(beta_hat) / abs(b)) * 100
        pwr = conditional_power(z) * 100
        p_dir = p_direction_correct(z) * 100
        ci_lo, ci_hi = shrinkage_interval(z, s)

        label = f"{name} ({typ})"
        p_str = ">99%" if p_dir > 99 else f"{p_dir:.0f}%"

        print(f"{label:<28} {est:8.2f} {est_shrunk:8.2f} "
              f"{shrink_pct:6.0f}% {pwr:6.0f}% {p_str:>7}")

    print()
    print("All calculations use van Zwet's conditional() function with mixture")
    print("parameters from Table 1 of van Zwet, Schwab & Senn (2021),")
    print("Statistics in Medicine, 40:6107-6117.")

    # ---- Detailed output ----

    print()
    print("=" * 80)
    print("Detailed calculations")
    print("=" * 80)

    for name, typ, est, lo, hi in TRIALS:
        b = np.log(est)
        s = (np.log(hi) - np.log(lo)) / 3.92
        z = b / s

        beta_hat = shrinkage_estimate(z, s)
        est_shrunk = np.exp(beta_hat)
        shrink_pct = (1 - abs(beta_hat) / abs(b)) * 100
        pwr = conditional_power(z) * 100
        p_dir = p_direction_correct(z) * 100
        ci_lo, ci_hi = shrinkage_interval(z, s)

        print(f"\n{name} ({typ})")
        print(f"  Original: {typ} = {est:.2f} (95% CI: {lo:.2f} to {hi:.2f})")
        print(f"  b = log({est:.2f}) = {b:.4f}")
        print(f"  s = (log({hi:.2f}) - log({lo:.2f})) / 3.92 = {s:.4f}")
        print(f"  z = b/s = {z:.2f}")
        print(f"  Shrinkage estimate: {typ} = {est_shrunk:.2f}"
              f"  (beta_hat = {beta_hat:.4f})")
        print(f"  Shrinkage: {shrink_pct:.0f}%")
        print(f"  Conditional 95% interval: "
              f"({np.exp(ci_lo):.2f}, {np.exp(ci_hi):.2f})")
        print(f"  Conditional power: {pwr:.0f}%")
        print(f"  P(direction correct): {p_dir:.0f}%")

    print()
