"""Independent Python validation of the imatinib shrinkage analysis."""

import csv
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.special import expit, logit
from scipy.stats import beta

SCRIPT_DIR = Path(__file__).resolve().parent

subtypes = [
    "Angiosarcoma", "Ewing", "Fibrosarcoma", "Leiomyosarcoma",
    "Liposarcoma", "MFH", "Osteosarcoma", "MPNST",
    "Rhabdomyosarcoma", "Synovial",
]
y = np.array([2, 0, 1, 6, 7, 3, 5, 1, 0, 3], dtype=float)
n = np.array([15, 13, 12, 28, 29, 29, 26, 5, 2, 20], dtype=float)
J = len(y)
l0 = logit(0.30)

alpha = 0.05
strat_pt = y / n
cp_lo = beta.ppf(alpha / 2, y, n - y + 1)
cp_lo[y == 0] = 0.0
cp_hi = beta.ppf(1 - alpha / 2, y + 1, n - y)
cp_hi[y == n] = 1.0

N_CHAINS = 4
N_ITER = 40_000
BURN = 10_000
A_SIGMA = 3.0
STEP_GAMMA = 0.6
STEP_LOG_SIGMA = 0.4
CHAIN_SEEDS = [20260707 + i for i in range(N_CHAINS)]


def loglik(gamma):
    p = np.clip(expit(l0 + gamma), 1e-12, 1 - 1e-12)
    return y * np.log(p) + (n - y) * np.log(1 - p)


def run_chain(seed):
    rng = np.random.default_rng(seed)
    gamma = np.zeros(J)
    mu = 0.0
    sigma = 1.0
    ll = loglik(gamma)
    keep_pi = np.empty((N_ITER - BURN, J))
    keep_sigma = np.empty(N_ITER - BURN)
    gamma_accept = 0
    sigma_accept = 0
    kept = 0

    for iteration in range(N_ITER):
        proposal = gamma + rng.normal(0, STEP_GAMMA, J)
        ll_proposal = loglik(proposal)
        lp_old = -0.5 * ((gamma - mu) / sigma) ** 2
        lp_new = -0.5 * ((proposal - mu) / sigma) ** 2
        accept = np.log(rng.random(J)) < (ll_proposal + lp_new) - (ll + lp_old)
        gamma_accept += int(accept.sum())
        gamma = np.where(accept, proposal, gamma)
        ll = np.where(accept, ll_proposal, ll)

        precision = 1 / 1e4 + J / sigma**2
        mu = rng.normal((gamma.sum() / sigma**2) / precision, np.sqrt(1 / precision))

        log_sigma = np.log(sigma)
        proposal_log_sigma = log_sigma + rng.normal(0, STEP_LOG_SIGMA)
        proposal_sigma = np.exp(proposal_log_sigma)

        def log_target(value, log_value):
            return (
                -J * np.log(value)
                - 0.5 * np.sum((gamma - mu) ** 2) / value**2
                - value**2 / (2 * A_SIGMA**2)
                + log_value
            )

        accept_sigma = np.log(rng.random()) < (
            log_target(proposal_sigma, proposal_log_sigma)
            - log_target(sigma, log_sigma)
        )
        if accept_sigma:
            sigma = proposal_sigma
            sigma_accept += 1

        if iteration >= BURN:
            keep_pi[kept] = expit(l0 + gamma)
            keep_sigma[kept] = sigma
            kept += 1

    return {
        "pi": keep_pi,
        "sigma": keep_sigma,
        "gamma_accept": gamma_accept / (N_ITER * J),
        "sigma_accept": sigma_accept / N_ITER,
    }


def gelman_rhat(chain_values):
    values = np.asarray(chain_values)
    n_draws = values.shape[1]
    within = np.var(values, axis=1, ddof=1).mean()
    between = n_draws * np.var(values.mean(axis=1), ddof=1)
    variance_plus = ((n_draws - 1) / n_draws) * within + between / n_draws
    return float(np.sqrt(variance_plus / within))


chains = [run_chain(seed) for seed in CHAIN_SEEDS]
draws = np.concatenate([chain["pi"] for chain in chains], axis=0)
rhat_pi = np.array([
    gelman_rhat([chain["pi"][:, j] for chain in chains]) for j in range(J)
])
rhat_sigma = gelman_rhat([chain["sigma"] for chain in chains])

post_mean = draws.mean(axis=0)
post_lo = np.percentile(draws, 2.5, axis=0)
post_hi = np.percentile(draws, 97.5, axis=0)
w_strat = np.mean(cp_hi - cp_lo)
w_borrow = np.mean(post_hi - post_lo)
width_reduction = 1 - w_borrow / w_strat

print(f"Overall observed rate: {y.sum() / n.sum() * 100:.1f}%   (paper: 15.6%)")
print(f"{'Subtype':<16}{'y/n':>7}{'strat%':>8}{'CP 95% CI':>14} | {'borrow%':>8}{'95% CrI':>14}")
for j in range(J):
    strat_ci = f"[{cp_lo[j] * 100:.0f},{cp_hi[j] * 100:.0f}]"
    borrow_ci = f"[{post_lo[j] * 100:.0f},{post_hi[j] * 100:.0f}]"
    print(
        f"{subtypes[j]:<16}{int(y[j]):>3}/{int(n[j]):<3}{strat_pt[j] * 100:>7.1f}"
        f"{strat_ci:>14} | {post_mean[j] * 100:>7.1f}{borrow_ci:>14}"
    )
print(
    f"\nMean 95% interval width  stratified: {w_strat * 100:.1f} pts   "
    f"borrowing: {w_borrow * 100:.1f} pts"
)
print(f"Interval shrinkage: {width_reduction * 100:.0f}% narrower on average")
print(
    f"Four-chain diagnostics: max R-hat {max(rhat_pi.max(), rhat_sigma):.4f}; "
    f"gamma acceptance {np.mean([c['gamma_accept'] for c in chains]) * 100:.1f}%; "
    f"sigma acceptance {np.mean([c['sigma_accept'] for c in chains]) * 100:.1f}%"
)

with (SCRIPT_DIR / "imatinib_validation.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow([
        "subtype", "y", "n", "stratified_rate", "stratified_ci_low",
        "stratified_ci_high", "posterior_mean", "posterior_ci_low",
        "posterior_ci_high", "rhat",
    ])
    for j in range(J):
        writer.writerow([
            subtypes[j], int(y[j]), int(n[j]), strat_pt[j], cp_lo[j], cp_hi[j],
            post_mean[j], post_lo[j], post_hi[j], rhat_pi[j],
        ])

diagnostics = {
    "implementation": "Python",
    "chains": N_CHAINS,
    "iterations_per_chain": N_ITER,
    "burn_in_per_chain": BURN,
    "retained_draws": int(draws.shape[0]),
    "mean_gamma_acceptance": float(np.mean([c["gamma_accept"] for c in chains])),
    "mean_sigma_acceptance": float(np.mean([c["sigma_accept"] for c in chains])),
    "max_response_rate_rhat": float(rhat_pi.max()),
    "sigma_rhat": rhat_sigma,
    "interval_width_reduction": float(width_reduction),
    "python_version": platform.python_version(),
    "numpy_version": np.__version__,
    "scipy_version": scipy.__version__,
}
(SCRIPT_DIR / "imatinib_validation_diagnostics.json").write_text(
    json.dumps(diagnostics, indent=2) + "\n"
)
