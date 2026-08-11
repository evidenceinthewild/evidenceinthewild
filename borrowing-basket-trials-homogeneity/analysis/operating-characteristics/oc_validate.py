"""Independent Python validation of the basket-trial OC simulation."""

import csv
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit, logit

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PATH = SCRIPT_DIR / "oc_simulation_inputs.csv"

l0 = logit(0.20)
N = 20
NB = 4
R_CAL = 2_500
R_EVAL = 5_000
N_ITER = 8_000
BURN = 3_000

SCEN = {
    1: [.20, .20, .20, .20],
    2: [.35, .35, .35, .35],
    3: [.20, .35, .35, .35],
    4: [.20, .20, .35, .35],
    5: [.10, .20, .30, .40],
    6: [.20, .20, .20, .35],
}
PROMISING = {s: [j for j, p in enumerate(SCEN[s]) if p > 0.20] for s in SCEN}
NULLSET = {s: [j for j, p in enumerate(SCEN[s]) if p <= 0.20] for s in SCEN}
Q_PAPER = {
    "weak": {"No": 0.982, "Moderate": 0.946, "Strong": 0.964},
    "strong": {"No": 0.982, "Moderate": 0.996, "Strong": 0.984},
}
A_SIGMA = {"Moderate": 3.0, "Strong": 0.3}


def binom_ll(y, gamma):
    p = np.clip(expit(l0 + gamma), 1e-12, 1 - 1e-12)
    return y * np.log(p) + (N - y) * np.log(1 - p)


def probpos_noborrow_table():
    """Pr(gamma > 0 | y) under N(0,100^2), integrated on (-inf, inf)."""

    def kernel(gamma, y_value):
        p = float(expit(l0 + gamma))
        prior = math.exp(-0.5 * (gamma / 100) ** 2) / (100 * math.sqrt(2 * math.pi))
        return prior * p**y_value * (1 - p) ** (N - y_value)

    table = np.zeros(N + 1)
    for y_value in range(N + 1):
        numerator = quad(
            kernel, 0, np.inf, args=(y_value,), epsabs=1e-12,
            epsrel=1e-10, limit=2_000,
        )[0]
        denominator = quad(
            kernel, -np.inf, np.inf, args=(y_value,), epsabs=1e-12,
            epsrel=1e-10, limit=2_000,
        )[0]
        table[y_value] = numerator / denominator
    return table


PP_NB = probpos_noborrow_table()


def run_hier(y, a_sigma, seed, n_iter=N_ITER, burn=BURN):
    """Vectorized Metropolis-within-Gibbs across simulated trials."""
    rng = np.random.default_rng(seed)
    n_trials = y.shape[0]
    gamma = np.zeros((n_trials, NB))
    mu = np.zeros(n_trials)
    sigma = np.full(n_trials, 0.5)
    ll = binom_ll(y, gamma)
    step_gamma, step_log_sigma = 0.5, 0.35
    posterior_positive = np.zeros((n_trials, NB))
    gamma_accept = 0
    sigma_accept = 0
    kept = 0

    for iteration in range(n_iter):
        proposal = gamma + rng.normal(0, step_gamma, (n_trials, NB))
        ll_proposal = binom_ll(y, proposal)
        lp_old = -0.5 * ((gamma - mu[:, None]) / sigma[:, None]) ** 2
        lp_new = -0.5 * ((proposal - mu[:, None]) / sigma[:, None]) ** 2
        accept = np.log(rng.random((n_trials, NB))) < (
            (ll_proposal + lp_new) - (ll + lp_old)
        )
        gamma_accept += int(accept.sum())
        gamma = np.where(accept, proposal, gamma)
        ll = np.where(accept, ll_proposal, ll)

        precision = 1 / 1e4 + NB / sigma**2
        mu = rng.normal(
            (gamma.sum(axis=1) / sigma**2) / precision,
            np.sqrt(1 / precision),
        )

        log_sigma = np.log(sigma)
        proposal_log_sigma = log_sigma + rng.normal(0, step_log_sigma, n_trials)
        proposal_sigma = np.exp(proposal_log_sigma)

        def log_target(value, log_value):
            return (
                -NB * np.log(value)
                - 0.5 * np.sum((gamma - mu[:, None]) ** 2, axis=1) / value**2
                - value**2 / (2 * a_sigma**2)
                + log_value
            )

        accept_sigma = np.log(rng.random(n_trials)) < (
            log_target(proposal_sigma, proposal_log_sigma)
            - log_target(sigma, log_sigma)
        )
        sigma_accept += int(accept_sigma.sum())
        sigma = np.where(accept_sigma, proposal_sigma, sigma)

        if iteration >= burn:
            posterior_positive += gamma > 0
            kept += 1

    return {
        "probability": posterior_positive / kept,
        "gamma_acceptance": gamma_accept / (n_iter * n_trials * NB),
        "sigma_acceptance": sigma_accept / (n_iter * n_trials),
    }


def load_shared_inputs():
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            "oc_simulation_inputs.csv is missing; run oc_simulation.R first"
        )
    grouped = {
        phase: {scenario: [] for scenario in SCEN}
        for phase in ("calibration", "evaluation")
    }
    with INPUT_PATH.open(newline="") as handle:
        for row in csv.DictReader(handle):
            phase = row["phase"]
            scenario = int(row["scenario"])
            grouped[phase][scenario].append(
                [int(row[f"b{basket}"]) for basket in range(1, NB + 1)]
            )
    output = {
        phase: {scenario: np.asarray(rows, dtype=float) for scenario, rows in scenarios.items()}
        for phase, scenarios in grouped.items()
    }
    for scenario in SCEN:
        if output["calibration"][scenario].shape != (R_CAL, NB):
            raise ValueError(f"unexpected calibration input shape for scenario {scenario}")
        if output["evaluation"][scenario].shape != (R_EVAL, NB):
            raise ValueError(f"unexpected evaluation input shape for scenario {scenario}")
    return output


inputs = load_shared_inputs()
prob_cal = {}
prob_eval = {}
diagnostics = []

for scenario in SCEN:
    y_all = np.vstack([
        inputs["calibration"][scenario],
        inputs["evaluation"][scenario],
    ])
    no_borrow = PP_NB[y_all.astype(int)]
    prob_cal[(scenario, "No")] = no_borrow[:R_CAL]
    prob_eval[(scenario, "No")] = no_borrow[R_CAL:]

    for borrow, seed_base in (("Moderate", 51_000), ("Strong", 61_000)):
        fit = run_hier(y_all, A_SIGMA[borrow], seed=seed_base + scenario)
        prob_cal[(scenario, borrow)] = fit["probability"][:R_CAL]
        prob_eval[(scenario, borrow)] = fit["probability"][R_CAL:]
        diagnostics.append({
            "scenario": scenario,
            "borrow": borrow,
            "trials": int(y_all.shape[0]),
            "iterations": N_ITER,
            "burn_in": BURN,
            "retained_draws_per_trial": N_ITER - BURN,
            "gamma_acceptance": fit["gamma_acceptance"],
            "sigma_acceptance": fit["sigma_acceptance"],
        })
    print(f"posterior complete for scenario {scenario}", flush=True)


def empirical_cutoff(values, target=0.05):
    return float(np.quantile(values, 1 - target, method="inverted_cdf"))


thresholds = {"weak": {}, "strong": {}}
for borrow in ("No", "Moderate", "Strong"):
    thresholds["weak"][borrow] = empirical_cutoff(
        prob_cal[(1, borrow)].max(axis=1)
    )
    candidates = []
    for scenario in SCEN:
        nulls = NULLSET[scenario]
        if nulls:
            candidates.append(empirical_cutoff(
                prob_cal[(scenario, borrow)][:, nulls].max(axis=1)
            ))
    thresholds["strong"][borrow] = max(candidates)

print("\nCalibrated q (Python, R-shared calibration inputs, paper reference):")
for control in ("weak", "strong"):
    for borrow in ("No", "Moderate", "Strong"):
        print(
            f"  {control:<6} {borrow:<9} python={thresholds[control][borrow]:.4f} "
            f"paper={Q_PAPER[control][borrow]:.3f}"
        )


def wilson(indicator, z=1.959963984540054):
    indicator = np.asarray(indicator, dtype=bool)
    n_trials = indicator.size
    estimate = indicator.mean()
    denominator = 1 + z**2 / n_trials
    center = (estimate + z**2 / (2 * n_trials)) / denominator
    half_width = z * np.sqrt(
        estimate * (1 - estimate) / n_trials + z**2 / (4 * n_trials**2)
    ) / denominator
    return {
        "estimate": float(estimate * 100),
        "low": float(max(0, center - half_width) * 100),
        "high": float(min(1, center + half_width) * 100),
    }


def metrics(reject, scenario):
    active = PROMISING[scenario]
    nulls = NULLSET[scenario]
    return {
        "basket": [wilson(reject[:, basket]) for basket in range(NB)],
        "FWER": wilson(reject[:, nulls].any(axis=1)) if nulls else None,
        "FWPD": wilson(reject[:, active].any(axis=1)) if active else None,
        "FWPC": wilson(reject[:, active].all(axis=1)) if active else None,
    }


results = []
for control in ("weak", "strong"):
    for scenario in SCEN:
        for borrow in ("No", "Moderate", "Strong"):
            reject = prob_eval[(scenario, borrow)] > thresholds[control][borrow]
            result = metrics(reject, scenario)
            results.append({
                "control": control,
                "scenario": scenario,
                "borrow": borrow,
                "threshold": thresholds[control][borrow],
                "calibration_trials": R_CAL,
                "evaluation_trials": R_EVAL,
                **result,
            })

output = {
    "metadata": {
        "calibration_trials_per_scenario": R_CAL,
        "evaluation_trials_per_scenario": R_EVAL,
        "mcmc_iterations": N_ITER,
        "mcmc_burn_in": BURN,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "input_file": INPUT_PATH.name,
    },
    "thresholds": thresholds,
    "diagnostics": diagnostics,
    "results": results,
}
(SCRIPT_DIR / "oc_results.json").write_text(json.dumps(output, indent=2) + "\n")

for control in ("weak", "strong"):
    print(f"\n===== {control.upper()} FWER TARGET: INDEPENDENT EVALUATION =====")
    print(f"{'Scen':<5}{'Borrow':<10}{'%reject (b1,b2,b3,b4)':<28}{'FWER':>8}{'FWP-D':>8}{'FWP-C':>8}")
    for result in results:
        if result["control"] != control:
            continue
        basket = ",".join(f"{x['estimate']:4.1f}" for x in result["basket"])
        fmt = lambda value: "   -   " if value is None else f"{value['estimate']:6.1f}"
        print(
            f"S{result['scenario']:<4}{result['borrow']:<10}{basket:<28}"
            f"{fmt(result['FWER']):>8}{fmt(result['FWPD']):>8}{fmt(result['FWPC']):>8}"
        )
