# Validation report

Validated August 10, 2026.

## Scope

All five scripts were reviewed for data transcription, basket indexing,
likelihood and prior implementation, MCMC updates, threshold calibration,
operating-characteristic definitions, reproducibility, and figure rendering.
Every script was then run from a clean staged copy.

## Corrective changes

The original operating-characteristics scripts calibrated decision thresholds
and evaluated FWER and power on the same 1,000 simulated trials. That design
could not independently demonstrate error control. The corrected pipeline:

- estimates thresholds from 2,500 calibration trials per scenario;
- evaluates performance on 5,000 separate trials per scenario;
- reports Wilson 95% Monte Carlo intervals;
- saves the thresholds, simulation inputs, and sampler acceptance rates;
- uses adaptive integration over the full real line for the no-borrowing prior;
- supplies the same simulated basket counts to R and Python; and
- preserves the paper's exact reported values separately from the seeded
  reimplementation.

The imatinib analyses now use four chains, report R-hat and acceptance rates,
and write readable CSV/JSON output instead of a pickle-backed NumPy object.

## Results

### Imatinib shrinkage

| Implementation | Mean interval-width reduction | Maximum reported R-hat | Gamma acceptance | Sigma acceptance |
|---|---:|---:|---:|---:|
| R | 50.5% | 1.0046 | 44.0% | 55.8% |
| Python | 51.1% | 1.0046 | 43.6% | 55.7% |

Both implementations support the post's “about half” characterization.

### Headline operating characteristics

Values are percentages in no-, moderate-, and strong-borrowing order.

| Quantity | Published | R evaluation | Python evaluation |
|---|---:|---:|---:|
| Homogeneous scenario: disjunctive power | 69.6 / 87.0 / 90.3 | 65.5 / 86.6 / 91.5 | 65.5 / 86.1 / 90.3 |
| One inactive basket: weak-control FWER | 1.1 / 8.7 / 33.9 | 0.9 / 7.2 / 31.6 | 0.9 / 7.2 / 29.5 |
| One active basket: strong-control power | 27.0 / 16.9 / 4.0 | 23.9 / 18.2 / 3.8 | 23.9 / 18.8 / 4.4 |

The median absolute R/Python difference across all reported evaluation metrics
was 0.14 percentage points. Their largest difference was 3.54 points for a
low-probability conjunctive-power metric; the headline results differed by no
more than 2.1 points. Full Monte Carlo intervals are in `oc_results.csv` and
`oc_results.json`.

## Checked environment

- R 4.6.1
- Python 3.14.6
- NumPy 2.5.2
- SciPy 1.18.0
- Matplotlib 3.11.1

## Remaining limitations

The hierarchical OC simulation runs one MCMC stream per simulated trial, so
per-trial R-hat and effective-sample-size diagnostics would be prohibitively
large. The repository instead retains aggregate acceptance rates and uses a
second implementation on identical data as a sensitivity check. The estimates
remain Monte Carlo results, and the biological exchangeability assumption is
not validated by computation.
