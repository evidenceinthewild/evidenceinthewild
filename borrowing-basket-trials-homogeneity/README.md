# Borrowing in Basket Trials: A Bet on Homogeneity

Reproducible companion analyses for the Evidence in the Wild post
“Borrowing in Basket Trials: A Bet on Homogeneity.”

These are independent, seeded reimplementations of examples discussed in
Zhou and Ji (2024), not the original authors' code.

See [`VALIDATION.md`](VALIDATION.md) for the independent audit, corrected
calibration design, execution environment, and cross-language results.

## Repository structure

```text
analysis/
├── borrowing-dial/
│   ├── figure_borrowing_dial.py
│   └── figure_borrowing_dial.png
├── imatinib-shrinkage/
│   ├── imatinib_shrinkage.R
│   ├── imatinib_validate.py
│   ├── imatinib_results.csv
│   ├── imatinib_diagnostics.csv
│   ├── imatinib_validation.csv
│   ├── imatinib_validation_diagnostics.json
│   └── figure_imatinib_shrinkage.png
└── operating-characteristics/
    ├── oc_simulation.R
    ├── oc_validate.py
    ├── oc_results.csv
    ├── oc_results.json
    ├── oc_thresholds.csv
    ├── oc_diagnostics.csv
    ├── oc_simulation_inputs.csv
    ├── figure_oc_borrowing.png
    └── figure_oc_borrowing_reimplementation.png
```

Each analysis is self-contained. Generated outputs live beside the scripts
that create them, which allows every script to be run from any working
directory without path configuration.

## Analyses

### Borrowing dial

`figure_borrowing_dial.py` generates the conceptual continuum from full
pooling to separate analyses.

### Imatinib shrinkage

`imatinib_shrinkage.R` recreates the basket-level shrinkage illustration and
writes the figure, numerical results, and four-chain diagnostics.
`imatinib_validate.py` provides an independent Python implementation. Both
implementations run four chains and report Gelman-Rubin R-hat and sampler
acceptance rates.

### Operating characteristics

`oc_simulation.R` uses 2,500 calibration trials and a separate 5,000-trial
evaluation set for each of six scenarios. It saves the calibrated thresholds,
shared simulated basket counts, sampler diagnostics, estimates, Wilson 95%
Monte Carlo intervals, and both operating-characteristic figures.
`oc_validate.py` independently repeats the posterior calculations in Python
using the exact same simulated counts and writes `oc_results.json`.

The publication-facing `figure_oc_borrowing.png` uses the exact operating
characteristics reported by Zhou and Ji. The separate
`figure_oc_borrowing_reimplementation.png` records the seeded Monte Carlo
reimplementation. These two evidentiary layers are deliberately distinct.

## Reproduce the outputs

The R analyses use base R only. Run the R operating-characteristics script
before its Python validator because R creates the shared simulation-input file:

```sh
Rscript analysis/imatinib-shrinkage/imatinib_shrinkage.R
Rscript analysis/operating-characteristics/oc_simulation.R
```

For the Python analyses:

```sh
python3 -m pip install -r requirements.txt
python3 analysis/imatinib-shrinkage/imatinib_validate.py
python3 analysis/operating-characteristics/oc_validate.py
python3 analysis/borrowing-dial/figure_borrowing_dial.py
```

All stochastic analyses use fixed seeds. The R and Python OC implementations
use the same calibration and evaluation counts but separate MCMC streams.
Their posterior summaries should agree within Monte Carlo variation. The full
OC scripts are intentionally computation-heavy and may take several minutes.

## Statistical safeguards

- Threshold calibration and operating-characteristic evaluation use disjoint
  simulated trials.
- Every evaluation probability includes a Wilson 95% Monte Carlo interval.
- The no-borrowing posterior is integrated adaptively over the full real line.
- Calibrated thresholds, sampler acceptance rates, simulation inputs, and
  software versions are retained with the outputs.
- The publication-facing figure remains separate from the reimplementation.

## Software versions used for the checked run

- R 4.6.1
- Python 3.14.6
- NumPy 2.5.2
- SciPy 1.18.0
- Matplotlib 3.11.1

No open-source license has been selected for this repository.

## Reference

Zhou T, Ji Y. Bayesian methods for information borrowing in basket trials: an
overview. *Cancers*. 2024;16:251. <https://doi.org/10.3390/cancers16020251>
