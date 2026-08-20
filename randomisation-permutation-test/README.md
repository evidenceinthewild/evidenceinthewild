# The Permutation Test Nobody Ran

Reconstruction-based permutation analyses for Study 201 (NCT01396239), the pivotal eteplirsen trial in Duchenne muscular dystrophy. Companion analysis to [Randomisation in the Wild Part 4](https://evidenceinthewild.com/randomisation-permutation-test/).

> **Scope correction (14 August 2026):** This repository does not contain the original individual patient data and does not recover a uniquely determined permutation p-value from the published record. The reported p = 0.029 is conditional on one constructed Week 24 outcome configuration. See [Provenance and inferential scope](PROVENANCE.md) and the [publication amendment](AMENDMENT-2026-08-14.md).

## What this is

Mendell et al. (2013) randomised 12 boys to three arms (30 mg/kg, 50 mg/kg, placebo) and analysed the six-minute walk test (6MWT) using a mixed model for repeated measures (MMRM). The MMRM found no significant difference for 50 mg/kg vs placebo at Week 24 (approximate two-sided p ≈ 0.56, derived from published adjusted means and standard errors).

For one constructed patient-level configuration compatible with selected published summaries, exhaustive enumeration of all 70 four-versus-four allocations gives a two-sided p-value of 0.029 for a raw Week 24 difference in means.

That result is not the same analysis as the published MMRM. The procedures differ in arm set, time structure, estimand, test statistic, and assumptions. The MMRM uses repeated measurements and all three randomised arms; the permutation calculation here uses a pairwise Week 24 contrast. It should be read as a reconstruction-based sensitivity analysis, not as a recovered analysis of the original patient-level data or a replacement for the published model.

## Important: reconstructed data

**The patient-level data in this repository are not original clinical data.** Individual Week 24 6MWT values were not published in tabular form. The baseline values in `data/reconstructed_patient_data.csv` were selected to match the arm-level means, standard deviations, minima, and maxima in Table 1 of Mendell et al. (2013), subject to integer-metre and published-rounding assumptions. The paper directly states the baselines for patients 009 and 010. The remaining baselines are reconstructed rather than observed.

The Week 24 values are also reconstructed. Figure 6 reports a combined six-patient ambulation-evaluable eteplirsen trajectory rather than a separate 50 mg/kg raw mean, so the 50 mg/kg mean and its four individual changes must not be described as values reported by Figure 6. They are analysis inputs used to construct one configuration. The row-level `source` field in the CSV has been replaced by separate provenance fields so directly reported and reconstructed elements are not conflated.

The sensitivity analysis (`code/sensitivity_analysis.py`) varies individual outcomes while preserving selected group-mean targets. Its ranges are deliberately broad but ad hoc. The resulting percentages describe the output of that generator; they are not probabilities that the unknown original-data p-value lies below a threshold. For the 50 mg/kg comparison, the median p-value is 0.086 and 25% of generated configurations fall below 0.05, so p = 0.029 is not stable across the generated configurations. With four participants per group, achievable two-sided p-values occur in steps of 2/70: 0.0286, 0.0571, 0.0857, and so on. Thus, 25% below 0.05 means exactly that 25% reached the floor of 2/70; no value between 0.0286 and 0.0571 is possible.

The mITT calculation is exploratory. It excludes two participants according to post-randomisation loss of ambulation. A simple reassignment test within that selected set is not automatically justified by the original randomisation and should not be called design-exact without a selection-aware analysis.

## What “exact” means here

The enumeration is computationally exact: all 70 allocations of the eight selected outcome values into groups of four are evaluated. The inferential interpretation is conditional on:

- the reconstructed outcome vector;
- the raw difference in means as the test statistic;
- a sharp null of no individual treatment effect;
- the selected 50 mg/kg and placebo participants; and
- the assumption that the relevant four-of-eight allocations represent the trial's assignment mechanism.

Exact enumeration does not make reconstructed outcomes exact, and it does not make two analyses with different questions interchangeable.

## Reproducing the results

```bash
pip install -r code/requirements.txt
python code/eteplirsen_permutation_analysis.py   # main analysis → results/
python code/sensitivity_analysis.py              # generator-based sensitivity → results/
python code/create_visualizations.py             # figures → figures/png/, figures/pdf/
python code/create_blog_figures.py                # title-free figures → figures/blog/
```

The analysis and sensitivity scripts must run before the figure renderers. Paths are resolved relative to each script, so the commands may be invoked from any working directory. Python 3.8+ required.

## Repository structure

```
data/
  reconstructed_patient_data.csv       Reconstructed 12-patient dataset with provenance
code/
  eteplirsen_permutation_analysis.py   Data reconstruction, verification, permutation tests
  sensitivity_analysis.py              Generator-based sensitivity across 1,000 configurations
  create_visualizations.py             Publication figures
  create_blog_figures.py               Title-free figures for the blog post
  requirements.txt                     Python dependencies
results/                               (generated by scripts)
  permutation_results.json             Machine-readable results
  perm_dist_*.npy                      Permutation distributions
  sensitivity_*.npy                    Sensitivity analysis p-value arrays
figures/                               (generated by scripts)
  png/                                 Publication-quality PNG figures
  pdf/                                 Vector PDF figures
  blog/                                Title-free PNG figures for Ghost captions
```

## Key results

| Comparison | Displayed difference | Permutation p (two-sided) | MMRM p |
|---|---|---|---|
| 50 mg/kg vs placebo | +67.5m | 0.029, conditional on the displayed reconstruction | ≈ 0.56 (derived; different analysis) |
| mITT eteplirsen (n=6) vs placebo | +70.0m | 0.005, exploratory after post-randomisation selection | -- |
| All eteplirsen (ITT) vs placebo | +2.9m | 0.97 | -- |

Under the ad hoc sensitivity generator, 90% of mITT configurations yield two-sided p < 0.05. That does not resolve the post-randomisation-selection problem. For the dose-arm comparison, 25% fall below p < 0.05 and 67% below p < 0.10; the median is 0.086. On the 4-v-4 grid, every p-value below 0.05 is the minimum 2/70.

## References

Mendell JR, Rodino-Klapac LR, Sahenk Z, et al. Eteplirsen for the treatment of Duchenne muscular dystrophy. *Ann Neurol*. 2013;74(5):637-647. doi:10.1002/ana.23982

Sarepta Therapeutics. Presentations for the April 25, 2016 Meeting of the Peripheral and Central Nervous System Drugs Advisory Committee. In particular, slide CE-23 presents the combined ambulation-evaluable and lost-ambulation trajectories; it does not tabulate individual Week 24 outcomes or a separate 50 mg/kg raw mean.

US Food and Drug Administration. Clinical Review, NDA 206488. This source is listed in the analysis code and must be cited at page level wherever a value is taken from it.

Fisher RA. *The Design of Experiments*. Edinburgh: Oliver and Boyd; 1935.

Horn PS. Some easy t statistics. *J Am Stat Assoc*. 1983;78(384):930-936. doi:10.1080/01621459.1983.10477042.

## License

Analysis code is provided for transparency and reproducibility. The constructed patient records are analysis inputs derived from or constrained by published information; they are not original clinical data and must not be represented as recovered observations.
