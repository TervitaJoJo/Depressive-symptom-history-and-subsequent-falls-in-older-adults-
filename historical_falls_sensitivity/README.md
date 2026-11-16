# Historical falls and interval-restricted absolute risks

This module evaluates historical falls and observation intervals within the same two-cohort workflow. It preserves the fixed primary samples, outcomes, thresholds, selected windows and primary fits. These sensitivity analyses were specified after the primary results and are exploratory. The dated plan, source correction and their original hashes are retained unchanged in `provenance/analysis_plan.md`, `provenance/metadata_amendment.md` and their corresponding lock files.

## Reproduction

Start with the [root README](../workflow/README_overview.md). The full workflow runs `history_inputs` and `history_falls` between targeted estimation and prediction. To run these stages alone, a private work directory must already contain the original `prepared/` files and `cache/main/` imputation objects. From the `code/` directory:

```text
python -B workflow/run_workflow.py --raw-root ../../licensed_data --work-dir ../../reproduction_private --rscript Rscript --start-at history_inputs --stop-after history_falls
```

The input and output paths are examples relative to the current working directory. All generated records, IDs, caches and logs stay outside the distributed code folder. Authorized CHARLS and ELSA data are required. No participant-level records are distributed.

`build_history_falls.py` links falls from the historical visit to each original selected window and checks every original column for identity. `fit_historical_falls.R` imports the frozen estimator helpers without executing the original driver, verifies the original primary and interval-restriction estimates, then fits the paired sensitivity comparisons. Differences in original N, event counts or estimates stop the script.

## Models and interpretation

- `frozen_primary_check`: numeric reproduction only; not a new analysis.
- `both_intervals_18_30`: original interval restriction and imputations, adding standardized risks and risk differences.
- `same_sample_original` and `same_sample_plus_history_fall`: identical observed-history samples, identical 30 fresh completed datasets, original full model with/without the binary historical-fall variable.
- `no_falls_at_history_and_index`: subset of the same completed data, included only when the recorded pre-fit feasibility criteria are satisfied.

The estimator is modified Poisson with household-clustered HC1 variance and finite-sample Rubin pooling. Risks average fitted means; risk differences use delta-method variances. The same-sample log-RR change is descriptive, not a causal contribution or mediation effect. Symptoms and falls at the history visit need not have a known temporal order. Some CHARLS historical-fall-adjusted Poisson means exceed one; diagnostics are supplied and these are not individual prediction models.

## Source definitions

| Cohort/wave | Historical source | Recall period |
|---|---|---|
| CHARLS 2011 | Raw DA023, 1 yes / 2 no; IDs aligned through ID_w1 | Previous two years |
| CHARLS 2013 | Raw DA023, 1 yes / 2 no | Since previous interview for returning participants |
| ELSA wave 4 | Harmonized r4fall1y, 1 yes / 0 no | Previous year |
| ELSA waves 5-7 | Harmonized r5fall-r7fall | Previous two years |

There is no `r4fall` variable. The wave-4 mapping amendment was documented before the new association fits. The original ELSA index waves 5-8 and follow-up waves 6-9 retain their two-year outcome/eligibility indicators. The model retains the original analytic-window covariate.

CHARLS main index/follow-up questions for returning participants ask about falls since the preceding interview, including the approximately three-year 2015-2018 interval. ELSA main questions use a fixed two-year recall period. The 18-30-month interview restriction does not make these recall periods identical.

Sources: original CHARLS 2011 DA023 questionnaire; CHARLS 2013 Wave 2 PDF page 82; CHARLS 2015 PDF page 68; CHARLS 2018 household PDF page 81; Harmonized ELSA G3 documentation PDF page 356. The source variable/column inventory and hashes are generated in the private linkage output; source data are never redistributed.

## Frozen aggregate outputs

`aggregate_results/` contains numeric model results, sample/availability counts, comparison checks and model/imputation diagnostics. The cohort-specific sensitivity imputations use 30 datasets and 20 iterations with fixed seeds recorded in the script. Missing historical falls are excluded only in this sensitivity. `capacity_and_missingness.csv` records the restriction eligibility check. `paired_point_estimate_changes.csv` is a descriptive comparison, with no claimed inferential or causal decomposition.

The publication locations are Table S4 panels B-C and the wave-definition and sensitivity-methods sections of Online Resource 1. All fixed main and interval-check results must pass before sensitivity estimates are exported.
