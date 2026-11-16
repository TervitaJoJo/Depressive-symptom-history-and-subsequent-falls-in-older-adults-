# ELSA threshold decomposition (post hoc)

Run with R 4.6.0 and packages mice, sandwich, splines, jsonlite, parallel:

```text
Rscript threshold_decomposition.R PRIVATE_ELSA_INPUT.csv PRIVATE_CACHE PRIVATE_RESULTS 1000 10
```

These are argument templates; writable paths must be outside the submission package. The complete runner and preparation instructions are in the [root README](../workflow/README_overview.md); stage dependencies are in [workflow/README.md](../workflow/README.md). The five arguments are the authorized analytic CSV, private cache directory, aggregate-output directory, number of bootstrap draws per imputation, and worker count. Defaults for the last two are 1000 and 10. See the saved plan for timing, contrasts and boundaries. The paper uses the MI_Boot rows of pooled_contrasts.csv. Joint_HC1 rows are an independent covariance cross-check. Within/between/total covariance matrices are in A/B/C/D column order. exp_estimate for a contrast is a ratio of fitted RRs, not an individual fall risk ratio.

Input: the same selected-window ELSA CSV consumed by the existing primary model (see ../charls_elsa_models/DATA_DICTIONARY.md). Required fields: id, cluster, fall, cesd_current, cesd_prior, age, female, education, window, exposure_gap_months, followup_gap_months, partnered, smoking, alcohol, stroke, hearte, diabe, arthre, stair_difficulty, adl_difficulty. One selected window per person. Keep identifiers as character values. No new window selection or outcome construction is performed here. The exact observed primary broad population has 4779 participants; assertions intentionally check this paper-specific reproduction.

Missing covariates only are imputed; outcomes and scores are observed. Thirty imputations, twenty iterations, common seed 20261101. Bootstrap stream m uses seed 20261101+1000*m. Household frequencies weight model fits; they are identical across four cells within each draw. Both eligibility paths are reported. The main estimates and prior sensitivity results are unchanged.

Private caches must remain outside any shared package. Reuse caches only for identical data and unchanged methods; use a fresh cache folder after changing the script or settings. Runtime/version information is in aggregate_results/sessionInfo.txt. The supplied outputs contain aggregate estimates and diagnostics, not participant records. Original_pipeline_comparison.csv additionally quotes the frozen primary and original threshold-3 CSVs to distinguish them from common-pipeline estimates.

Method sources: Schomaker & Heumann (2018), DOI 10.1002/sim.7654; Shorrocks (2013), DOI 10.1007/s10888-011-9214-z. See supplementary methods for complete references. This is a descriptive decomposition of conditional associations, not a causal mediation analysis or attribution of clinical mechanisms.
