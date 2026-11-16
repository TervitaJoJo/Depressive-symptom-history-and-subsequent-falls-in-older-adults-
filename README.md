# CHARLS and ELSA analysis code: README

Code accompanying **Depressive symptom history and subsequent falls in older adults with low current symptoms and no recent falls: a longitudinal study in CHARLS and ELSA**.

This package provides one workflow from licensed source files to analytic samples, cohort-specific models, sensitivity analyses, table sources, figures and numerical checks. The two cohorts are analyzed separately. Included CSV files contain aggregate reference results and diagnostics; participant data must be obtained independently from the data providers.

**Start here.** Run the commands below from the `code/` directory. This README is stored in `workflow/` so that the code root contains only project folders. Detailed module guides are linked below.

## 1. Study and analysis scope

The study asks whether elevated depressive symptoms at the assessment preceding the index visit are associated with subsequent reported falls among adults aged 65 years or older who have low current symptoms and no recent falls.

| Feature | CHARLS | ELSA |
|---|---|---|
| Historical elevated symptoms | CES-D10 score of 10 or greater | CES-D8 score of 4 or greater |
| Current low symptoms | CES-D10 score below 10 | CES-D8 score below 4 |
| Assessment windows | 2011-2013-2015 and 2013-2015-2018 | Index waves 5-8, with the preceding and following waves |
| Primary observed sample | 2,475 participants; 406 subsequent fall reports | 4,779 participants; 1,134 subsequent fall reports |
| Main estimation | Modified Poisson models with household-clustered covariance and multiple imputation | Same modeling approach, fitted independently |

Each participant contributes one selected window. The workflow preserves the original ordering of eligibility and current-score restrictions; details are in the data-preparation guide. ELSA supplies independent evidence about the association. The primary estimates are not pooled, and the prediction comparisons are internal to each cohort.

Supporting analyses examine current symptom severity, function, historical falls, interview intervals, follow-up observation, medication and vulnerability measures, recorded death, symptom-score definitions, longer symptom history and incremental prediction. Dated analysis plans and source amendments are retained in the relevant modules. Their timing is unchanged by this integrated workflow; the package is not a preregistration.

## 2. Package contents

```text
code/
  data_preparation/                 Build analytic samples from licensed files
  charls_elsa_models/                Main, four-visit, targeted and prediction analyses
    aggregate_results/              Reference estimates and diagnostics
  historical_falls_sensitivity/      Historical falls and interval-restricted absolute risks
    aggregate_results/              Reference estimates and diagnostics
    provenance/                     Dated plans, source amendment and checksums
  elsa_threshold_decomposition/      ELSA eligibility/exposure-definition decomposition
    aggregate_results/              Reference estimates and diagnostics
  workflow/                         Run, export, render and verify the complete workflow
    README_overview.md              This guide
    run_workflow.py                 Single execution entry point
    OUTPUT_MAP.md                   Manuscript table/figure correspondence
    REPRODUCTION.md                 Scope of completed reproduction checks
```

| Module | Guide | Recursive file count |
|---|---|---:|
| `data_preparation/` | [Data preparation](../data_preparation/README.md) | 6 |
| `charls_elsa_models/` | [Models](../charls_elsa_models/README.md) | 40 |
| `historical_falls_sensitivity/` | [Historical-fall analyses](../historical_falls_sensitivity/README.md) | 19 |
| `elsa_threshold_decomposition/` | [Threshold decomposition](../elsa_threshold_decomposition/README.md) | 22 |
| `workflow/` | [Stage dependencies](README.md) | 10 |
| **Total** | Including all nested files | **97** |

Every immediate project folder contains at most 100 files recursively. Keep environments, generated outputs, logs and caches outside the submission package. The existing `aggregate_results/` folders are reference material for comparison and should remain unchanged.

## 3. Software environment

The recorded verification environment used **Python 3.14.5** and **R 4.6.0**. Python package versions are listed in [ENVIRONMENT.json](../charls_elsa_models/ENVIRONMENT.json). Recorded R package versions include `mice` 3.19.0, `sandwich` 3.1-3, `mitml` 0.4-5 and `jsonlite` 2.0.0. R session records are supplied with the aggregate diagnostics.

Use a Python environment outside this package, then install the listed dependencies:

```text
python -m pip install -r charls_elsa_models/requirements.txt
Rscript -e "install.packages(c('mice','sandwich','mitml','jsonlite'))"
```

The Python requirements include NumPy, pandas, Matplotlib, Pillow, scikit-learn, threadpoolctl and pyreadstat. `splines` and `parallel` are supplied with R. These installation commands do not pin versions; consult the recorded environment when investigating numerical differences. `Rscript` must be available on `PATH`, or supplied with the runner's `--rscript` option.

## 4. Prepare licensed source files

Create a private input directory using the following provider filenames and layout:

```text
licensed_data/
  CHARLS/
    Harmonized_CHARLS/
      H_CHARLS_D_Data.dta
    2011/Health_Status_and_Functioning.dta
    2013/Health_Status_and_Functioning.dta
    2015/Health_Status_and_Functioning.dta
    2018/Health_Status_and_Functioning.dta
  ELSA/
    Raw_data/UKDA-5050-stata/stata/stata13_se/
      h_elsa_g3.dta
      wave6/wave_6_elsa_nurse_data_v2.dta
```

Use the displayed filename capitalization on case-sensitive systems. The workflow starts from Harmonized CHARLS D, Harmonized ELSA G3 and selected raw survey modules. It neither downloads restricted data nor reconstructs the providers' harmonized products. Code access does not grant data access.

For construction and variable definitions, see:

- [Eligibility, linkage and source preparation](../data_preparation/README.md)
- [Primary data dictionary](../charls_elsa_models/DATA_DICTIONARY.md)
- [Medication, vulnerability and recorded-death variables](../charls_elsa_models/EXTENDED_DATA_DICTIONARY.md)
- [Historical-fall sources and recall periods](../historical_falls_sensitivity/README.md)

## 5. Run the complete workflow

The examples assume that `licensed_data/` and a new `reproduction_private/` directory are siblings of `submission_package/`:

```text
project/
  licensed_data/
  reproduction_private/             Created by the workflow
  submission_package/
    code/                            Run the commands here
```

First inspect the planned stages without reading data or creating outputs:

```text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --dry-run
```

Then run all stages:

```text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private"
```

Replace the two relative paths with your private input and output locations when using a different layout. Generated outputs must remain outside the submission package; the scripts enforce this boundary. A dry run checks stage selection and script paths, but does not verify data availability or installed dependencies.

The default sequence is:

| Stage | Purpose | Output within the private work directory |
|---|---|---|
| `prepare` | Derive eligible and observed samples | `prepared/` |
| `enrich` | Add medication, vulnerability and recorded-death measures | `extended/` |
| `main` | Fit primary and standard sensitivity models | `results/main/`, `cache/main/` |
| `four` | Fit the four-visit CHARLS analysis | `results/four/` |
| `targeted` | Fit directed sensitivities, continuous-score and spline models | `results/targeted/` |
| `history_inputs` | Link historical-fall responses to selected windows | `history_inputs/` |
| `history_falls` | Fit paired historical-fall models and interval-restricted risks | `results/historical_falls/` |
| `prediction` | Compare internally cross-validated prediction models | `results/prediction/` |
| `decomposition` | Separate descriptive eligibility and history-definition changes | `results/decomposition/` |
| `tables` | Export permitted aggregate summaries and table sources | Three aggregate folders and `tables/` |
| `figures` | Render Figures 1-2 and Figure S1 | `figures/` |
| `compare` | Compare exported results with 53 reference CSVs | `numerical_comparison.json` |

A failed stage stops the workflow. Its output is recorded in `logs/<stage>.log` under the private work directory. Successful completion of the default full run includes the numerical comparison stage.

## 6. Computation settings and partial runs

The analysis uses 30 imputations with 20 iterations. Default runner settings are 500 prediction-bootstrap replicates per cohort, 1,000 paired household-bootstrap replicates per imputation for threshold decomposition, and 10 decomposition workers. A fresh decomposition requires 120,000 bootstrap model fits across four cells. The complete run can therefore require substantial computation.

| Option | Default | Purpose |
|---|---:|---|
| `--prediction-bootstrap` | 500 | Prediction uncertainty calculations |
| `--decomposition-bootstrap` | 1000 | Paired bootstrap replicates per imputation |
| `--workers` | 10 | Decomposition parallel workers |
| `--start-at` | `prepare` | First stage to execute |
| `--stop-after` | `compare` | Last stage to execute |
| `--rscript` | `Rscript` | R executable |

Reduced bootstrap counts can help check execution but do not reproduce the reported uncertainty. Seeds are defined in the analysis scripts.

Once all model outputs exist, regenerate tables and figures and run the numerical comparison:

```text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --start-at tables
```

With prepared samples and compatible main-model caches available, run only historical-fall preparation and estimation:

```text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --start-at history_inputs --stop-after history_falls
```

Check an already exported run directly:

```text
python -B workflow/compare_results.py "../../reproduction_private"
```

Stage selection does not construct omitted prerequisites. Reuse caches only when inputs, model definitions, seeds and settings match. Use a new private directory for a complete rerun; the four-visit model refuses to overwrite a completed result directory.

## 7. Find the manuscript outputs

The [complete output map](OUTPUT_MAP.md) identifies the scripts and exported source files for each display.

| Display | Information |
|---|---|
| Figure 1; Table 1 | Assessment design, sample flow and participant characteristics |
| Table 2; Figure 2 | Sequentially adjusted risk ratios and standardized absolute risks |
| Table S1 | ELSA threshold decomposition |
| Table S2 | Missingness and follow-up observation |
| Table S3 | Conditional-information tests and descriptive cohort comparison |
| Table S4 | Restrictions, weights, secondary outcomes, interval-restricted risks and historical falls |
| Table S5 | Medication, measured vulnerability and E-values |
| Table S6 | Recorded deaths, unavailable outcomes and composite-outcome sensitivity |
| Table S7; Figure S1 | Continuous historical scores, association shape and prediction increment |
| Table S8 | Four-visit CHARLS symptom history |

The `tables/` outputs are numerical sources for manuscript assembly. The workflow does not rebuild the Word manuscript or its administrative documents.

Exported main and targeted summaries are in `aggregate_results/`; historical-fall summaries are in `historical_falls_aggregate/`; threshold summaries are in `decomposition_aggregate/`. These paths are relative to the private work directory. `results/`, `prepared/`, `extended/`, `history_inputs/`, caches and logs can contain participant information and must remain private.

## 8. Verify and interpret a run

Expected observed sample counts are shown in Section 1. The four-visit CHARLS sample contains 1,443 participants and 230 subsequent fall reports. Larger eligibility extracts are expected before current-score and follow-up restrictions. Failed count assertions should be traced to source releases, definitions or linkage rather than disabled.

`compare_results.py` checks 53 reference CSVs using `rtol=1e-8` and `atol=1e-10`, trims text-cell whitespace and excludes wall-clock duration columns. A mismatch produces a nonzero exit status and a JSON report. The check covers the included aggregate estimates and diagnostics, rather than every manuscript sentence or methodological assumption.

[REPRODUCTION.md](REPRODUCTION.md) records what has actually been re-estimated, which computations reused matching caches, and which checks verified export and integration. The recorded combined export check passed all 53 comparisons. It did not freshly regenerate every imputation and decomposition bootstrap. This documentation update also does not constitute a new model run.

Interpret estimates as conditional associations. Cohort recall periods differ; the recorded-death analysis is a composite-outcome sensitivity; prediction comparisons use internal cross-validation. Consult the manuscript and module guides for these distinctions. Dated provenance files and inconclusive results remain part of the release.

## 9. Troubleshooting and anonymous sharing

| Issue | Action |
|---|---|
| Missing source file | Check provider release, directory layout and filename capitalization |
| `Rscript` unavailable | Configure `PATH` or pass the installed executable through `--rscript` |
| Output-path rejection | Choose a private directory outside the submission package |
| Failure after skipping stages | Restore the prerequisites listed in the [stage reference](README.md) |
| Existing completed four-visit output | Use a new work directory for a full rerun |
| Numerical comparison failure | Inspect the JSON report, input releases, environment and cache compatibility |
| Unexpected fit or imputation messages | Inspect private stage logs and supplied model diagnostics before interpreting results |

The release contains code, documentation and aggregate outputs for **CHARLS and ELSA only**. It includes no restricted source data, participant records, identifiers, imputations or individual prediction caches. Review newly generated aggregate files under the providers' disclosure rules before sharing.

Preserve anonymity by excluding author details, institutional affiliations, contact information, machine-specific absolute paths, credentials and unredacted logs. Keep the release below the recursive file limit, and store backups and generated files separately. Changes to this local package do not automatically update a remote repository.
