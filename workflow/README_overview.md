# CHARLS and ELSA analysis workflow

Code accompanying **Depressive symptom history and subsequent falls in older adults with low current symptoms and no recent falls: a longitudinal study in CHARLS and ELSA**.

This is one workflow for source preparation, cohort-specific estimation, sensitivity analyses, table and figure generation, and numerical verification. Start here for installation and execution; module guides document implementation details. This overview is stored inside workflow/ so that code/ contains only project directories. Run all commands below from the code/ directory.

## Research scope

The study examines whether elevated depressive symptoms at the assessment preceding the index visit are associated with subsequent reported falls among adults aged 65 years or older who have low current symptoms and no recent falls. CHARLS and ELSA models are fitted separately. ELSA provides an independent assessment of the association; the main results are not pooled and this is not external validation of a prediction model.

Primary definitions, selected assessment windows and estimates are fixed. Supporting analyses address current symptom severity, functional status, historical falls, interview intervals, symptom-score thresholds, observation, recorded death and prediction increment. Analyses specified after the primary results remain exploratory. The dated plans and source correction are retained under the relevant modules; this repository is not a preregistration.

## Repository structure

~~~text
code/
  data_preparation/                    Licensed sources to analytic samples
  charls_elsa_models/                   Main, cumulative-history and targeted models
    aggregate_results/                 Frozen aggregate reference outputs
  historical_falls_sensitivity/         Historical falls and interval-restricted risks
    aggregate_results/                 Frozen aggregate reference outputs
    provenance/                        Unmodified dated plan, correction and hashes
  elsa_threshold_decomposition/         Eligibility and exposure-definition contrasts
    aggregate_results/                 Frozen aggregate reference outputs
  workflow/                            Runner, export, rendering and verification
    README_overview.md                 Complete workflow guide (this file)
~~~

| Module | Details |
|---|---|
| Data preparation | [Input layout and sample construction](../data_preparation/README.md) |
| Main and targeted estimation | [Model guide](../charls_elsa_models/README.md), [base variables](../charls_elsa_models/DATA_DICTIONARY.md), [extended variables](../charls_elsa_models/EXTENDED_DATA_DICTIONARY.md) |
| Historical falls | [Source definitions, paired comparisons and diagnostics](../historical_falls_sensitivity/README.md) |
| ELSA thresholds | [Four-cell decomposition](../elsa_threshold_decomposition/README.md) |
| Execution and manuscript outputs | [Stage dependencies](README.md), [table/figure map](OUTPUT_MAP.md), [verification record](REPRODUCTION.md) |

The release contains **97 files in total**, including nested and hidden files. Recursive counts are 40 for charls_elsa_models/, 6 for data_preparation/, 19 for historical_falls_sensitivity/, 22 for elsa_threshold_decomposition/ and 10 for workflow/, including this overview. Thus the root and every nested directory contain at most 100 files recursively. Store environments, generated outputs and caches outside the package; do not add archives or duplicate releases inside it.

## Software and licensed data

The recorded verification environment used Python 3.14.5 and R 4.6.0. Python packages and versions are recorded in [ENVIRONMENT.json](../charls_elsa_models/ENVIRONMENT.json). The recorded R packages were mice 3.19.0, sandwich 3.1-3, mitml 0.4-5 and jsonlite 2.0.0; splines and parallel ship with R. Session records are retained with the aggregate diagnostics.

Install dependencies in an environment **outside this package**. From the code directory:

~~~text
python -m pip install -r charls_elsa_models/requirements.txt
Rscript -e "install.packages(c('mice','sandwich','mitml','jsonlite'))"
~~~

These install commands are not version-pinned. Use the recorded environment when investigating numerical differences.

A complete run requires authorized access to Harmonized CHARLS D, CHARLS health-status/functioning modules for 2011, 2013, 2015 and 2018, Harmonized ELSA G3, and the specified wave-6 ELSA nurse module. The exact provider filenames and layout are in the [input guide](../data_preparation/README.md). The workflow starts from provider-released harmonized files and selected raw modules; it does not rebuild the harmonized products or download data. Data access is separate from code access.

No licensed source data, participant records, identifiers, imputations or bootstrap/prediction caches are distributed. The included CSVs are aggregate reference outputs and model diagnostics.

## Run the complete workflow

Use a layout such as:

~~~text
reproduction_workspace/
  licensed_data/                 Authorized local inputs
  reproduction_private/         Generated files; never upload wholesale
  submission_package/
    code/                       This release
~~~

Run from **submission_package/code/**:

~~~text
python -B workflow/run_workflow.py --help
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --dry-run
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private"
~~~

The full command runs all 12 stages, including historical-fall analyses and the final comparison of **53 aggregate CSV files**. The dry run checks selected script paths and lists their order; it does not read data, create outputs or validate data availability. Paths above are examples resolved from the current directory. Supply local paths at runtime rather than writing them into shared scripts. If Rscript is not on PATH, supply its executable through the --rscript option.

Choose a **new private work directory** for a full run. Output guards protect the parent directory of code/ and reject generated files within it. A failed stage stops the run; its log stays under the private work directory. The four-visit stage refuses to overwrite a completed result directory. The runner does not rebuild Word manuscript files or publish repository contents.

## Execution order and dependencies

| Order | Stage | Operation |
|---|---|---|
| 1 | prepare | Construct both cohort samples and the CHARLS four-visit sample |
| 2 | enrich | Derive medication, functional-domain and recorded-death variables |
| 3 | main | Fit primary, sequential, restricted, weighted and secondary-outcome models |
| 4 | four | Estimate CHARLS cumulative symptom-history associations |
| 5 | targeted | Estimate medication/domain sensitivities, E-values, continuous scores, splines and mortality-related sensitivities |
| 6 | history_inputs | Link historical falls to the same selected assessment windows |
| 7 | history_falls | Estimate paired historical-fall comparisons, two-assessment no-fall restrictions and interval-restricted absolute risks |
| 8 | prediction | Compare models using household-grouped cross-validation and bootstrap uncertainty |
| 9 | decomposition | Estimate ELSA threshold cells and descriptive component contrasts |
| 10 | tables | Collect allowlisted aggregate outputs and export Table 1–2 and S1–S8 numerical sources |
| 11 | figures | Render Figure 1–2 and Figure S1 |
| 12 | compare | Compare regenerated aggregate CSVs with frozen references; fail on discrepancies |

All sensitivity components participate in the default run. Historical-fall estimation uses the prepared samples, linked historical-fall inputs and matching main-model imputations. Table export requires all estimation stages. The [execution guide](README.md) records each stage's prerequisites and outputs.

Default settings are 30 imputations with 20 iterations, 500 prediction bootstrap replicates per cohort, and 1,000 paired household-bootstrap replicates per imputation for threshold decomposition. The decomposition entails 120,000 four-cell fits in a fresh run; default parallelism is 10 workers. Fresh runs can be computationally substantial. Options --workers, --prediction-bootstrap and --decomposition-bootstrap control resources, but reduced replicate counts do not reproduce the reported uncertainty.

Run a contiguous subset only after its prerequisites exist. With all models already estimated, export, render and verify together:

~~~text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --start-at tables
~~~

With prepared samples and main-model caches available, run the historical-fall stages:

~~~text
python -B workflow/run_workflow.py --raw-root "../../licensed_data" --work-dir "../../reproduction_private" --start-at history_inputs --stop-after history_falls
~~~

An existing exported run can also be checked directly:

~~~text
python -B workflow/compare_results.py "../../reproduction_private"
~~~

Reuse caches only when inputs, specifications, seeds and settings match. Caches are not distributed; use fresh private directories after substantive changes.

## Outputs and manuscript correspondence

All generated paths below are relative to the private work directory.

| Path | Contents |
|---|---|
| prepared/, extended/, history_inputs/, cache/ | Participant-level inputs, linkage information and computational caches; private |
| results/ | Estimation outputs mixed with participant audits, predictions or nested caches; do not share wholesale |
| logs/ | Execution logs; private and potentially identifying |
| aggregate_results/ | Allowlisted main, four-visit, targeted and prediction summaries |
| historical_falls_aggregate/ | Allowlisted historical-fall summaries and diagnostics |
| decomposition_aggregate/ | Allowlisted threshold-decomposition summaries and diagnostics |
| tables/ | Numerical sources for main and supplementary tables, plus interpretation notes |
| figures/ | Main and supplementary figures in vector/raster formats and editable sources where supported |
| numerical_comparison.json | Comparison results for the 53 reference CSVs |

The [complete output map](OUTPUT_MAP.md) connects every display to its script and source. Table S4 includes the main restrictions and weights, interval-restricted standardized risks, and paired historical-fall analyses as one sensitivity table. The reference folders in this release are read-only comparison targets, not destinations for new runs. Generated aggregates and figures still require review under provider disclosure rules before sharing.

## Verification and interpretation

Expected observed primary samples are 2,475 participants with 406 subsequent fall reports in CHARLS, and 4,779 participants with 1,134 reports in ELSA. The CHARLS four-visit observed sample contains 1,443 participants and 230 reports. Prepared eligibility extracts are larger. Assertions are specific to the recorded releases and definitions; investigate failures rather than disabling them.

The comparison uses rtol=1e-8 and atol=1e-10, trims text-cell whitespace and excludes wall-clock duration columns. It checks scientific estimates and other included diagnostics, writes a JSON report and returns a nonzero status on failure. It does not verify every manuscript cell or every methodological assumption.

The [verification record](REPRODUCTION.md) distinguishes model re-estimation, reuse of matching imputations/bootstrap caches, and export/integration checks. All 53 aggregate comparisons passed in the recorded combined export check. The code-integration check does not represent a fresh rerun of every model or all decomposition bootstraps. The historical plans in [provenance](../historical_falls_sensitivity/provenance/analysis_plan.md) and their hashes are preserved unchanged.

The primary estimates are conditional associations. Historical-fall adjustment is a paired sensitivity comparison, not causal mediation. CHARLS returning participants report falls since the preceding interview; ELSA main measures use two-year recall, and ELSA historical wave 4 uses one-year recall. Restricting interview intervals does not make these recall periods identical. The fall-or-recorded-death outcome is a composite sensitivity, not a time-to-event competing-risk model. Prediction comparison is internally cross-validated; clinical utility remains unestablished. Attenuated threshold findings and inconclusive cumulative-history results are retained.

Share only reviewed code, documentation and permitted aggregate outputs. Keep identities, contact information, local absolute paths, credentials, participant data and unredacted logs out of an anonymous repository. This local folder does not synchronize with a remote repository automatically.
