# Prepared input dictionary and derivation contract

This package constructs authorized participant extracts from the specified provider-released harmonized files and raw survey modules, then re-estimates the models. See ../data_preparation/README.md for the executable source-preparation workflow. It does not redistribute participant-level data. Treat input CSVs, generated individual outputs and MICE objects as restricted local working data. Provider access permissions remain required.

## Main input files and eligibility

- `CHARLS_analysis.csv`: one first eligible window per participant; 4,062 rows before the current-symptom restriction, 2,848 with CES-D10 <10, and 2,475 with an observed subsequent fall response (406 reports).
- `ELSA_analysis.csv`: one first eligible window per participant; 6,213 rows before the current-symptom restriction, 5,489 with CES-D8 <4, and 4,779 with an observed subsequent fall response (1,134 reports).
- These files must retain eligible people with unavailable outcomes for the observation-weighting model. A dataset restricted to complete follow-up is insufficient.
- Eligibility requires age at least 65 at index, complete item-based symptom scores at the preceding and index non-proxy interviews, no recent falls at index, and cohort-specific interview-date checks described in ../data_preparation/README.md. First-window selection precedes restriction to low current symptoms; do not search later windows for a preferred symptom category. CHARLS uses the community-based study sample. ELSA excludes recorded institutional residence or unknown institutional status. The supplied analysis files already implement these choices.
- CHARLS triplets are 2011–2013–2015 and 2013–2015–2018. ELSA index waves are 5–8, each paired with its preceding and following waves. Follow-up outcomes require a valid subsequent non-proxy response. Blank/NA outcomes remain missing, not zero.

## Required fields

| Field | Coding and purpose |
|---|---|
| `id` | Unique within-cohort participant key; no duplicate analysis rows. Keep locally. |
| `cluster` | Household connected-component key. Participants sharing a recorded household across the included index waves are linked; transitive links define one variance cluster. Keep locally. |
| `window` | Categorical eligible triplet/index-wave identifier. |
| `age` | Age in years at index; natural cubic spline with 3 degrees of freedom. |
| `female` | 0 male, 1 female, as recorded in the source survey. |
| `education` | Harmonized ordered source categories 1/2/3, fitted categorically; do not equate educational credentials across countries. |
| `cesd_prior` | Complete-item score at the preceding visit: 0–30 in CHARLS, 0–8 in ELSA. |
| `cesd_current` | Complete-item score at index. Main analysis requires <10 in CHARLS or <4 in ELSA. Fitted as exact score categories. |
| `fall` | Any fall reported at follow-up: 0 no, 1 yes, missing if unavailable. Returning CHARLS participants report falls since the previous interview; ELSA main outcome waves report the previous two years. These are retrospective reports, not uniform-duration incidence or first lifetime falls. |
| `partnered` | 0 not partnered, 1 partnered under the source marital-status mapping. |
| `smoking` | 0 never, 1 former, 2 current; categorical. |
| `alcohol` | 0/1 source-defined reported alcohol use. No standardized dose is inferred. |
| `stroke`, `hearte`, `diabe`, `arthre` | 0/1 self-reported stroke, heart disease, diabetes and arthritis, respectively. |
| `stair_difficulty` | 0 no, 1 reported stair-climbing difficulty; source items differ across cohorts. |
| `adl_difficulty` | 1 if any difficulty in dressing, bathing, eating, bed transfer or toileting; 0 only if all five are observed without difficulty; otherwise missing. |
| `exposure_gap_months` | History-to-index interval in months, calculated from interview dates (days/30.4375). |
| `followup_gap_months` | Index-to-follow-up interval in months. Not used as a predictor of observation because it is unavailable for some nonresponders. |
| `survey_weight` | ELSA index-wave cross-sectional response weight. Positive and zero weights are retained in nominal N; zero weights make no likelihood contribution. No matching CHARLS weight is supplied. |
| `recurrent_fall` | Optional ELSA 0/1 indicator for at least two reported falls at follow-up; missing when unavailable. |
| `injurious_fall` | Optional ELSA 0/1 derived reported-fall/medical-treatment-injury indicator. See timing caveat below. |

Other supplied columns (`cohort`, `wave`, `group`, response/status fields) document construction; model scripts derive history from scores and observation from nonmissing outcomes. Provider-specific missing-value codes must already be converted to missing. Outcome, history and current symptom scores are not imputed. Covariates may contain NA.

The ELSA harmonized medical-treatment-injury item is routed among fall reporters. The harmonized summary uses “ever injured seriously enough to need medical treatment”; wave-specific routing and raw questionnaire wording must be retained when rebuilding it. It is analyzed as a reported falls/injury indicator and must not be relabelled as adjudicated or precisely timed incident injury. Non-fallers are coded 0 only where the questionnaire routing supports this; an unknown injury response among fallers is missing.

## Modeling and cache contract

`fit_models.R` forms `hist=1` for prior CES-D10 ≥10 or CES-D8 ≥4 and restricts the index score below the same threshold. Main-model caches must match the complete input matrix, column types, factor levels and row order. The script checks equality before fitting. Caches use 30 imputations, 20 iterations; seeds are CHARLS 20261331 and ELSA 20261230 + threshold + endpoint index (any/recurrent/injury = 1/2/3). Numeric covariates use predictive mean matching; incomplete binary/multicategory factors use logistic/multinomial regression.

Modified Poisson models use household-clustered HC1 covariance, Rubin pooling and finite-sample t intervals with complete-data degrees of freedom equal to clusters minus model parameters. The main model conditions on categorical current score, function, demographics, health characteristics, study window and both interview intervals. Standardized risks average predictions over each cohort's covariate distribution; they are adjusted observational contrasts.

Restricted analyses use the original primary imputations. The no-ADL analysis requires originally observed ADL=0 and removes the ADL term. Interval restrictions retain both intervals in 18–30 months; a constant window term is removed. Complete-case estimates use only complete observed covariates. Follow-up observation weights are stabilized by the observed fraction and estimated among all eligible current-low participants; categorical missingness is represented explicitly and numeric missingness uses medians in the observation model. The observation weighting is a sensitivity analysis, not proof that attrition bias is eliminated.

## Four-visit input

`CHARLS_four_visit_observed.csv` contains 1,443 unique participants, 230 reported falls and 1,242 clusters. Required model fields are `id`, `household`, `cluster`, `y`, `history`, `age`, `female`, `education`, `current_score`, `partnered`, `smoking`, `alcohol`, `stroke`, `hearte`, `diabe`, `arthre`, `stair_difficulty`, `adl_difficulty`, `history_gap_months`, `landmark_gap_months`, and `followup_gap_months`.

`history` is the count (0/1/2) of CES-D10 ≥10 assessments in 2011 and 2013. Participants have low symptoms (<10) and no recent falls in 2015; `current_score` is the 2015 score, and `y` is the 2018 fall report. The three intervals join 2011–2013, 2013–2015 and 2015–2018. Age spline knots are fixed from this analytic sample. Optional `primary_mi_private.rds` must match the model matrix; otherwise MICE uses seed 2026103001. The script writes participant-level audit matrices only into the authorized input directory, never into the aggregate submission folder.

The model script reproduces estimates and observed-group counts. The upstream prepare_four_visit.py now constructs the eligible extract, and workflow/export_tables.py regenerates eligible-group counts including unavailable follow-up. This overlapping subset is not an independent replication cohort.

## Aggregate result files

The supplied CSVs contain cohort-, group- or model-level summaries only. `risk0` and `risk1` refer to standardized risks under prior-low and prior-elevated symptom histories; `risk_difference` is risk1 minus risk0. Confidence limits are decimal proportions, not percentages. Non-significant results are retained. The cross-cohort fixed-effect comparison is descriptive and does not establish transportability.
