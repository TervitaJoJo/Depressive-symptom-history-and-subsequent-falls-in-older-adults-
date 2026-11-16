# Extended input contract

The base input contract remains in DATA_DICTIONARY.md. Enrichment preserves its first-window selection and original rows. Identifier values are not supplied. `wave` identifies the index wave. Added variables are not part of the frozen primary imputation matrix.

| Field | Definition and provenance |
|---|---|
| recorded_death | Next-wave harmonized `r{w+1}iwstat` in {5,6}; other/unknown status does not mean confirmed survival. |
| target_status | Unmodified harmonized next-wave interview status. |
| recall_total | Index immediate plus delayed word recall, each restricted to 0–10; missing if either missing. |
| poor_eyesight | Fair/poor near or distance vision; CHARLS DA033/DA034 ≥4, ELSA dsight/nsight ≥4 including blind6. One observed impaired measure suffices for1; both must be observed and unimpaired for0. |
| log_walk_time | Natural log of positive harmonized `r{w}wspeed`, representing mean walking-test time in seconds, not metres/second. Protocols differ; cohorts are modeled separately. |
| balance_score | CHARLS index `r{w}balance`, categories1–4, modeled as a factor. Not used in ELSA. |
| antidepressant_report | CHARLS DA020S2 =2 at either preceding or index wave. A zero means no report, not confirmed nonuse. |
| sedative_report | CHARLS DA020S3 =3 (tranquilizers/sleeping pills) at either wave. |
| medication_item_completed | Any DA020 response at both assessments; descriptive coverage flag, not used to assert nonuse. |
| elsa_prior_med_observed | Index wave7 only: wave6 nurse date strictly precedes index month, plus explicit no medicines or at least one classifiable six-digit drug reported taken in previous7days. |
| elsa_prior_antidepressant | Any prior nurse drug code beginning0403 with MedBIA/MedBIAj=1, among observed inventories. Missing otherwise. |

CHARLS raw2011 IDs use the harmonized ID_w1 crosswalk. ELSA DrC1–DrC27 drug classes follow the supplied Wave6 Drug Codebook (04.03 antidepressants). Full wave8 inventory classes require different access and are not inferred from available variables. Incomplete lists and use outside the nurse snapshot remain unresolved.

Memory, vision and walking are observed-domain subset restrictions; CHARLS additionally requires balance. Both original and augmented models use the same newly imputed original covariates within each subset. Historical-score models include the observed score and spline basis in new covariate imputations. The composite outcome includes valid fall respondents plus recorded deaths; missing other outcomes are never set to zero. All added outcome analyses are exploratory.

Aggregate outputs: targeted_models summarizes14 fits; targeted_diagnostics summarizes 18 model types (including linear and spline score models), each across30 imputations; continuous_score_models has2 linear contrasts; spline_tests4 tests; spline_curves302 grid estimates; spline_knots2 knot sets; score_counts are frequency tables; mortality_ascertainment has16 status cells; domain_missingness7 counts; evalues2 rows; prediction_increment18 metric estimates. None is person-level data.
