# Manuscript output map

The tables stage exports numerical source files from the complete analysis. It does not rebuild the Word manuscript. Paths in the last column are relative to the private work directory. The three aggregate folders have distinct roles: aggregate_results/ for main and targeted models, historical_falls_aggregate/ for historical-fall sensitivity, and decomposition_aggregate/ for ELSA threshold decomposition.

| Display | Computation | Generated sources |
|---|---|---|
| Table 1 | Observed prepared samples, variable-specific denominators | tables/Table_1.csv; tables/Table_1_numeric.csv |
| Table 2 | charls_elsa_models/fit_models.R | tables/Table_2_source.csv; aggregate_results/sequential_models.csv |
| Figure 1 | data_preparation/build_inputs.py; charls_elsa_models/gen_fig_flow.py | aggregate_results/cohort_flow.csv; figures/ |
| Figure 2 | charls_elsa_models/fit_models.R; gen_fig_risks.py | aggregate_results/sequential_models.csv, fully adjusted rows; figures/ |
| Table S1 | elsa_threshold_decomposition/threshold_decomposition.R | tables/Table_S1_cells.csv; tables/Table_S1_contrasts.csv |
| Table S2 | Observed samples and observation summaries | tables/Table_S2_missingness.csv; tables/Table_S2_observation.csv |
| Table S3 | Joint-information tests and descriptive cross-cohort comparison | tables/Table_S3_joint_tests.csv; tables/Table_S3_descriptive_comparison.csv |
| Table S4A and secondary-outcome rows | Main restrictions, weights and ELSA secondary outcomes | tables/Table_S4_restrictions.csv; Table_S4_observation_weights.csv; Table_S4_response_weights.csv; Table_S4_secondary_outcomes.csv (all under tables/) |
| Table S4B | historical_falls_sensitivity/fit_historical_falls.R, both_intervals_18_30 | tables/Table_S4B_interval_absolute_risks.csv |
| Table S4C | historical_falls_sensitivity/fit_historical_falls.R, paired models and two-assessment restriction | tables/Table_S4C_historical_falls.csv; tables/Table_S4C_sample_and_missingness.csv |
| Table S5 | charls_elsa_models/targeted_analyses.R | tables/Table_S5_directed_sensitivity.csv, medication/domain rows; tables/Table_S5_evalues.csv |
| Table S6 | charls_elsa_models/targeted_analyses.R | tables/Table_S6_mortality.csv; tables/Table_S6_composite_models.csv, restricted/composite rows |
| Table S7 | charls_elsa_models/targeted_analyses.R; prediction_increment.py | tables/Table_S7_continuous_score.csv; tables/Table_S7_spline_tests.csv; tables/Table_S7_prediction.csv |
| Table S8 | charls_elsa_models/fit_four_visit.R; eligible/observed summaries | tables/Table_S8_counts.csv; tables/Table_S8_models.csv; tables/Table_S8_joint_test.csv |
| Figure S1 | charls_elsa_models/gen_fig_score_shape.py | aggregate_results/spline_curves.csv; cohort score-count CSVs in the same folder; figures/ |

The historical-fall source is historical_falls_aggregate/historical_falls_models.csv. capacity_and_missingness.csv, sample_counts.csv and linkage_summary.csv describe available samples and linkage; model_diagnostics.csv and the imputation CSVs describe numerical behavior. locked_result_verification.csv verifies the fixed primary and interval-restricted results. These diagnostics support Table S4 and the methods in Online Resource 1.

Standardized risks describe adjusted associations. Interval restriction does not unify cohort recall periods. The paired historical-fall log-RR change is descriptive, not causal attribution. The four-visit sample overlaps the primary CHARLS sample. Prediction estimates use internal household-grouped cross-validation. The fall-or-recorded-death outcome is a composite sensitivity, not a time-to-event competing-risk model. Threshold decomposition remains exploratory; both switching orders and uncertainty are retained.

The tables/README.json generated during export identifies row-selection and interpretation details. All 53 reference CSVs are checked after the complete workflow; only wall-clock duration columns are excluded from the numerical comparison.
