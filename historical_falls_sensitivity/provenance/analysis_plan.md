# Post hoc analysis addendum: historical falls and interval-restricted absolute risk

Locked on 2026-11-15, before extracting historical falls or fitting the new association models. This addendum follows knowledge of the original primary and supplementary results and is not a preregistration.

## Locked primary analysis
Retain the original outcome, thresholds, first-eligible-window selection before the low-current-score restriction, current-age criterion, covariates, household clustering, and published primary estimates. Do not select new windows or thresholds. All new results are exploratory and will be reported irrespective of significance.

## Historical falls
Retrieve the recorded binary falls response at the assessment that supplied the historical CES-D score. CHARLS history waves are 2011 and 2013; ELSA history waves are 4-7. Preserve original IDs, row order and selected windows. Do not reconstruct eligibility from the added field. Document source questions and recall periods.

Report missing historical falls, available-sample N/events/history groups, and the historical-falls distribution. Participants with missing historical falls are excluded from this sensitivity only; history falls is not imputed. Fit the original full model and the full model plus historical falls in exactly the same available-history sample and the same completed data sets. Use 30 chained-equation imputations, 20 iterations, with historical falls included as a fully observed auxiliary predictor. This paired imputation model avoids attributing a change between models to a change in sample or imputed values. Record fixed seeds 202611151 (CHARLS), 202611152 (ELSA). Retain the original factor levels, predictors, modified-Poisson log link, HC1 household-cluster covariance and finite-sample Rubin pooling. Report paired models, their point log-RR difference and the RR ratio descriptively; do not call this a causal decomposition. Separately reproduce the frozen primary model using the original matching imputation objects.

Restriction to no falls at both the historical and index assessments will be included only if, before fitting, each symptom-history group has at least 100 participants and at least 20 follow-up falls, with at least 200 total falls and 40 households. These fixed capacity rules are not a power calculation. Use the same completed data sets restricted to observed historical nonfallers. Do not call these participants lifetime nonfallers.

## Interval-restricted absolute risks
Use the existing restriction that both adjacent interview intervals lie within 18-30 months. Reuse the frozen primary imputations and model specification. Estimate standardized risks for historical low and elevated scores over the restricted sample, and the adjusted risk difference with the existing delta-method/Rubin CI. Verify N/events and RR against the existing restricted-analysis output. Do not compare absolute risks across cohorts as a common annual incidence.

## Diagnostics and deliverables
Check data linkage, uniqueness, missingness, imputation events/convergence, model convergence, finite covariance, fitted probabilities, exact primary-result invariance and restricted-model reproduction. Save person-level data/caches/logs only in private work directories. Release only anonymous scripts, the dated plan, aggregate outputs, source-question mapping and reproducibility notes. Add findings to an explicitly post hoc section and supplementary table S4; retain all prior sensitivity results. No new prediction, machine learning, mediation, or time-to-event analysis is planned.
