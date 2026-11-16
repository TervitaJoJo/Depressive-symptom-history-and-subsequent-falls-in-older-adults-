# ELSA threshold decomposition: post hoc analysis plan

> Editorial numbering note (2026-11-02): the decomposition display originally called Table S8 in this saved post hoc plan is Table S1 in the current manuscript. The analytical plan below is unchanged; only its display number has changed to follow first citation order.

Created 2026-11-01, after existing primary/sensitivity results and before fitting B/C or the common-pipeline decomposition. Explicitly authorized by the user. This is not preregistration and does not replace any existing primary analysis.

## Question and fixed scope
Among the existing first eligible ELSA windows, how does the fitted historical-symptom association change when current-score eligibility and historical-score dichotomization are changed separately? Keep observed follow-up falls, age >=65, no index falls, person/window selection, and original full covariate adjustment unchanged.

A: current CES-D8 <4, historical CES-D8 >=4.
B: current CES-D8 <3, historical CES-D8 >=4.
C: current CES-D8 <4, historical CES-D8 >=3.
D: current CES-D8 <3, historical CES-D8 >=3.

## Controlled estimation
Use the broad A/C population with observed falls (expected N=4779). One common 30-imputation,20-iteration MICE system, seed 20261101, includes exact observed historical score as a categorical predictor replacing the original binary history variable; derive both history variables afterward. Current score is categorical. Use PMM/logistic/polytomous regression for incomplete covariates as in original analysis. Do not impute outcomes or symptom scores. Keep all four fits on these same completed covariate data, with age spline knots/boundaries fixed from the broad sample. The narrower analyses subset completed data without reselecting windows or rerunning imputation.

Original full model: modified Poisson log link; history + 3-df natural spline of age + sex + education + window + both interview intervals + categorical current score + partnership + smoking + alcohol + stroke + heart disease + diabetes + arthritis + stair difficulty + ADL difficulty. Current score3 is absent from B/D and its column is removed.

## Estimands and uncertainty
Let a,b,c,d be logRRs. Total=d-a; eligibility=((b-a)+(d-c))/2; history-definition=((c-a)+(d-b))/2; path interaction=d-b-c+a. Report both ordered paths as well as the average. This is an algebraic, order-averaged allocation of changes in fitted associations, not a causal decomposition or an invariant percentage attribution.

Primary uncertainty: MI Boot as in Schomaker & Heumann (2018), extended to a four-estimate vector. Within each completed dataset, draw 1000 household bootstrap samples from all broad-sample households, pair each draw across all four cells, use bootstrap covariance as within-imputation uncertainty, and combine it with between-imputation covariance using multivariate Rubin rules. Point estimates come from unbootstrapped completed data. Scalar contrasts use t intervals with Barnard-Rubin degrees of freedom and conservative complete-data df=minimum(households-model parameters). Independently cross-check against joint household HC1 sandwich covariance. Independent bootstrap streams across imputations, paired within each four-cell fit. No percentile pooling of imputed estimates.

Report all fits and failed resamples. If nonconvergence or singularity occurs, retain audit and assess validity before inference; do not silently filter failures. Check decomposition identity, bootstrap vs sandwich precision, MI diagnostics, and original-vs-controlled endpoint estimates. Original A/D are not overwritten. Differences from original pipeline are reported separately and are not assigned to either rule change.

## Deliverables and boundaries
Add Table S8 and a reproducible code/aggregate-results folder. Limited replacement/addition in Methods, Results, Discussion and cover letter. Preserve main two tables/two figures, outcomes, original main statistics and contradictory evidence. No data-driven new threshold, no component percentage headline, no causal interpretation, no claim of threshold invariance. Private data and caches stay outside the submission package.
