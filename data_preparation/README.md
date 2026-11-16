# Licensed input preparation

Start with the [root README](../workflow/README_overview.md). This module constructs every prepared input used by the delivered two-cohort analysis. No participant data are included.

## Required layout under the supplied RAW_ROOT

```text
CHARLS/
  Harmonized_CHARLS/H_CHARLS_D_Data.dta
  2011/Health_Status_and_Functioning.dta
  2013/Health_Status_and_Functioning.dta
  2015/Health_Status_and_Functioning.dta
  2018/Health_Status_and_Functioning.dta
ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/
  h_elsa_g3.dta
  wave6/wave_6_elsa_nurse_data_v2.dta
```

CHARLS health filenames are resolved case-insensitively by the initial preparation module. The enrichment script uses the displayed capitalization; on a case-sensitive system, use the displayed spelling. These are the provider-release filenames used in this analysis. Keep accompanying provider documentation and access licenses locally; no access authorization is granted by this code.

Run from the submission-package directory:

```text
python -B code/data_preparation/build_inputs.py "../licensed_data" "../reproduction_private/prepared"
```

All three constructors run sequentially. `paths.py` rejects output inside the submission package. The existing model scripts consume the generated `CHARLS_analysis.csv`, `ELSA_analysis.csv`, and `CHARLS_four_visit_observed.csv`. The four-visit eligible extract is retained for missing-follow-up and group summaries.

## Construction sequence

`prepare_charls.py` reads Harmonized CHARLS D and the 2011, 2013, 2015 and 2018 health modules. It recalculates complete CES-D10 scores, reverses the two positive items, verifies agreement with the harmonized totals, aligns baseline IDs through ID_w1, and derives self-response and falls from the raw health records. The index windows are 2011–2013–2015 and 2013–2015–2018. Age, community residence, no index falls, complete scores and ordered dates determine eligibility before selecting the first window. Restriction to low current symptoms occurs later. Household links across index waves define clusters.

`prepare_elsa.py` reads Harmonized ELSA G. It validates complete CES-D8 items and reversed positive items for each index wave 5–8 and its preceding wave. It selects each person's first eligible window before applying the original recorded-institution/unknown-status filter and before the current-score restriction. It then constructs connected household clusters using index-wave household identifiers and retains the index response weight. It preserves the original order of these operations and verifies positive intervals in the final observed primary sample; it does not search later windows after exclusions. No-fall responses map to zero for secondary outcomes; missing information among fall reporters remains missing.

`prepare_four_visit.py` reconstructs the fixed-2015 CHARLS index analysis with symptom assessments in 2011, 2013 and 2015 and follow-up in 2018. It retains both the eligible and observed extracts, counts 0/1/2 elevated historical assessments, verifies score recalculation and date ordering, and checks common index-covariate coding against the primary CHARLS extract.

Expected primary inputs: CHARLS 4,062 eligible rows, 2,848 currently low, 2,475 observed, 406 fall reports; ELSA 6,213 eligible rows, 5,489 currently low, 4,779 observed, 1,134 fall reports. The four-visit extract has 1,718 eligible and 1,443 observed participants with 230 fall reports. Assertions deliberately check this manuscript's releases and definitions; a failed assertion must be investigated rather than disabled.

The covariate mappings and outcome cautions are documented in `../charls_elsa_models/DATA_DICTIONARY.md`. Directed medication, recall, vision, walking, balance and recorded-death enrichment is implemented in `../charls_elsa_models/build_extended_inputs.py` and documented in `../charls_elsa_models/EXTENDED_DATA_DICTIONARY.md`. Historical falls are linked to these same selected windows by `../historical_falls_sensitivity/build_history_falls.py`; see its module guide for wave-specific recall periods.

Generated participant identifiers are necessary to link visits and account for household dependence. They remain exclusively in the caller's private output area and are not anonymous shareable data. Only aggregate results belong in a public code release.
