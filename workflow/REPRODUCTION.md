# Reproduction checks and boundaries

Computational checks completed on 2026-11-01; release review completed on 2026-11-02. This record describes verification of the anonymous code release; it is not a preregistration or a claim of a new independent study.

## Input construction

The workflow was executed locally using authorized provider-released harmonized CHARLS and ELSA files and the raw modules specified in `../data_preparation/README.md`. It does not reconstruct the data providers' harmonized products. Neither input data nor participant-level derived files are included in this release.

Freshly reconstructed extracts were compared with the original analysis extracts using exact value, row-order and column-order comparisons:

| Extract | Rows | Columns | Check |
|---|---:|---:|---|
| CHARLS main analysis input | 4,062 | 24 | Exact match |
| ELSA main analysis input | 6,213 | 27 | Exact match |
| CHARLS four-visit eligible sample | 1,718 | 30 | Exact match |
| CHARLS four-visit observed sample | 1,443 | 30 | Exact match |

Extended inputs for the targeted analyses were also reconstructed and matched exactly. The observed primary samples contained 2,475 participants and 406 fall reports in CHARLS, and 4,779 participants and 1,134 fall reports in ELSA. Eligible input extracts include participants outside these final observed samples.

## Execution and numerical checks

- Main models, supplementary restrictions, weighting analyses, alternative ELSA thresholds and secondary outcomes were re-estimated.
- The CHARLS four-visit models and targeted supplementary models were re-estimated.
- Prediction comparisons were rerun, including 500 bootstrap replicates per cohort.
- The ELSA threshold-decomposition stage was executed using the matching original imputation objects and paired bootstrap caches; its pooled results and diagnostics were regenerated.
- Table-source export and rendering of both main figures and the supplementary figure completed.
- All 42 aggregate CSV files passed comparison with the frozen reference outputs using numeric tolerances `rtol=1e-8` and `atol=1e-10`; text columns were compared after whitespace normalization. The 41 previously supplied aggregate CSV files are unchanged, and one aggregate sample-flow file was added.

The main, four-visit and targeted checks reused matching original multiple-imputation objects. The decomposition reused 30 completed data sets and their 1,000 paired household-bootstrap replicates per imputation. Consequently, this verification is **not** a fresh regeneration of every imputation and all 120,000 decomposition model fits. The scripts support a fresh run without caches, but that full uncached run was not performed for this code-release check. The prediction bootstraps were freshly run.

## Repeating the checks

Follow the root `../README.md` to run from licensed inputs in a new private working directory. Then run `compare_results.py` against that directory. Caches are deliberately not distributed. Fresh runs require substantial computation and may differ numerically across software versions or platforms; inspect any differences before interpreting new results. Do not change the analysis to force statistical significance or reference agreement.

The checked software versions are recorded in `../charls_elsa_models/ENVIRONMENT.json` and the aggregate session information. Some installed R packages reported that they were built under a later patch version; execution completed and the aggregate comparisons passed in the checked environment. This is not a guarantee of identical results on every installation.

## Release contents

Only CHARLS and ELSA code, documentation and aggregate reference outputs are included. Participant identifiers, restricted data, derived participant records, imputations, individual predictions, bootstrap caches, unredacted run logs, production scripts and private audit reports remain outside the release. Output-path guards require generated files to stay outside the submission package.

## Historical-fall verification, 15 November 2026

Original primary N/events/RRs/standardized risks/RDs and existing interval-restriction N/events/RRs passed numeric reproduction at tolerance 1e-8. The historical-fall sensitivity used fresh 30-by-20 imputation per cohort; all 300 completed-data model fits (including reproduction checks) converged, with no logged imputation events. Model and imputation diagnostics are included. No original threshold, primary outcome or selected window was changed. The complete original computationally expensive prediction/decomposition suite was not rerun for a manuscript-format revision; their previously verified frozen outputs are retained.

## Integrated release verification, 15 November 2026

The complete release organizes analyses by purpose and uses one runner for preparation, estimation, table export, figure rendering and comparison. Historical-fall sensitivity is part of the default sequence before prediction and threshold decomposition. The documented plan and source amendment have been relocated to `../historical_falls_sensitivity/provenance/` without changing their contents or hashes; the timing of exploratory analyses is unchanged.

This integration check verifies Python/R syntax, stage order and script paths, partial-run selection, failure propagation, output-path guards, local documentation links and file inventories. The table export is executed using the existing model outputs and historical-fall results, and all 53 aggregate CSV comparisons are checked. All frozen aggregate files are byte-identical to the preceding release. The integration does not rerun statistical models, prediction bootstraps or decomposition imputations/bootstrap fits. It does not alter the manuscript, supplements or figures.

The code directory contains 97 files in total. Every directory, including the root, has at most 100 files recursively. The anonymous release excludes participant records, caches, private audit reports, author details, personal paths and credentials.
