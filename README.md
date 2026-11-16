# Workflow execution reference

The [workflow overview](README_overview.md) is the installation and complete-run guide. This file documents the dependency contract for the single entry point, run_workflow.py. All paths in the output column are relative to a private work directory outside the submission package.

| Stage | Required inputs | Main outputs |
|---|---|---|
| prepare | Authorized harmonized products and survey modules | prepared/ |
| enrich | Authorized sources; prepared/ | extended/ |
| main | prepared/ | results/main/; cache/main/ |
| four | prepared/ four-visit sample | results/four/ |
| targeted | extended/; matching cache/main/ | results/targeted/ |
| history_inputs | Authorized sources; prepared/ | history_inputs/ |
| history_falls | prepared/; history_inputs/; matching cache/main/ | results/historical_falls/; cache/historical_falls/ |
| prediction | extended/ | results/prediction/ |
| decomposition | prepared/ELSA_analysis.csv | results/decomposition/; cache/decomposition/ |
| tables | All preceding estimation outputs and prepared samples | aggregate_results/; historical_falls_aggregate/; decomposition_aggregate/; tables/ |
| figures | Exported aggregate summaries | figures/ |
| compare | Exported aggregates and repository reference CSVs | numerical_comparison.json |

Stages execute in this order by default. The --start-at and --stop-after options select a contiguous range; they do not create omitted prerequisites. The --dry-run option lists that range and checks the selected script paths without reading data or creating files. It does not establish input completeness or software availability.

The runner uses the active Python interpreter with -B, enables all main-model sensitivity outputs, and sets PYTHONDONTWRITEBYTECODE=1. Rscript must be on PATH or supplied through --rscript. Each stage writes a private log and a nonzero subprocess status stops execution. compare is the default final stage, so success of a complete run includes the aggregate reference check. To resume after a failure, resolve its cause, then start from the failed stage with compatible preceding outputs.

Use a new private directory for a full run. In particular, the four-visit model refuses to overwrite a completed result directory. Cache reuse requires matching data, model matrices, methods, seeds and settings. Neither the runner nor the cache checks replace review of changed model specifications.

export_tables.py selects aggregate files explicitly from each estimation component. The results/ directory itself is not a public export. render_figures.py renders Figure 1, Figure 2 and Figure S1 from aggregate estimates. Word manuscripts and administrative submission files are outside the analysis workflow.

compare_results.py compares 53 frozen CSVs across the model, historical-fall and threshold modules with rtol=1e-8 and atol=1e-10. Text is compared after whitespace trimming; wall-clock duration columns are excluded. Differences stop the run and should be investigated without altering analyses to obtain a preferred result.

See [OUTPUT_MAP.md](OUTPUT_MAP.md) for manuscript display sources and [REPRODUCTION.md](REPRODUCTION.md) for the scope of completed verification.
