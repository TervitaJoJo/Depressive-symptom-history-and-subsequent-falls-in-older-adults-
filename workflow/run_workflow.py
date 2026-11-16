"""Run the integrated CHARLS/ELSA analysis with authorized local data."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
from runtime import external_path

CODE = Path(__file__).resolve().parents[1]
STAGES = (
    "prepare", "enrich", "main", "four", "targeted",
    "history_inputs", "history_falls", "prediction", "decomposition",
    "tables", "figures", "compare",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--start-at", choices=STAGES, default=STAGES[0])
    parser.add_argument("--stop-after", choices=STAGES, default=STAGES[-1])
    parser.add_argument("--prediction-bootstrap", type=int, default=500)
    parser.add_argument("--decomposition-bootstrap", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--dry-run", action="store_true",
                        help="List selected stages and check script paths without reading data or creating outputs.")
    args = parser.parse_args()
    if min(args.prediction_bootstrap, args.decomposition_bootstrap, args.workers) < 1:
        parser.error("Counts must be positive")
    start, end = STAGES.index(args.start_at), STAGES.index(args.stop_after)
    if end < start:
        parser.error("stop-after precedes start-at")
    raw = args.raw_root.expanduser().resolve()
    work = external_path(args.work_dir)
    prepared, extended = work / "prepared", work / "extended"
    results, cache = work / "results", work / "cache/main"
    models = CODE / "charls_elsa_models"
    history = CODE / "historical_falls_sensitivity"
    decomposition = CODE / "elsa_threshold_decomposition"
    python, rscript = [sys.executable, "-B"], [args.rscript]
    # Source preparation, estimation, export and reference checks share one entry point.
    # Each tuple names the stage, interpreter, script and arguments explicitly.
    stages = [
        ("prepare", python, CODE / "data_preparation/build_inputs.py", [raw, prepared]),
        ("enrich", python, models / "build_extended_inputs.py", [raw, prepared, extended]),
        ("main", rscript, models / "fit_models.R", [prepared, results / "main", cache]),
        ("four", rscript, models / "fit_four_visit.R", [prepared, results / "four"]),
        ("targeted", rscript, models / "targeted_analyses.R", [extended, results / "targeted", cache]),
        ("history_inputs", python, history / "build_history_falls.py", [raw, prepared, work / "history_inputs"]),
        ("history_falls", rscript, history / "fit_historical_falls.R",
         [prepared, work / "history_inputs", cache, results / "historical_falls", work / "cache/historical_falls"]),
        ("prediction", python, models / "prediction_increment.py", [extended, results / "prediction", args.prediction_bootstrap]),
        ("decomposition", rscript, decomposition / "threshold_decomposition.R",
         [prepared / "ELSA_analysis.csv", work / "cache/decomposition", results / "decomposition",
          args.decomposition_bootstrap, args.workers]),
        ("tables", python, CODE / "workflow/export_tables.py", [work]),
        ("figures", python, CODE / "workflow/render_figures.py", [work]),
        ("compare", python, CODE / "workflow/compare_results.py", [work]),
    ]
    assert tuple(stage[0] for stage in stages) == STAGES
    selected = stages[start:end + 1]
    for name, _, script, _ in selected:
        if not script.is_file():
            parser.error(f"Missing script for {name}: {script.relative_to(CODE)}")
    if args.dry_run:
        for name, _, script, _ in selected:
            print(f"{name}: {script.relative_to(CODE).as_posix()}")
        print(f"Planned {len(selected)} stages; no data read or outputs created. Prerequisites are not validated.")
        return
    if not raw.is_dir() and any(name in {"prepare", "enrich", "history_inputs"} for name, *_ in selected):
        parser.error("raw-root must be an existing authorized-data directory for the selected stages")
    work.mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    logs.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment["RUN_EXTRA_MODELS"] = "true"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    for name, interpreter, script, arguments in selected:
        print(f"Starting {name}; details are in the private log directory.", flush=True)
        command = [str(value) for value in interpreter + [script] + arguments]
        with (logs / f"{name}.log").open("w", encoding="utf8") as log:
            subprocess.run(command, cwd=work, env=environment,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        print(f"Completed {name}.", flush=True)
    if selected[-1][0] == "compare":
        print("Requested stages and aggregate reference checks completed successfully.")
    else:
        print("Requested stages completed. Run the compare stage after exporting tables to check numerical agreement.")


if __name__ == "__main__":
    main()
