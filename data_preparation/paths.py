"""Paths and output guards; no participant data are distributed."""
from pathlib import Path
import argparse
import pyreadstat

CODE_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = CODE_ROOT.parent

def private_output(path):
    path = Path(path).expanduser().resolve()
    if path == PACKAGE_ROOT or PACKAGE_ROOT in path.parents:
        raise ValueError("Generated files must be outside the submission package.")
    path.mkdir(parents=True, exist_ok=True)
    return path

def arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("raw_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    return args.raw_root.expanduser().resolve(), private_output(args.output)

def charls_sources(root):
    harmonized = root / "CHARLS/Harmonized_CHARLS/H_CHARLS_D_Data.dta"
    raw = {}
    for year in (2011, 2013, 2015, 2018):
        folder = root / "CHARLS" / str(year)
        candidates = [f for f in folder.iterdir()
                      if f.name.lower() == "health_status_and_functioning.dta"]
        if len(candidates) != 1:
            raise ValueError(f"Expected one health file for CHARLS {year}")
        path = candidates[0]
        _, md = pyreadstat.read_dta(str(path), metadataonly=True)
        raw[f"CHARLS_{year}"] = {"path": str(path),
            "columns": [{"name": n} for n in md.column_names]}
    return {"CHARLS": {"path": str(harmonized)}}, raw
