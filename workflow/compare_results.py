"""Compare regenerated aggregate CSVs with frozen publication references."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from runtime import external_path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('work_directory',type=Path);p.add_argument('--rtol',type=float,default=1e-8);p.add_argument('--atol',type=float,default=1e-10)
a=p.parse_args();w=external_path(a.work_directory);code=Path(__file__).resolve().parents[1];checks=[]
for project,generated in [('charls_elsa_models','aggregate_results'),('elsa_threshold_decomposition','decomposition_aggregate'),('historical_falls_sensitivity','historical_falls_aggregate')]:
 for reference in sorted((code/project/'aggregate_results').glob('*.csv')):
  target=w/generated/reference.name
  result=dict(project=project,file=reference.name,passed=False)
  if not target.exists():result['reason']='Missing regenerated file';checks.append(result);continue
  try:
   try:left=pd.read_csv(reference)
   except pd.errors.EmptyDataError:left=pd.DataFrame()
   try:right=pd.read_csv(target)
   except pd.errors.EmptyDataError:right=pd.DataFrame()
   if list(left.columns)!=list(right.columns):
    if set(left.columns)!=set(right.columns):raise AssertionError('Column names differ')
    right=right[left.columns]
   assert len(left)==len(right),'Row count differs'
   # Wall-clock timings are diagnostic metadata, not reproducibility targets.
   ignored=[c for c in left if c in {'seconds','elapsed_seconds'}]
   result['ignored_runtime_columns']=ignored
   left=left.drop(columns=ignored);right=right.drop(columns=ignored)
   numeric=[c for c in left if pd.api.types.is_numeric_dtype(left[c]) and pd.api.types.is_numeric_dtype(right[c])]
   labels=[c for c in left if c not in numeric]
   for c in labels:
    left[c]=left[c].astype('string').str.strip();right[c]=right[c].astype('string').str.strip()
   if labels:left=left.sort_values(labels,na_position='last').reset_index(drop=True);right=right.sort_values(labels,na_position='last').reset_index(drop=True)
   pd.testing.assert_frame_equal(left[labels],right[labels],check_dtype=False)
   diffs={}
   for c in numeric:
    l=left[c].to_numpy(float);r=right[c].to_numpy(float)
    assert np.allclose(l,r,rtol=a.rtol,atol=a.atol,equal_nan=True),f'Numeric mismatch: {c}'
    diff=np.abs(l-r);diff=diff[np.isfinite(diff)];diffs[c]=float(diff.max()) if len(diff) else 0.0
   result.update(passed=True,rows=len(left),numeric_columns=len(numeric),max_absolute_difference=max(diffs.values(),default=0))
  except Exception as exc:result['reason']=str(exc)
  checks.append(result)
report=dict(passed=all(c['passed'] for c in checks),rtol=a.rtol,atol=a.atol,checks=checks)
(w/'numerical_comparison.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps({'passed':report['passed'],'tables_compared':len(checks),'failures':[c for c in checks if not c['passed']]},indent=2))
raise SystemExit(0 if report['passed'] else 1)
