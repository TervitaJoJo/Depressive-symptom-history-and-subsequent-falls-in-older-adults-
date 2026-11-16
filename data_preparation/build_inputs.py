"""Build all prepared input files required by this two-cohort analysis."""
from pathlib import Path
import subprocess
import sys
import json
import pandas as pd
from paths import arguments

root, output = arguments("Build CHARLS and ELSA inputs from licensed releases.")
for script in ("prepare_charls.py", "prepare_elsa.py", "prepare_four_visit.py"):
    subprocess.run([sys.executable, "-B", str(Path(__file__).with_name(script)),
                    str(root), str(output)], check=True)
flow=[]
for cohort, cutoff in (("CHARLS",10),("ELSA",4)):
    d=pd.read_csv(output/f"{cohort}_analysis.csv")
    low=d.cesd_current.lt(cutoff); observed=low & d.fall.notna()
    flow.append(dict(cohort=cohort, eligible_all=len(d), low_current=int(low.sum()),
        observed=int(observed.sum()), falls=int(d.loc[observed,"fall"].sum()),
        prior_low=int((observed & d.cesd_prior.lt(cutoff)).sum()),
        prior_elevated=int((observed & d.cesd_prior.ge(cutoff)).sum())))
pd.DataFrame(flow).to_csv(output/"cohort_flow.csv",index=False)
print(json.dumps({"completed":True,"cohort_flow":flow},indent=2))
