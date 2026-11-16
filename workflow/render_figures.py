"""Render the two main figures and one supplementary figure."""
from pathlib import Path
import subprocess,sys
from runtime import external_path
if len(sys.argv)!=2:raise SystemExit("Usage: python render_figures.py PRIVATE_WORK_DIRECTORY")
w=external_path(sys.argv[1]);m=Path(__file__).resolve().parents[1]/"charls_elsa_models"
a=w/"aggregate_results";f=w/"figures";f.mkdir(exist_ok=True)
for cmd in [
 [m/"gen_fig_flow.py","--source",a/"cohort_flow.csv","--output",f/"Figure_1_design_flow"],
 [m/"gen_fig_risks.py","--source",a/"sequential_models.csv","--output",f/"Figure_2_standardized_risks"],
 [m/"gen_fig_score_shape.py",a,f]]:
 subprocess.run([sys.executable,"-B",*map(str,cmd)],check=True)
