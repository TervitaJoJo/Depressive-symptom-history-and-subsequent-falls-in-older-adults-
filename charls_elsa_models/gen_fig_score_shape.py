"""Adapted from vivid basic.line, using fitted historical-score associations.
Final width 6.3 in / native width 6.8 in = 0.926;
11-point labels render at 10.2 points. No outcome-derived cutoff selection.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from figure_style import setup_style, PALETTE, COLORS
setup_style()
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 10})
inp = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent/'aggregate_results'
out = external_path(sys.argv[2] if len(sys.argv)>2 else Path.cwd())
out.mkdir(parents=True, exist_ok=True)
curves = pd.read_csv(inp/'spline_curves.csv')
assert set(curves.cohort)=={'CHARLS','ELSA'}

fig = plt.figure(figsize=(6.8, 4.8))
gs = fig.add_gridspec(2, 2, height_ratios=[3.8, 1], hspace=0.10, wspace=0.23)
ymin = min(0.8, float(curves.ci_low.min())*0.94)
ymax = float(curves.ci_high.max())*1.06
for i, (co, cutoff, scale) in enumerate([('CHARLS',10, 'CES-D10'), ('ELSA',4, 'CES-D8')]):
    ax = fig.add_subplot(gs[0,i])
    histax = fig.add_subplot(gs[1,i], sharex=ax)
    d = curves.loc[curves.cohort.eq(co)]
    x = d.score.to_numpy()
    y = d.rr.to_numpy()
    color = PALETTE[0 if co=='CHARLS' else 2]
    # Template background, one genuine 95% confidence band, fitted line/markers.
    # Arbitrary demo bands and best-method extrema do not apply to this analysis.
    ax.axhspan(ymin, ymax, alpha=0.02, color=color, zorder=0)
    ax.fill_between(x, d.ci_low.to_numpy(), d.ci_high.to_numpy(),
                    alpha=0.15, color=color, linewidth=0)
    ax.plot(x, y, 'o-', color=color, linewidth=1.8, markersize=5,
            markevery=30, markeredgecolor='white', markeredgewidth=1.2, zorder=3)
    ax.axhline(1, color=COLORS['ref_line'], linewidth=1, linestyle='--')
    ax.axvline(cutoff, color=COLORS['gray'], linewidth=1, linestyle=':')
    ax.set_title(f"{'A' if i==0 else 'B'}  {co}", loc='left', fontsize=12, fontweight='bold', pad=9)
    ax.set_ylim(ymin,ymax)
    ax.set_xlim(-.45, x.max()+.45)
    if i==0: ax.set_ylabel('Adjusted risk ratio', fontsize=11)
    else: ax.tick_params(labelleft=False)
    ax.tick_params(labelbottom=False, labelsize=9)
    ax.grid(alpha=0.12, linestyle='--', color=COLORS['grid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    counts = pd.read_csv(inp/f'{co}_score_counts.csv')
    histax.bar(counts.score, counts.Freq, width=0.78, color=color, alpha=0.65, edgecolor='white', linewidth=0.3)
    histax.axvline(cutoff, color=COLORS['gray'], linewidth=1, linestyle=':')
    histax.set_xlabel(f'Prior {scale} score', fontsize=11)
    if i==0: histax.set_ylabel('Count', fontsize=10)
    histax.tick_params(labelsize=9)
    histax.set_xticks([0,5,10,15,20,25] if co=='CHARLS' else np.arange(0,9))
    histax.set_xlim(-.45,x.max()+.45)
    histax.spines['top'].set_visible(False)
    histax.spines['right'].set_visible(False)
fig.legend(handles=[Line2D([0],[0],color=COLORS['dark'], marker='o', markersize=4, label='Spline estimate'),
                    Patch(facecolor=COLORS['gray'],alpha=.4,label='95% CI'),
                    Line2D([0],[0],color=COLORS['gray'],linestyle=':',label='Primary cutoff')],
           loc='lower center', ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(.51,.01))
fig.subplots_adjust(left=.11,right=.98,top=.91,bottom=.23)
for ext in ['png','pdf','tiff']:
    kw={'dpi':600, 'facecolor':'white'}
    if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
    if ext=='pdf':kw['metadata']={'Creator':'','Author':''}
    fig.savefig(out/f'Figure_S1_history_score_shape.{ext}',**kw)
plt.close(fig)

