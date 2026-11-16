from figure_style import setup_style, PALETTE, COLORS
setup_style()
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path

import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--source', default=str(Path(__file__).resolve().parent/'aggregate_results/sequential_models.csv'))
parser.add_argument('--output', default='exports/Figure_2_standardized_risks')
args = parser.parse_args()
data = pd.read_csv(args.source).query("model == 'full_adjustment'").set_index('cohort').loc[['CHARLS','ELSA']]
metrics = data.index.tolist()
low_history = (100*data.risk0).tolist()
elevated_history = (100*data.risk1).tolist()
# Native width 7 in; manuscript width 6.35 in: 10 pt becomes 9.07 pt.
plt.rcParams.update({'font.family':'Arial','font.size':10,'text.color':'#242424',
                     'axes.labelcolor':'#242424','xtick.color':'#242424','ytick.color':'#242424',
                     'pdf.fonttype':42,'svg.fonttype':'none'})

# ── 自适应高度
_fig_h = 3.25
fig, (ax, axrd) = plt.subplots(1, 2, figsize=(7, _fig_h), gridspec_kw={'width_ratios':[1.45,1]})
y = np.arange(len(metrics))

# Subtle grid
ax.grid(axis='x', alpha=0.15, linestyle='-', color=COLORS['grid'])
ax.set_axisbelow(True)

# Nondirectional connectors compare two history groups.
for i in range(len(metrics)):
    # Simple connector line
    ax.plot([low_history[i], elevated_history[i]], [y[i], y[i]],
            color='#adad9a', linewidth=2.0, solid_capstyle='round')

    # History-group estimates
    ax.scatter(low_history[i], y[i], color=PALETTE[0], s=80, zorder=3,
               edgecolors='white', linewidths=1.0, label='Prior low' if i == 0 else '')
    ax.scatter(elevated_history[i], y[i], color=PALETTE[2], s=80, zorder=3,
               edgecolors='white', linewidths=1.0, label='Prior elevated' if i == 0 else '')

    # Display exact endpoint values; no invented significance stars or endpoint CIs.
    ax.text(low_history[i], y[i]-0.23, f'{low_history[i]:.1f}%', ha='center', fontsize=9, color=PALETTE[0])
    ax.text(elevated_history[i], y[i]+0.28, f'{elevated_history[i]:.1f}%', ha='center', fontsize=9, color=PALETTE[2])
    row = data.iloc[i]
    rd, low, high = 100*row.risk_difference, 100*row.rd_low, 100*row.rd_high
    axrd.errorbar(rd, y[i], xerr=[[rd-low],[high-rd]], fmt='D', ms=6,
                  color=PALETTE[0], ecolor=PALETTE[0], capsize=4, lw=1.6)
    axrd.text(rd, y[i]-0.24, f'{rd:.1f} ({low:.1f} to {high:.1f})', ha='center', fontsize=9)

ax.set_yticks(y)
ax.set_yticklabels(metrics, fontsize=10)
ax.set_xlabel('Standardized fall risk (%)', fontsize=10)
ax.set_xlim(0,40)
ax.set_xticks([0,10,20,30,40])
ax.set_ylim(-0.6,1.6)
ax.legend(loc='upper center', bbox_to_anchor=(0.55,-0.28), ncol=2, frameon=False, fontsize=9, columnspacing=1)
ax.invert_yaxis()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
axrd.set_xlim(0,12)
axrd.set_xticks([0,4,8,12])
axrd.set_ylim(1.6,-0.6)
axrd.set_yticks([])
axrd.set_xlabel('Risk difference (percentage points)', fontsize=9)
axrd.axvline(0,color='#777777',lw=0.8,ls='--')
axrd.spines[['top','right','left']].set_visible(False)
axrd.grid(axis='x',alpha=0.15)
ax.set_title('A  Adjusted absolute risks', loc='left', fontsize=11, fontweight='bold', pad=15)
axrd.set_title('B  Difference and 95% CI', loc='left', fontsize=11, fontweight='bold', pad=15)
fig.subplots_adjust(left=0.12,right=0.98,bottom=0.31,top=0.83,wspace=0.3)
out=external_path(args.output)
out.parent.mkdir(parents=True,exist_ok=True)
for ext in ['png','pdf','svg']:
    fig.savefig(out.with_suffix('.'+ext),dpi=600,facecolor='white')
fig.savefig(out.with_suffix('.tiff'),dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print('Created',out)

