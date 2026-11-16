"""Deterministic design/selection diagram, plus editable Draw.io geometry.
Native width 7 in, delivered width 6.35 in; 10 pt -> 9.07 pt.
The installed host has no Draw.io/TikZ export executable. Matplotlib renders
these same explicit nodes/edges to vector PDF/SVG; Draw.io source is retained.
"""
from pathlib import Path
import argparse, json
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path

import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from figure_style import setup_style, PALETTE
setup_style()
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
p=argparse.ArgumentParser();p.add_argument('--output',default='exports/Figure_1_design_flow');p.add_argument('--source',default=str(Path(__file__).resolve().parent/'aggregate_results/cohort_flow.csv'));a=p.parse_args()
out=external_path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
fig,ax=plt.subplots(figsize=(7,5.8));ax.set_xlim(0,7);ax.set_ylim(0,5.8);ax.axis('off')
color=PALETTE[0];txt='#202020';nodes=[];edges=[];labels=[]
def text(x,y,t,size=9,bold=False,ha='center'):
    ax.text(x,y,t,ha=ha,va='center',fontsize=size,color=txt,fontweight='bold' if bold else 'normal')
    labels.append((x,y,t,size,bold,ha))
def box(key,x,y,w,h,t,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=0.03',facecolor='white',edgecolor=color,lw=1.0))
    text(x+w/2,y+h/2,t,9,bold);nodes.append((key,x,y,w,h,t))
def arrow(x1,y1,x2,y2):
    ax.annotate('',xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle='-|>',color=color,lw=1.1));edges.append((x1,y1,x2,y2))
text(.08,5.63,'A  Temporal design',11,True,'left')
box('prior',.08,4.66,2.03,.66,'History assessment\nPrior elevated or low symptoms')
box('index',2.48,4.66,2.03,.66,'Index assessment\nLow symptoms; no recent falls')
box('follow',4.88,4.66,2.03,.66,'Following assessment\nSelf-reported falls')
arrow(2.13,4.99,2.43,4.99);arrow(4.53,4.99,4.83,4.99)
text(3.5,4.34,'CHARLS: 2011 → 2013 → 2015 or 2013 → 2015 → 2018\nELSA: index waves 5–8, each with preceding and following waves',8.8)
text(.08,3.85,'B  Participant selection',11,True,'left')
text(3.5,3.51,'Age ≥65; complete symptom scores at two non-proxy interviews;\nno recent falls at index; one first eligible window per participant',8.8)
flow=pd.read_csv(a.source).set_index('cohort')
for x,co in [(.12,'CHARLS'),(3.73,'ELSA')]:
    row=flow.loc[co]
    n,lo,obs,events,hist=(int(row[k]) for k in ['eligible_all','low_current','observed','falls','prior_elevated'])
    hi,miss=n-lo,lo-obs
    w=3.13
    text(x+w/2,3.10,co,11,True)
    box(co+'eligible',x,2.44,w,.46,f'First eligible participants: {n:,}')
    arrow(x+.08,2.41,x+.08,2.04)
    text(x+w/2,2.23,f'Current elevated symptoms excluded: {hi:,}',8.0)
    box(co+'low',x,1.52,w,.46,f'Low current symptoms: {lo:,}')
    arrow(x+.08,1.49,x+.08,1.12)
    text(x+w/2,1.31,f'Follow-up fall response unavailable: {miss:,}',8.0)
    box(co+'observed',x,.41,w,.64,f'Primary sample: {obs:,}\nSubsequent fall reports: {events:,}',True)
    text(x+w/2,.16,f'Prior elevated symptoms: {hist:,}',8.4)
fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01)
for ext in ['png','pdf','svg']:fig.savefig(out.with_suffix('.'+ext),dpi=600,facecolor='white')
fig.savefig(out.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'},facecolor='white')
plt.close(fig)
# Editable source: labels and geometry are shared with the rendered diagram.
mx=ET.Element('mxfile',host='app.diagrams.net');di=ET.SubElement(mx,'diagram',name='Study design and participant selection')
model=ET.SubElement(di,'mxGraphModel',page='1',pageWidth='700',pageHeight='580');root=ET.SubElement(model,'root')
ET.SubElement(root,'mxCell',id='0');ET.SubElement(root,'mxCell',id='1',parent='0')
for key,x,y,w,h,t in nodes:
 c=ET.SubElement(root,'mxCell',id=key,value='',style=f'rounded=1;whiteSpace=wrap;html=0;fillColor=#ffffff;strokeColor={color};fontFamily=Arial;fontSize=12;',vertex='1',parent='1')
 ET.SubElement(c,'mxGeometry',x=str(x*100),y=str((5.8-y-h)*100),width=str(w*100),height=str(h*100),attrib={'as':'geometry'})
for i,(x1,y1,x2,y2) in enumerate(edges):
 c=ET.SubElement(root,'mxCell',id=f'edge{i}',style=f'endArrow=block;strokeColor={color};',edge='1',parent='1')
 g=ET.SubElement(c,'mxGeometry',relative='1',attrib={'as':'geometry'})
 ET.SubElement(g,'mxPoint',x=str(x1*100),y=str((5.8-y1)*100),attrib={'as':'sourcePoint'})
 ET.SubElement(g,'mxPoint',x=str(x2*100),y=str((5.8-y2)*100),attrib={'as':'targetPoint'})
for j,(x,y,t,size,bold,ha) in enumerate(labels):
    width=min(680,max(100,max(len(v) for v in t.split('\n'))*size*0.65))
    xx=x*100 if ha=='left' else x*100-width/2
    c=ET.SubElement(root,'mxCell',id=f'label{j}',value=t,style=f'text;html=0;whiteSpace=wrap;align={ha};fontFamily=Arial;fontSize={size*1.39};fontStyle={1 if bold else 0};',vertex='1',parent='1')
    ET.SubElement(c,'mxGeometry',x=str(xx),y=str((5.8-y)*100-22),width=str(width),height='44',attrib={'as':'geometry'})
ET.ElementTree(mx).write(out.with_suffix('.drawio'),encoding='utf-8',xml_declaration=True)
print('Created',out)
