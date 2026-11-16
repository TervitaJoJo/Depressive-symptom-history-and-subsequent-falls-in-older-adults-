from pathlib import Path
from datetime import datetime
import json,re,sys,hashlib
import numpy as np,pandas as pd,pyreadstat
sys.stdout.reconfigure(encoding='utf-8')
from paths import arguments, charls_sources
root,p=arguments('Build CHARLS first eligible windows.')
meta,rm=charls_sources(root)
checks=[];flow=[];sources=[];capacity=[]
def read(f,cols=None):
 f=Path(f);sources.append(str(f));d,_=pyreadstat.read_dta(str(f),usecols=cols);return d
ndef=lambda s:pd.to_numeric(s,errors='coerce')
def bin01(s):return ndef(s).where(ndef(s).isin([0,1]))
def yes12(s):return ndef(s).map({1:1,2:0})
def yes15(s):return ndef(s).map({1:1,5:0})
def date(y,m):return pd.to_datetime(dict(year=ndef(y),month=ndef(m),day=15),errors='coerce')
def adl(a):
 z=pd.Series(np.nan,index=a.index);z[a.eq(0).all(axis=1)]=0;z[a.eq(1).any(axis=1)]=1;return z
ndefcheck=lambda name,n,bad:checks.append({'check':name,'n':int(n),'failed':int(bad)})
def smoke(ever,now):
 z=pd.Series(np.nan,index=ever.index);z[ever.eq(0)]=0;z[ever.eq(1)&now.eq(0)]=1;z[now.eq(1)]=2;return z

def save_cohort(co,waves,windows,cut):
 used=set();out=[];house=[]
 for w,d in waves.items():
  assert d.index.is_unique
  house.append(d[['household']].assign(wave=w))
 for wi,(a,b,c) in enumerate(windows,1):
  prior,index,follow=waves[a],waves[b],waves[c]
  z=index.join(prior[['score','self','date']].add_prefix('prior_'),how='left').join(follow[['fall','self','date']].add_prefix('future_'),how='left')
  valid=z.age.ge(65)&z.self.eq(1)&z.prior_self.eq(1)&z.score.notna()&z.prior_score.notna()&z.community.eq(1)&z.fall.eq(0)
  valid&=z.date.notna()&z.prior_date.notna()&(z.date>z.prior_date)
  flow.append({'cohort':co,'index_wave':b,'step':'all_eligible_before_first_window_selection','n':int(valid.sum())})
  valid&=~z.index.isin(used);q=z.loc[valid].copy();used.update(q.index)
  q['cesd_prior']=q.prior_score;q['cesd_current']=q.score;q['group']=(q.prior_score.ge(cut).astype(int).astype(str)+q.score.ge(cut).astype(int).astype(str))
  q['followup_observed']=q.future_self.eq(1)&q.future_date.notna()&(q.future_date>q.date)&q.future_fall.notna()
  q['fall']=q.future_fall.where(q.followup_observed)
  q['exposure_gap_months']=(q.date-q.prior_date).dt.days/30.4375;q['followup_gap_months']=(q.future_date-q.date).dt.days/30.4375
  q['wave']=b;q['window']=wi;q['cohort']=co
  out.append(q)
 d=pd.concat(out);assert d.index.is_unique
 # Household connected components across index visits among all eligible participants.
 parent={i:i for i in d.index}
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 def union(a,b):
  x,y=find(a),find(b)
  if x!=y:parent[y]=x
 for w in sorted({w[1] for w in windows}):
  hs=waves[w].reindex(d.index).dropna(subset=['household'])
  for _,g in hs.groupby('household'):
   ids=g.index.tolist()
   for i in ids[1:]:union(ids[0],i)
 d['cluster']=[str(find(i)) for i in d.index];d['id']=d.index.astype(str)
 cols=['id','cohort','cluster','wave','window','age','female','education','cesd_prior','cesd_current','group','fall','followup_observed','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty','exposure_gap_months','followup_gap_months']
 d[cols].to_csv(p/f'{co}_analysis.csv',index=False)
 low=d.cesd_current.lt(cut);obs=low&d.fall.notna();a=d.loc[obs]
 model=['age','female','education','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty','exposure_gap_months','followup_gap_months']
 groups=a.group.value_counts().to_dict();events=int(a.fall.sum())
 ok=min([groups.get('00',0),groups.get('10',0)])>=50 and events>=100
 row={'cohort':co,'eligible_all_current_scores':len(d),'current_low_eligible':int(low.sum()),'n_primary':len(a),'history_group_n':groups.get('10',0),'reference_group_n':groups.get('00',0),'events_both_groups_combined':events,'proceed_capacity':ok,'window_counts':a.window.value_counts().to_dict(),'groups_eligible':d.loc[low,'group'].value_counts().to_dict(),'covariate_missing':{x:int(a[x].isna().sum()) for x in model},'exposure_gap_quantiles':a.exposure_gap_months.quantile([0,.25,.5,.75,1]).to_dict(),'followup_gap_quantiles':a.followup_gap_months.quantile([0,.25,.5,.75,1]).to_dict(),'n_clusters':a.cluster.nunique()}
 capacity.append(row);print(json.dumps(row,ensure_ascii=False))
 ndefcheck(co+'_unique_id',len(d),d.index.duplicated().sum());ndefcheck(co+'_positive_observed_intervals',len(a),((a.exposure_gap_months<=0)|(a.followup_gap_months<=0)).sum())

# CHARLS: scores/covariates from official harmonized file, falls/proxy directly from raw.
if True:
 items=['depresl','effortl','sleeprl','whappyl','flonel','botherl','goingl','mindtsl','fhopel','fearll']
 cov=['agey','iwy','iwm','iwstat','nhmliv','mstat','smokev','smoken','drinkl','stroke','hearte','diabe','arthre','climsa','dressa','batha','eata','beda','toilta','cesd10','cesd10m']
 cols=['ID','ID_w1','householdID','ragender','raeducl']+[f'r{w}{s}' for w in [1,2,3,4] for s in cov+items]
 h=read(meta['CHARLS']['path'],cols).set_index('ID');waves={}
 for w,yr in enumerate([2011,2013,2015,2018],1):
  k=f'CHARLS_{yr}';names=[r['name'] for r in rm[k]['columns']];fallvar='da023_w4' if w==4 else 'da023';uid='ID' if 'ID' in names else 'id'
  r=read(rm[k]['path'],[uid,fallvar,'db032']).set_index(uid);r.index=r.index.astype(str)
  if w==1:
   r=r.reindex(h.ID_w1.astype(str));r.index=h.index
  else:r=r.reindex(h.index)
  a=h[[f'r{w}{s}' for s in items]].apply(ndef);complete=a.isin([1,2,3,4]).all(axis=1)&ndef(h[f'r{w}cesd10m']).eq(0)
  a=a-1
  for s in ['whappyl','fhopel']:a[f'r{w}{s}']=3-a[f'r{w}{s}']
  score=a.sum(axis=1,min_count=10).where(complete)
  ndefcheck(f'CHARLS_score_w{w}',complete.sum(),(complete&score.ne(ndef(h[f'r{w}cesd10']))).sum())
  d=pd.DataFrame(index=h.index);d['score']=score;d['date']=date(h[f'r{w}iwy'],h[f'r{w}iwm']);d['age']=ndef(h[f'r{w}agey']);d['female']=ndef(h.ragender).map({1:0,2:1});d['education']=ndef(h.raeducl).map({0:1,1:1,2:2,3:3})
  d['self']=(ndef(h[f'r{w}iwstat']).eq(1)&ndef(r.db032).isin([1,2,3])).astype(int);d['community']=ndef(h.get(f'r{w}nhmliv',pd.Series(np.nan,index=h.index))).eq(0).astype(int);d['fall']=yes12(r[fallvar]);d['household']=h.householdID
  d['partnered']=ndef(h[f'r{w}mstat']).map({1:1,3:1,4:0,5:0,7:0,8:0});d['smoking']=smoke(bin01(h[f'r{w}smokev']),bin01(h[f'r{w}smoken']));d['alcohol']=bin01(h[f'r{w}drinkl'])
  for v in ['stroke','hearte','diabe','arthre']:d[v]=bin01(h[f'r{w}{v}'])
  d['stair_difficulty']=bin01(h[f'r{w}climsa']);d['adl_difficulty']=adl(h[[f'r{w}{v}' for v in ['dressa','batha','eata','beda','toilta']]].apply(bin01));waves[w]=d
 save_cohort('CHARLS',waves,[(1,2,3),(2,3,4)],10)


pd.DataFrame(checks).to_csv(p/'CHARLS_measurement_checks.csv',index=False)
pd.DataFrame(flow).to_csv(p/'CHARLS_window_flow.csv',index=False)
(p/'CHARLS_capacity_summary.json').write_text(json.dumps(capacity,indent=2),encoding='utf8')
assert all(x['failed']==0 for x in checks)
assert capacity[0]['eligible_all_current_scores']==4062
assert capacity[0]['current_low_eligible']==2848
assert capacity[0]['n_primary']==2475
assert capacity[0]['events_both_groups_combined']==406
