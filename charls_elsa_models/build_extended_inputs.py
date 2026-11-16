"""Enrich authorized prepared CHARLS/ELSA extracts with targeted sensitivity variables.
Run data_preparation/build_inputs.py first to construct the original eligibility windows. All generated rows
are individual-level data: keep INPUT and OUTPUT in authorized private storage.
"""
from pathlib import Path
import json,re,sys,pandas as pd,numpy as np,pyreadstat
sys.stdout.reconfigure(encoding='utf8')
if len(sys.argv)!=4:
 raise SystemExit('Usage: python build_extended_inputs.py RAW_ROOT BASE_INPUT_DIRECTORY PRIVATE_OUTPUT_DIRECTORY')
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path
r=Path(sys.argv[1]);base=Path(sys.argv[2]);p=external_path(sys.argv[3]);private=p
p.mkdir(parents=True,exist_ok=True)
hpaths={'CHARLS':r/'CHARLS/Harmonized_CHARLS/H_CHARLS_D_Data.dta',
        'ELSA':r/'ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta'}
logs=[]
for co in ['CHARLS','ELSA']:
 raw=base/f'{co}_analysis.csv'
 d=pd.read_csv(raw,dtype={'id':str,'cluster':str});d['id']=d.id.str.replace(r'\.0$','',regex=True)
 f=str(hpaths[co]);_,md=pyreadstat.read_dta(f,metadataonly=True); uid='ID' if co=='CHARLS' else 'idauniq';waves=[1,2,3,4] if co=='CHARLS' else list(range(4,10))
 cols=[uid]+(['ID_w1'] if co=='CHARLS' else [])
 suffix=['iwstat','imrc','dlrc','wspeed','balance','iwy','iwm','iwindy','iwindm','psyche','dsight','nsight','rxpsych']
 cols += [f'r{w}{s}' for w in waves for s in suffix if f'r{w}{s}' in md.column_names]
 h,_=pyreadstat.read_dta(f,usecols=cols);h[uid]=h[uid].astype(str).str.replace(r'\.0$','',regex=True);h=h.set_index(uid);assert h.index.is_unique; h.loc[:, [c for c in h if c != "ID_w1"]] = h.loc[:, [c for c in h if c != "ID_w1"]].apply(pd.to_numeric,errors="coerce"); h = h.astype({c:float for c in h if c != "ID_w1"})
 for c in ['recorded_death','target_status','recall_total','poor_eyesight','log_walk_time','balance_score','antidepressant_report','sedative_report','medication_item_completed','elsa_prior_med_observed','elsa_prior_antidepressant'] : d[c]=np.nan
 for w in sorted(d.wave.unique()):
  take=d.wave.eq(w);ii=d.index[take];q=h.reindex(d.loc[ii,'id']).reset_index(drop=True);q.index=ii
  status=q[f'r{w+1}iwstat'];d.loc[ii,'target_status']=status;d.loc[ii,'recorded_death']=status.isin([5,6]).astype(int)
  a=q[f'r{w}imrc'].where(q[f'r{w}imrc'].between(0,10));c=q[f'r{w}dlrc'].where(q[f'r{w}dlrc'].between(0,10));d.loc[ii,'recall_total']=a+c
  walk=q[f'r{w}wspeed'];d.loc[ii,'log_walk_time']=np.log(walk.where(walk.gt(0)))
  if co=='CHARLS':
   d.loc[ii,'balance_score']=q[f'r{w}balance'].where(q[f'r{w}balance'].isin([1,2,3,4]))
  else:
   a=q[f'r{w}dsight'].where(q[f'r{w}dsight'].isin([1,2,3,4,5,6]));c=q[f'r{w}nsight'].where(q[f'r{w}nsight'].isin([1,2,3,4,5,6]));v=pd.Series(np.nan,index=ii);v[a.notna()&c.notna()]=0;v[a.ge(4)|c.ge(4)]=1;d.loc[ii,'poor_eyesight']=v
 if co=='CHARLS':
  rawwaves={}
  for w,year in [(1,2011),(2,2013),(3,2015)]:
   ff=str(r/f'CHARLS/{year}/Health_Status_and_Functioning.dta');_,mm=pyreadstat.read_dta(ff,metadataonly=True);ruid='ID' if 'ID' in mm.column_names else 'id';cc=[ruid,'da033','da034']+[n for n in mm.column_names if re.fullmatch('da020s[1-4]',n)]
   z,_=pyreadstat.read_dta(ff,usecols=cc);z[ruid]=z[ruid].astype(str).str.replace(r'\.0$','',regex=True);z=z.set_index(ruid)
   if w==1:z=z.reindex(h.ID_w1.astype(str));z.index=h.index
   rawwaves[w]=z
  for w in [2,3]:
   ii=d.index[d.wave.eq(w)];prior=rawwaves[w-1].reindex(d.loc[ii,'id']);current=rawwaves[w].reindex(d.loc[ii,'id']);prior.index=current.index=ii
   d.loc[ii,'antidepressant_report']=(prior.da020s2.eq(2)|current.da020s2.eq(2)).astype(int)
   d.loc[ii,'sedative_report']=(prior.da020s3.eq(3)|current.da020s3.eq(3)).astype(int)
   d.loc[ii,'medication_item_completed']=(prior.filter(regex='da020s').notna().any(axis=1)&current.filter(regex='da020s').notna().any(axis=1)).astype(int)
   a=current.da033.where(current.da033.isin([1,2,3,4,5]));c=current.da034.where(current.da034.isin([1,2,3,4,5]));v=pd.Series(np.nan,index=ii);v[a.notna()&c.notna()]=0;v[a.ge(4)|c.ge(4)]=1;d.loc[ii,'poor_eyesight']=v
 else:
  nf=r/'ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/wave6/wave_6_elsa_nurse_data_v2.dta';n,nm=pyreadstat.read_dta(str(nf));n['idauniq']=n.idauniq.astype(str).str.replace(r'\.0$','',regex=True);n=n.set_index('idauniq')
  ii=d.index[d.wave.eq(7)];z=n.reindex(d.loc[ii,'id']);z.index=ii;q=h.reindex(d.loc[ii,'id']);q.index=ii
  before=(z.visyear*12+z.vismon)<(q.r7iwindy*12+q.r7iwindm)
  drugcols=[c for c in z if re.fullmatch('DrC[0-9]+',c)]
  yes=pd.Series(False,index=ii);valid=pd.Series(False,index=ii)
  for c in drugcols:
   j=int(c[3:]);taken='MedBIA' if j==1 else f'MedBIA{j}';codes=z[c].astype('string').str.strip()
   valid |= codes.str.fullmatch('[0-9]{6}',na=False)&z[taken].eq(1)
   yes |= codes.str.startswith('0403',na=False)&z[taken].eq(1)
  observed=before&(z.MedCNJD.eq(2)|(z.MedCNJD.eq(1)&valid))
  d.loc[ii,'elsa_prior_med_observed']=observed.astype(int);d.loc[ii,'elsa_prior_antidepressant']=yes.astype(float).where(observed)
  logs.append({'cohort':co,'drug_definition':'Wave 6 nurse: BNF 0403*, taken in previous 7 days; nurse month strictly precedes wave 7 index month; at least one classifiable current drug or explicit no medicine. Snapshot covers only prior wave, not all historical/current use.','codes':sorted({x for c in drugcols for x in z[c].dropna().astype(str).unique() if x.startswith('0403')})})
 low=d.cesd_current.lt(10 if co=='CHARLS' else 4);main=low&d.fall.notna()
 assert not ((d.recorded_death==1)&d.fall.notna()).any(), 'Death and observed falls overlap: investigate timing'
 d.to_csv(private/f'{co}_extended.csv',index=False)
 info={'cohort':co,'eligible_low':int(low.sum()),'main_n':int(main.sum()),'recorded_deaths':int(d.loc[low,'recorded_death'].sum()),'target_status':d.loc[low,'target_status'].value_counts(dropna=False).to_dict(),'antidepressant_reports':int(d.loc[main,'antidepressant_report'].sum()),'sedative_reports':int(d.loc[main,'sedative_report'].sum()),'domain_missing':d.loc[main,['recall_total','poor_eyesight','log_walk_time','balance_score']].isna().sum().to_dict(),'elsa_prior_med_observed':int(d.loc[main,'elsa_prior_med_observed'].sum()),'elsa_prior_antidepressants':int(d.loc[main,'elsa_prior_antidepressant'].sum())}
 logs.append(info);print(json.dumps(info,ensure_ascii=False))
(p/'data_extension_audit.json').write_text(json.dumps(logs,ensure_ascii=False,indent=2),encoding='utf8')
