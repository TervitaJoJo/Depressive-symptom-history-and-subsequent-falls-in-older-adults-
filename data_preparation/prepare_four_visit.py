"""Construct the fixed-2015 CHARLS landmark sample; no association fitting."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import pyreadstat
sys.stdout.reconfigure(encoding='utf-8')
from paths import arguments, charls_sources
root, BASE = arguments('Build the four-visit CHARLS sample.')
PRIVATE = BASE
meta, rawmeta = charls_sources(root)
items = ['depresl','effortl','sleeprl','whappyl','flonel','botherl','goingl','mindtsl','fhopel','fearll']
core = ['agey','iwstat','nhmliv','iwy','iwm','cesd10','cesd10m']
cov = ['mstat','smokev','smoken','drinkl','stroke','hearte','diabe','arthre','climsa','dressa','batha','eata','beda','toilta']
cols = ['ID','ID_w1','householdID','ragender','raeducl'] + [f'r{w}{v}' for w in range(1,5) for v in core+items] + [f'r3{v}' for v in cov]
_, hm = pyreadstat.read_dta(meta['CHARLS']['path'], metadataonly=True)
cols = [v for v in cols if v in hm.column_names]
h, _ = pyreadstat.read_dta(meta['CHARLS']['path'], usecols=cols)
h = h.set_index('ID'); h.index = h.index.astype(str)
assert h.index.is_unique
for col in h.columns:
    if col not in ['ID_w1','householdID']:
        h[col] = pd.to_numeric(h[col], errors='coerce')
raw, score, self_response, date = {}, {}, {}, {}
checks, sources = [], [meta['CHARLS']['path']]
def check(name, passed, **detail):
    checks.append({'check': name, 'passed': bool(passed), **detail})
    assert passed, (name, detail)
for w, year in enumerate([2011,2013,2015,2018], 1):
    info = rawmeta[f'CHARLS_{year}']; names = [v['name'] for v in info['columns']]
    idcol = 'ID' if 'ID' in names else 'id'
    fallvar = 'da023_w4' if w == 4 else 'da023'
    z, _ = pyreadstat.read_dta(info['path'], usecols=[idcol,fallvar,'db032'])
    sources.append(info['path']); z = z.set_index(idcol); z.index = z.index.astype(str)
    check(f'raw_id_unique_{year}', z.index.is_unique)
    if w == 1:
        z = z.reindex(h.ID_w1.astype(str)); z.index = h.index
    else:
        z = z.reindex(h.index)
    z = z.apply(pd.to_numeric, errors='coerce'); raw[w] = z
    a = h[[f'r{w}{v}' for v in items]].copy()
    complete = a.isin([1,2,3,4]).all(axis=1) & h[f'r{w}cesd10m'].eq(0)
    a = a - 1
    for v in ['whappyl','fhopel']:
        a[f'r{w}{v}'] = 3 - a[f'r{w}{v}']
    score[w] = a.sum(axis=1, min_count=10).where(complete)
    check(f'cesd10_recalculation_{year}', score[w].dropna().eq(h.loc[score[w].notna(),f'r{w}cesd10']).all(), n=int(complete.sum()))
    self_response[w] = h[f'r{w}iwstat'].eq(1) & z.db032.isin([1,2,3])
    date[w] = pd.to_datetime(dict(year=h[f'r{w}iwy'], month=h[f'r{w}iwm'], day=15), errors='coerce')
flow = []
keep = pd.Series(True, index=h.index)
def flow_step(label, condition):
    global keep
    old = int(keep.sum()); keep &= condition
    flow.append({'step':label, 'remaining_n':int(keep.sum()), 'excluded_n':old-int(keep.sum())})
flow_step('all_harmonized_participants', keep)
flow_step('age_65_or_older_in_2015', h.r3agey.ge(65))
flow_step('community_dwelling_2015', h.r3nhmliv.eq(0))
flow_step('self_participation_all_three_symptom_visits', self_response[1] & self_response[2] & self_response[3])
flow_step('complete_cesd10_all_three_symptom_visits', score[1].notna() & score[2].notna() & score[3].notna())
flow_step('current_cesd10_below_10', score[3].lt(10))
flow_step('no_recent_falls_2015', raw[3].da023.eq(2))
flow_step('valid_ordered_2011_2013_2015_dates', date[1].notna() & date[2].notna() & date[3].notna() & (date[2]>date[1]) & (date[3]>date[2]))
eligible = keep.copy()
observed = self_response[4] & raw[4].da023_w4.isin([1,2]) & date[4].notna() & (date[4]>date[3])
flow_step('self_reported_followup_with_valid_date_and_fall_response', observed)
d = pd.DataFrame(index=h.index)
d['id'] = h.index
d['household'] = h.householdID.astype('string').str.strip()
check('eligible_household_id_complete', (d.loc[eligible,'household'].notna() & d.loc[eligible,'household'].ne('')).all())
d['cluster'] = 'hh_' + d.household
for w, yr in enumerate([2011,2013,2015,2018],1):
    d[f'date_{yr}'] = date[w]
    if w < 4: d[f'cesd_{yr}'] = score[w]
d['history'] = score[1].ge(10).astype(int) + score[2].ge(10).astype(int)
d['age'] = h.r3agey
d['female'] = h.ragender.map({1:0,2:1})
d['education'] = h.raeducl.map({0:1,1:1,2:2,3:3})
d['current_score'] = score[3]
d['partnered'] = h.r3mstat.map({1:1,3:1,4:0,5:0,7:0,8:0})
bin01 = lambda s: s.where(s.isin([0,1]))
ever, now = bin01(h.r3smokev), bin01(h.r3smoken)
d['smoking'] = np.nan
d.loc[ever.eq(0),'smoking'] = 0
d.loc[ever.eq(1)&now.eq(0),'smoking'] = 1
d.loc[now.eq(1),'smoking'] = 2
d['alcohol'] = bin01(h.r3drinkl)
for v in ['stroke','hearte','diabe','arthre']:
    d[v] = bin01(h[f'r3{v}'])
d['stair_difficulty'] = bin01(h.r3climsa)
adl = h[[f'r3{v}' for v in ['dressa','batha','eata','beda','toilta']]].apply(bin01)
d['adl_difficulty'] = np.nan
d.loc[adl.eq(0).all(axis=1),'adl_difficulty'] = 0
d.loc[adl.eq(1).any(axis=1),'adl_difficulty'] = 1
for v, a, b in [('history_gap_months',1,2),('landmark_gap_months',2,3),('followup_gap_months',3,4)]:
    d[v] = (date[b]-date[a]).dt.days/30.4375
d['followup_observed'] = observed.astype(int)
d['y'] = raw[4].da023_w4.map({1:1,2:0}).where(observed)
d['followup_iwstat'] = h.r4iwstat
el = d.loc[eligible].copy(); obs = el.loc[el.followup_observed.eq(1)].copy()
check('eligible_n_matches_feasibility', len(el)==1718)
check('observed_n_matches_feasibility', len(obs)==1443)
check('total_events_matches_feasibility', int(obs.y.sum())==230)
check('eligible_group_sizes_match_feasibility', el.history.value_counts().sort_index().to_dict()=={0:1171,1:423,2:124})
check('observed_group_sizes_match_feasibility', obs.history.value_counts().sort_index().to_dict()=={0:981,1:356,2:106})
check('unique_people', el.id.is_unique)
check('positive_all_observed_time_intervals', obs[['history_gap_months','landmark_gap_months','followup_gap_months']].gt(0).all().all())
check('y_missing_if_followup_unobserved', el.loc[el.followup_observed.eq(0),'y'].isna().all())
model_cols = ['y','history','age','female','education','current_score','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty','history_gap_months','landmark_gap_months','followup_gap_months']
el.to_csv(PRIVATE/'CHARLS_four_visit_eligible.csv', index=False)
obs.to_csv(PRIVATE/'CHARLS_four_visit_observed.csv', index=False)
pd.DataFrame(flow).to_csv(BASE/'primary_cohort_flow.csv', index=False)
missing = [{'variable':v,'missing_n':int(obs[v].isna().sum()),'missing_percent':float(obs[v].isna().mean()*100),'distinct_observed_values':int(obs[v].nunique())} for v in model_cols]
pd.DataFrame(missing).to_csv(BASE/'covariate_missingness.csv', index=False)
capacity = []
for g in [0,1,2]:
    a = el.loc[el.history.eq(g)]
    capacity.append({'history':g,'eligible_n':len(a),'observed_n':int(a.followup_observed.sum())})
# Only exposure-stratified covariate support, no exposure-stratified outcome values.
support = []
for g, a in obs.groupby('history'):
    for v in model_cols:
        if v in ['y','history','age','history_gap_months','landmark_gap_months','followup_gap_months']: continue
        for value, n in a[v].value_counts(dropna=False).items():
            support.append({'history':int(g),'variable':v,'value':None if pd.isna(value) else float(value),'n':int(n)})
pd.DataFrame(support).to_csv(BASE/'pre_association_covariate_support.csv', index=False)
# Compare against the already-audited 2015 window of the original three-visit data.
old = pd.read_csv(PRIVATE/'CHARLS_analysis.csv',dtype={'id':str,'cluster':str,'group':str}).set_index('id')
old = old.loc[old.wave.eq(3)]
common = obs.index.intersection(old.index)
for v in ['age','female','education','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty','followup_gap_months']:
    left, right = obs.loc[common,v].to_numpy(float), old.loc[common,v].to_numpy(float)
    check(f'original_2015_mapping_agrees_{v}', np.allclose(left,right,equal_nan=True), n=len(common))
summary = {'created_at':datetime.now().astimezone().isoformat(),'new_associations_fitted':False,'eligible_n':len(el),'observed_n':len(obs),'events_combined':int(obs.y.sum()),'history_capacity':capacity,'n_clusters':int(obs.cluster.nunique()),'cluster_size_max':int(obs.groupby('cluster').size().max()),'complete_covariate_rows':int(obs[model_cols].notna().all(axis=1).sum()),'missingness':missing,'interval_quantiles':{v:{str(k):float(val) for k,val in obs[v].quantile([0,.25,.5,.75,1]).items()} for v in ['history_gap_months','landmark_gap_months','followup_gap_months']},'checks':checks,'private_files_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in PRIVATE.glob('CHARLS_four_visit_*.csv')}}
(BASE/'primary_pre_association_audit.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({k:summary[k] for k in ['eligible_n','observed_n','events_combined','n_clusters','cluster_size_max','complete_covariate_rows','missingness','interval_quantiles']},indent=2))
print('All construction checks passed:',len(checks))