"""ELSA primary-window construction; restricted inputs/outputs remain local."""
from pathlib import Path
from datetime import datetime
import json, hashlib, sys, re
import numpy as np, pandas as pd, pyreadstat
sys.stdout.reconfigure(encoding='utf-8')
from paths import arguments
root, p = arguments('Build the first eligible ELSA window from licensed harmonized data.')
items = ['depres', 'effort', 'sleepr', 'whappy', 'flone', 'enlife', 'fsad', 'going']
res = []
flows = []
groups = []
miss = []
checks = []
rawqc = []
for co in ['ELSA']:
    h = False
    uid = 'idauniq'
    age = 'agey'
    idxs = [5, 6, 7, 8]
    f = root / 'ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta'
    _, md = pyreadstat.read_dta(str(f), metadataonly=True)
    names = set(md.column_names)
    cov = ['smoken', 'smokev', 'drink', 'mstat', 'hibpe', 'diabe', 'hearte', 'stroke', 'arthre', 'shlt', 'clim1a', 'dressa', 'batha', 'eata', 'beda', 'toilta']
    need = [uid, 'ragender', 'raeducl']
    for w in range(min(idxs) - 1, max(idxs) + 2):
        suff = ['cesd', 'cesdm', 'iwstat', 'proxy', age] + items + cov + ['iwindm', 'iwindy', 'fall', 'fallnum', 'fallinj', 'walk100', 'bmi']
        need += [f'r{w}{x}' for x in suff if f'r{w}{x}' in names]
    d, m = pyreadstat.read_dta(str(f), usecols=list(dict.fromkeys(need)))
    d[uid] = pd.to_numeric(d[uid])
    assert d[uid].is_unique
    used = set()
    assembled = []
    for iw, w in enumerate(idxs):
        w0 = w - 1
        w2 = w + 1
        stat = [1]
        valid = d[f'r{w}{age}'].ge(65)
        flows.append({'cohort': co, 'wave': w, 'step': 'age_ge65', 'n': int(valid.sum())})
        for v in [w0, w]:
            valid &= d[f'r{v}iwstat'].isin(stat) & d[f'r{v}proxy'].eq(0)
        flows.append({'cohort': co, 'wave': w, 'step': 'two_self_interviews', 'n': int(valid.sum())})
        for v in [w0, w]:
            cols = [f'r{v}{s}' for s in items]
            assert all((c in d for c in cols)), (co, v, [c for c in cols if c not in d])
            complete = d[cols].isin([0, 1]).all(axis=1) & d[f'r{v}cesdm'].eq(0) & d[f'r{v}cesd'].between(0, 8)
            recomputed = d[cols].copy()
            for s in ['whappy', 'enlife']:
                recomputed[f'r{v}{s}'] = 1 - recomputed[f'r{v}{s}']
            agrees = recomputed.sum(axis=1).eq(d[f'r{v}cesd'])
            checks.append({'cohort': co, 'check': f'cesd_recalculation_w{v}_window{w}', 'n_checked': int(complete.sum()), 'n_failed': int((complete & ~agrees).sum())})
            valid &= complete
        flows.append({'cohort': co, 'wave': w, 'step': 'two_complete_cesd', 'n': int(valid.sum())})
        valid &= d[f'r{w}fall'].eq(0)
        flows.append({'cohort': co, 'wave': w, 'step': 'current_no_falls', 'n': int(valid.sum())})
        valid &= ~d[uid].isin(used)
        q = d.loc[valid].copy()
        used.update(q[uid].tolist())
        flows.append({'cohort': co, 'wave': w, 'step': 'first_eligible_unique', 'n': len(q)})
        y = pd.to_numeric(q[f'r{w2}fall'], errors='coerce').where(lambda s: s.isin([0, 1]))
        observed = q[f'r{w2}iwstat'].isin(stat) & q[f'r{w2}proxy'].eq(0)
        y = y.where(observed)
        num = pd.to_numeric(q[f'r{w2}fallnum'], errors='coerce')
        recur = pd.Series(np.nan, index=q.index)
        recur.loc[y.eq(0)] = 0
        recur.loc[y.eq(1) & num.ge(1)] = num.ge(2).astype(float)
        inj = pd.to_numeric(q[f'r{w2}fallinj'], errors='coerce').where(lambda s: s.isin([0, 1])).where(observed)
        inj.loc[y.eq(0)] = 0
        z = pd.DataFrame({'id': q[uid], 'cohort': co, 'window': iw, 'wave': w, 'age': q[f'r{w}{age}'], 'female': q.ragender.map({1: 0, 2: 1}), 'education_raw': q['raeducl'], 'cesd_prior': q[f'r{w0}cesd'], 'cesd_current': q[f'r{w}cesd'], 'fall': y, 'recurrent_fall': recur, 'injurious_fall': inj, 'followup_status': q[f'r{w2}iwstat']})
        z['group'] = z.cesd_prior.ge(4).astype(int).astype(str) + z.cesd_current.ge(4).astype(int).astype(str)
        for suf in cov:
            for v, tag in [(w0, 'prior'), (w, 'current')]:
                z[tag + '_' + suf] = q[f'r{v}{suf}'] if f'r{v}{suf}' in q else np.nan
        for v, tag in [(w0, 'prior'), (w, 'index'), (w2, 'followup')]:
            z[tag + '_date'] = pd.to_datetime(dict(year=q[f'r{v}iwindy'], month=q[f'r{v}iwindm'], day=15), errors='coerce')
        z['exposure_gap_months'] = (z.index_date - z.prior_date).dt.days / 30.4375
        z['followup_gap_months'] = (z.followup_date - z.index_date).dt.days / 30.4375
        assembled.append(z)
    out = pd.concat(assembled, ignore_index=True)
    assert out.id.is_unique
    out.to_csv(p / f'{co}_capacity.csv', index=False)
    for group, g in out.groupby('group'):
        groups.append({'cohort': co, 'group': group, 'n_eligible': len(g), 'n_outcome_known': int(g.fall.notna().sum())})
    low = out.cesd_current.lt(4)
    row = {'cohort': co, 'eligible': len(out), 'outcome_known': int(out.fall.notna().sum()), 'any_fall_events_overall': int(out.fall.sum()), 'unknown_outcome': int(out.fall.isna().sum()), 'recurrent_events_overall': int(out.recurrent_fall.sum()), 'recurrent_known': int(out.recurrent_fall.notna().sum()), 'injurious_events_overall': int(out.injurious_fall.sum()), 'injurious_known': int(out.injurious_fall.notna().sum()), 'current_low_n': int(low.sum()), 'current_low_outcome_known': int(out.loc[low, 'fall'].notna().sum()), 'current_low_events_combined_00_10': int(out.loc[low, 'fall'].sum()), 'exposure_gap_quantiles': out.exposure_gap_months.quantile([0, 0.25, 0.5, 0.75, 1]).to_dict(), 'followup_gap_quantiles_known': out.loc[out.fall.notna(), 'followup_gap_months'].quantile([0, 0.25, 0.5, 0.75, 1]).to_dict()}
    res.append(row)
    for c in out.columns:
        if c.startswith(('prior_', 'current_')) or c in ['age', 'female', 'education_raw']:
            miss.append({'cohort': co, 'variable': c, 'n': len(out), 'missing': int(out[c].isna().sum())})
    checks.append({'cohort': co, 'check': 'selected_person_id_unique', 'n_checked': len(out), 'n_failed': int(out.id.duplicated().sum())})
    print(co, json.dumps(row, ensure_ascii=False))
    print('EXPOSURE GROUP CAPACITY', [g for g in groups if g['cohort'] == co])
summary = {'created_at': datetime.now().astimezone().isoformat(), 'association_models_run': False, 'cohorts': res}
(p / 'ELSA_capacity_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
(p / 'ELSA_fall_coding_qc.json').write_text(json.dumps(rawqc, ensure_ascii=False, indent=2), encoding='utf-8')
for name, rows in [('flow', flows), ('exposure_group_capacity', groups), ('covariate_missingness', miss), ('capacity_checks', checks)]:
    pd.DataFrame(rows).to_csv(p / ('ELSA_' + name + '.csv'), index=False)
print('CHECKS', len(checks), 'FAILED', sum((c['n_failed'] for c in checks)))
from pathlib import Path
import pandas as pd, numpy as np, pyreadstat, json, sys
sys.stdout.reconfigure(encoding='utf-8')
log = []
for co in ['ELSA']:
    h = False
    uid = 'idauniq'
    waves = [5, 6, 7, 8]
    f = root / 'ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta'
    cols = [uid] + [c for w in waves for c in [f'r{w}nhmliv', f'hh{w}hhid', f'r{w}cwtresp']]
    z, _ = pyreadstat.read_dta(str(f), usecols=cols)
    z[uid] = pd.to_numeric(z[uid])
    z = z.rename(columns={uid: 'id'})
    z = z.apply(pd.to_numeric, errors='coerce')
    d = pd.read_csv(p / f'{co}_capacity.csv', dtype={'group': str}).merge(z, on='id', validate='1:1')
    d['institution'] = np.nan
    d['survey_weight'] = np.nan
    for w in waves:
        b = d.wave.eq(w)
        d.loc[b, 'institution'] = d.loc[b, f'r{w}nhmliv']
        d.loc[b, 'survey_weight'] = d.loc[b, f'r{w}cwtresp']
    oldn = len(d)
    d = d[d.institution.eq(0)].copy()
    parent = {i: i for i in d.id}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        aa, bb = (find(a), find(b))
        if aa != bb:
            parent[max(aa, bb)] = min(aa, bb)
    for w in waves:
        c = f'hh{w}hhid'
        for _, g in d.dropna(subset=[c]).groupby(c):
            ids = g.id.tolist()
            for j in ids[1:]:
                union(ids[0], j)
    d['cluster'] = [find(i) for i in d.id]
    d['education'] = d.education_raw.where(d.education_raw.isin([1, 2, 3]))
    d['partnered'] = d.current_mstat.map({1: 1, 2: 1, 3: 1, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0})
    d['smoking'] = np.nan
    d.loc[d.current_smoken.eq(0) & d.current_smokev.eq(0), 'smoking'] = 0
    d.loc[d.current_smoken.eq(0) & d.current_smokev.eq(1), 'smoking'] = 1
    d.loc[d.current_smoken.eq(1), 'smoking'] = 2
    d['alcohol'] = d.current_drink.where(d.current_drink.isin([0, 1]))
    for s in ['stroke', 'hearte', 'diabe', 'arthre']:
        d[s] = d['current_' + s].where(d['current_' + s].isin([0, 1]))
    d['stair_difficulty'] = d.current_clim1a.map({0: 0, 1: 1})
    ac = d[['current_' + s for s in ['dressa', 'batha', 'eata', 'beda', 'toilta']]]
    pos = [1]
    d['adl_difficulty'] = np.nan
    d.loc[ac.eq(0).all(axis=1), 'adl_difficulty'] = 0
    d.loc[ac.isin(pos).any(axis=1), 'adl_difficulty'] = 1
    cols = ['id', 'cohort', 'cluster', 'wave', 'window', 'age', 'female', 'education', 'cesd_prior', 'cesd_current', 'group', 'fall', 'recurrent_fall', 'injurious_fall', 'followup_status', 'partnered', 'smoking', 'alcohol', 'stroke', 'hearte', 'diabe', 'arthre', 'stair_difficulty', 'adl_difficulty', 'exposure_gap_months', 'followup_gap_months', 'survey_weight']
    d[cols].to_csv(p / f'{co}_analysis.csv', index=False)
    low = d.cesd_current.lt(4)
    obs = d.fall.notna()
    row = {'cohort': co, 'before_common_community_filter': oldn, 'excluded_institution_or_unknown': oldn - len(d), 'eligible_final': len(d), 'low_current_eligible': int(low.sum()), 'low_current_observed': int((low & obs).sum()), 'main_group_counts_observed': d.loc[low & obs, 'group'].value_counts().to_dict(), 'household_clusters': int(d.cluster.nunique()), 'max_cluster_size': int(d.groupby('cluster').size().max()), 'model_covariate_missingness': {c: int(d.loc[low & obs, c].isna().sum()) for c in cols if c not in ['fall', 'recurrent_fall', 'injurious_fall']}}
    log.append(row)
    print(json.dumps(row, ensure_ascii=False))
(p / 'ELSA_preparation_audit.json').write_text(json.dumps(log, ensure_ascii=False, indent=2), encoding='utf-8')

assert len(d) == 6213
assert int((low & obs).sum()) == 4779
assert int(d.loc[low & obs, "fall"].sum()) == 1134
assert d.loc[low & obs, ["exposure_gap_months", "followup_gap_months"]].gt(0).all().all()
assert all(x["n_failed"] == 0 for x in checks)
