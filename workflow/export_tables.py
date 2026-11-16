"""Export manuscript table source data and collect shareable aggregate results.
All computations use regenerated outputs; no participant rows enter this export.
"""
from pathlib import Path
import argparse,json,shutil
import numpy as np
import pandas as pd
from runtime import external_path

p=argparse.ArgumentParser(description=__doc__);p.add_argument('work_directory',type=Path)
a=p.parse_args();w=external_path(a.work_directory)
prep=w/'prepared';res=w/'results';agg=w/'aggregate_results';dec=w/'decomposition_aggregate';tables=w/'tables'
for d in [agg,dec,tables]:d.mkdir(parents=True,exist_ok=True)
def read(name,folder=agg):return pd.read_csv(folder/name)
def save(d,name,folder=agg):d.to_csv(folder/name,index=False)
# Explicit allowlist. Private matrices, predictions, IDs and imputation objects are never copied.
main_files=['sequential_models.csv','incremental_joint_wald.csv','followup_observation.csv','followup_ipw_results.csv','followup_ipw_diagnostics.csv','weighted_sensitivity.csv','charls_elsa_heterogeneity.csv']
for name in main_files:shutil.copy2(res/'main'/name,agg/name)
restricted=read('restricted_analyses.csv',res/'main').rename(columns={'model':'sensitivity'})
extra=read('elsa_extra_models.csv',res/'main')
threshold=extra[extra.model.eq('alternative_threshold_3')].rename(columns={'model':'sensitivity'})
robust=pd.concat([restricted,threshold],ignore_index=True);robust['endpoint']='fall'
cols=['cohort','sensitivity','endpoint','n','events','rr','ci_low','ci_high','p'];save(robust[cols],'robustness_sensitivity.csv')
secondary=extra[extra.model.isin(['recurrent_fall','injurious_fall'])].rename(columns={'model':'endpoint'}).copy();secondary['analysis']='secondary'
e=read('ELSA_analysis.csv',prep)
secondary['n_history']=[int((e.cesd_current.lt(4)&e.cesd_prior.ge(4)&e[ep].notna()).sum()) for ep in secondary.endpoint]
cols=['cohort','endpoint','analysis','n','events','n_history','rr','ci_low','ci_high','p'];save(secondary[cols],'elsa_secondary_outcomes.csv')
for old,new in [('rr.csv','charls_four_visit_rr.csv'),('joint.csv','charls_four_visit_joint.csv')]:
    data=read(old,res/'four')
    # The four-visit script also estimates complete-case checks; retain those in results/four.
    labels={'primary'}
    save(data[data.analysis.isin(labels)],new)
el=read('CHARLS_four_visit_eligible.csv',prep);obs=read('CHARLS_four_visit_observed.csv',prep);groups=[]
for history in range(3):
    e=el[el.history.eq(history)];d=obs[obs.history.eq(history)]
    groups.append(dict(history=history,eligible_n=len(e),observed_n=len(d),missing_outcome_n=len(e)-len(d),observed_fraction=len(d)/len(e),events=int(d.y.sum()),crude_risk=d.y.mean(),age_mean=d.age.mean(),female_fraction=d.female.mean(),current_score_mean=d.current_score.mean(),followup_months_median=d.followup_gap_months.median()))
save(pd.DataFrame(groups),'charls_four_visit_group_summary.csv')
for name in ['targeted_models.csv','targeted_diagnostics.csv','evalues.csv','continuous_score_models.csv','spline_tests.csv','spline_curves.csv','spline_knots.csv','mortality_ascertainment.csv','domain_missingness.csv','CHARLS_score_counts.csv','ELSA_score_counts.csv']:
    shutil.copy2(res/'targeted'/name,agg/name)
shutil.copy2(res/'prediction/prediction_increment.csv',agg/'prediction_increment.csv')
shutil.copy2(prep/'cohort_flow.csv',agg/'cohort_flow.csv')
# Historical-fall sensitivity uses the same allowlisted export path as other analyses.
history_agg=w/"historical_falls_aggregate";history_agg.mkdir(exist_ok=True)
for name in ['historical_falls_models.csv', 'model_diagnostics.csv', 'paired_point_estimate_changes.csv', 'capacity_and_missingness.csv', 'locked_result_verification.csv', 'CHARLS_imputation_convergence.csv', 'CHARLS_imputation_logged_events.csv', 'ELSA_imputation_convergence.csv', 'ELSA_imputation_logged_events.csv']:
    shutil.copy2(res/"historical_falls"/name,history_agg/name)
for name in ["sample_counts.csv","linkage_summary.csv"]:
    shutil.copy2(w/"history_inputs"/name,history_agg/name)

# Decomposition diagnostics are copied only by the exact, reviewed deliverable names.
refdec=Path(__file__).resolve().parents[1]/'elsa_threshold_decomposition/aggregate_results'
for ref in refdec.glob('*.csv'):
    if ref.name=='original_pipeline_comparison.csv':continue
    shutil.copy2(res/'decomposition'/ref.name,dec/ref.name)
cells=read('four_cell_results.csv',dec);seq=read('sequential_models.csv');comp=[]
for endpoint,row in [('A',seq[(seq.cohort=='ELSA')&(seq.model=='full_adjustment')].iloc[0]),('D',threshold.iloc[0])]:
    comp.append(dict(endpoint=endpoint,pipeline='original',rr=row.rr,ci_low=row.ci_low,ci_high=row.ci_high))
for endpoint in ['A','D']:
    row=cells[cells.cell.eq(endpoint)].iloc[0]
    comp.append(dict(endpoint=endpoint,pipeline='common',rr=row.exp_estimate,ci_low=row.exp_low,ci_high=row.exp_high))
save(pd.DataFrame(comp),'original_pipeline_comparison.csv',dec)
# Table 1 and S2 use observed, variable-specific denominators, without imputation.
specs=[('Participants','n',None),('Age, years','age','mean'),('Female sex','female',1),('Education: lower category','education',1),('Education: middle category','education',2),('Education: upper category','education',3),('Partnered','partnered',1),('Current smoker','smoking',2),('Alcohol use','alcohol',1),('Stroke','stroke',1),('Heart disease','hearte',1),('Diabetes','diabe',1),('Arthritis','arthre',1),('Stair-climbing difficulty','stair_difficulty',1),('Any ADL difficulty','adl_difficulty',1),('Current CES-D score','cesd_current','mean'),('Historical CES-D score','cesd_prior','mean'),('History-to-index interval, months','exposure_gap_months','mean'),('Index-to-follow-up interval, months','followup_gap_months','mean')]
cohorts={};missing=[]
missing_vars=['age','female','education','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty','cesd_current','exposure_gap_months','followup_gap_months','window']
for co,cut in [('CHARLS',10),('ELSA',4)]:
    raw=read(f'{co}_analysis.csv',prep);d=raw[raw.cesd_current.lt(cut)&raw.fall.notna()].copy();d['history']=d.cesd_prior.ge(cut).astype(int);cohorts[co]=d
    for variable in missing_vars:missing.append(dict(cohort=co,variable=variable,n=len(d),missing_n=int(d[variable].isna().sum()),missing_percent=100*d[variable].isna().mean()))
rows=[];numeric=[]
for label,var,kind in specs:
    row={'Characteristic':label}
    for co,d in cohorts.items():
        for history in [0,1]:
            z=d[d.history.eq(history)];key=f'{co}_history_{history}'
            if var=='n':value=f'{len(z):,}';numer=dict(n=len(z))
            elif kind=='mean':value=f'{z[var].mean():.1f} ({z[var].std():.1f})';numer=dict(n=int(z[var].notna().sum()),mean=z[var].mean(),sd=z[var].std())
            else:
                x=z[var].dropna();count=int(x.eq(kind).sum());percent=100*count/len(x)
                value=f'{count:,} ({percent:.1f})';numer=dict(n=len(x),count=count,percent=percent)
            row[key]=value;numeric.append(dict(cohort=co,history=history,characteristic=label,**numer))
    rows.append(row)
save(pd.DataFrame(rows),'Table_1.csv',tables);save(pd.DataFrame(numeric),'Table_1_numeric.csv',tables)
save(pd.DataFrame(missing),'Table_S2_missingness.csv',tables)
save(seq,'Table_2_source.csv',tables)
# One or more numeric source tables correspond to each displayed supplementary table.
map_tables={
'Table_S2_observation.csv':'followup_observation.csv',
'Table_S4_restrictions.csv':'robustness_sensitivity.csv',
'Table_S4_observation_weights.csv':'followup_ipw_results.csv',
'Table_S4_response_weights.csv':'weighted_sensitivity.csv',
'Table_S4_secondary_outcomes.csv':'elsa_secondary_outcomes.csv',
'Table_S3_joint_tests.csv':'incremental_joint_wald.csv',
'Table_S3_descriptive_comparison.csv':'charls_elsa_heterogeneity.csv',
'Table_S8_counts.csv':'charls_four_visit_group_summary.csv',
'Table_S8_models.csv':'charls_four_visit_rr.csv',
'Table_S8_joint_test.csv':'charls_four_visit_joint.csv',
'Table_S5_directed_sensitivity.csv':'targeted_models.csv',
'Table_S5_evalues.csv':'evalues.csv',
'Table_S7_continuous_score.csv':'continuous_score_models.csv',
'Table_S7_spline_tests.csv':'spline_tests.csv',
'Table_S7_prediction.csv':'prediction_increment.csv',
'Table_S6_mortality.csv':'mortality_ascertainment.csv',
'Table_S6_composite_models.csv':'targeted_models.csv'}
for dest,source in map_tables.items():shutil.copy2(agg/source,tables/dest)
for dest,source in [('Table_S1_cells.csv','four_cell_results.csv'),('Table_S1_contrasts.csv','pooled_contrasts.csv')]:shutil.copy2(dec/source,tables/dest)
# Table S4 panels B-C: interval-standardized risks and historical-fall sensitivity.
hist_models=read("historical_falls_models.csv",history_agg)
save(hist_models[hist_models.model.eq("both_intervals_18_30")],"Table_S4B_interval_absolute_risks.csv",tables)
save(hist_models[hist_models.model.isin(["same_sample_original","same_sample_plus_history_fall","no_falls_at_history_and_index"])],"Table_S4C_historical_falls.csv",tables)
save(read("capacity_and_missingness.csv",history_agg),"Table_S4C_sample_and_missingness.csv",tables)

notes={'Table_S4B_interval_absolute_risks':'Standardized risks and risk differences with both adjacent interviews restricted to 18-30 months; recall periods still differ by cohort.',
'Table_S4C_historical_falls':'Paired same-sample models and restriction to no recorded falls at history and index visits; use accompanying capacity/missingness counts.',
'Table_S5_directed_sensitivity':'Select medication exclusions and measured-domain rows; Table_S6 uses outcome-restriction/composite rows from the same source.',
'Table_S7_prediction':'Use point and paired-bootstrap intervals from the metric field; this is internally cross-validated comparison, not external prediction validation.',
'Table_S1_contrasts':'The displayed interval method is MI_Boot. joint_HC1 rows are a covariance check. Exponentiated component changes are ratios of RRs.',
'publication':'These CSVs contain numerical table sources, not a rebuilt Word submission. Counts use observed data. Generated outputs require a disclosure review before sharing.'}
(tables/'README.json').write_text(json.dumps(notes,indent=2),encoding='utf8')
print('Exported aggregate source tables and figure inputs.')
