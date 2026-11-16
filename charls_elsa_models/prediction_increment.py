"""Exploratory, paired, internally cross-validated prediction comparison.
No outcome-dependent imputation or cross-household leakage. Bootstrap resamples
households within fixed folds and refits both models and preprocessing each time.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path

import os
os.environ.setdefault('OMP_NUM_THREADS','1');os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('MKL_NUM_THREADS','1')
import sys,json,warnings,concurrent.futures,time
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score,log_loss
from threadpoolctl import threadpool_limits
from sklearn.exceptions import ConvergenceWarning
sys.stdout.reconfigure(encoding='utf8')
CAT=['female','education','window','cesd_current','partnered','smoking','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty']
NUM=['age','exposure_gap_months']

def weighted_median(x,w):
 ok=np.isfinite(x)&(w>0);a=x[ok];ww=w[ok];o=np.argsort(a);return a[o][np.searchsorted(np.cumsum(ww[o]),ww.sum()/2)]
def rcs(x):
 k=[65.,70.,80.,90.];cols=[x]
 for a in k[:-2]:cols.append((np.maximum(x-a,0)**3-np.maximum(x-k[-2],0)**3*(k[-1]-a)/(k[-1]-k[-2])+np.maximum(x-k[-1],0)**3*(k[-2]-a)/(k[-1]-k[-2]))/(k[-1]-k[0])**2)
 return np.column_stack(cols)
def design(d,train,weights,co):
 parts=[]
 for c in CAT:
  x=d[c].to_numpy(float);valid=train&np.isfinite(x)&(weights>0)
  levels=([1,2,3] if c=='education' else ([1,2] if co=='CHARLS' else [0,1,2,3]) if c=='window' else list(range(10 if co=='CHARLS' else 4)) if c=='cesd_current' else [0,1,2] if c=='smoking' else [0,1])
  counts=[weights[valid&(x==v)].sum() for v in levels];mode=levels[int(np.argmax(counts))];missing=~np.isfinite(x);a=np.where(missing,mode,x)
  parts.extend([(a==v).astype(float) for v in levels[1:]]);parts.append(missing.astype(float))
 num=[]
 for c in NUM:
  x=d[c].to_numpy(float);med=weighted_median(x[train],weights[train]);missing=~np.isfinite(x);x=np.where(missing,med,x)
  if c=='age':num.extend(rcs(x).T)
  else:num.append(x)
  parts.append(missing.astype(float))
 num=np.column_stack(num);mean=np.average(num[train],axis=0,weights=weights[train]);sd=np.sqrt(np.average((num[train]-mean)**2,axis=0,weights=weights[train]));sd[sd<1e-12]=1
 return np.column_stack([*parts,*((num-mean)/sd).T])

def metrics(y,p0,p1,w):
 a0=roc_auc_score(y,p0,sample_weight=w);a1=roc_auc_score(y,p1,sample_weight=w)
 b0=np.average((y-p0)**2,weights=w);b1=np.average((y-p1)**2,weights=w)
 l0=log_loss(y,p0,sample_weight=w,labels=[0,1]);l1=log_loss(y,p1,sample_weight=w,labels=[0,1])
 return dict(auc_base=a0,auc_history=a1,delta_auc=a1-a0,brier_base=b0,brier_history=b1,delta_brier=b1-b0,logloss_base=l0,logloss_history=l1,delta_logloss=l1-l0)

def pipeline(d,y,hist,fold,w,co):
 pp=np.full((len(d),2),np.nan);maxiter=0
 for f in range(5):
  tr=(fold!=f)&(w>0);te=fold==f;x=design(d,tr,w,co)
  for j in range(2):
   xx=np.column_stack([x,hist]) if j else x
   model=LogisticRegression(C=1.0,solver='lbfgs',max_iter=1000,tol=1e-7)
   with warnings.catch_warnings():
    warnings.simplefilter('error',ConvergenceWarning);model.fit(xx[tr],y[tr],sample_weight=w[tr])
   pp[te,j]=model.predict_proba(xx[te])[:,1];maxiter=max(maxiter,int(model.n_iter_[0]))
 assert np.isfinite(pp).all()
 return metrics(y,pp[:,0],pp[:,1],w),pp,maxiter

def main():
 if len(sys.argv)<3:raise SystemExit('Usage: python prediction_increment.py INPUT_DIRECTORY OUTPUT_DIRECTORY [BOOTSTRAPS=500]')
 inp=external_path(sys.argv[1]);out=external_path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True);nb=int(sys.argv[3]) if len(sys.argv)>3 else 500
 rows=[];logs=[]
 with threadpool_limits(limits=1):
  for co in ['CHARLS','ELSA']:
   d=pd.read_csv(inp/f'{co}_extended.csv',dtype={'id':str,'cluster':str});cut=10 if co=='CHARLS' else 4;d=d[d.cesd_current.lt(cut)&d.fall.notna()].reset_index(drop=True)
   y=d.fall.to_numpy(int);hist=d.cesd_prior.ge(cut).to_numpy(float);g=d.cluster.to_numpy();fold=np.full(len(d),-1,int)
   for f,(_,te) in enumerate(StratifiedGroupKFold(5,shuffle=True,random_state=20261031).split(d,y,g)):fold[te]=f
   assert pd.DataFrame({'g':g,'fold':fold}).groupby('g').fold.nunique().max()==1
   point,pp,mi=pipeline(d,y,hist,fold,np.ones(len(d)),co)
   pd.DataFrame({'id':d.id,'cluster':g,'fold':fold,'y':y,'base_probability':pp[:,0],'history_probability':pp[:,1]}).to_csv(inp/f'{co}_prediction_private.csv',index=False)
   groups,gi=np.unique(g,return_inverse=True);gf=np.array([fold[np.where(gi==i)[0][0]] for i in range(len(groups))]);foldgroups=[np.where(gf==f)[0] for f in range(5)]
   def boot(b):
    rng=np.random.default_rng(20261031+10000*(co=='ELSA')+b);gw=np.zeros(len(groups))
    for ids in foldgroups:gw[ids]=rng.multinomial(len(ids),np.repeat(1/len(ids),len(ids)))
    vals,_,iters=pipeline(d,y,hist,fold,gw[gi],co);return {'replicate':b,'max_iterations':iters,**vals}
   allboot=[];start=time.time();print(co,'point',point,flush=True)
   with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for vals in pool.map(boot,range(nb)):
     allboot.append(vals)
     if len(allboot)%50==0:print(co,'bootstrap',len(allboot),'/',nb,'seconds',round(time.time()-start),flush=True)
   bs=pd.DataFrame(allboot);bs.to_csv(out/f'{co}_prediction_bootstrap.csv',index=False)
   for metric,val in point.items():
    lo,hi=np.quantile(bs[metric],[.025,.975]);rows.append(dict(cohort=co,n=len(d),events=int(y.sum()),households=len(groups),metric=metric,estimate=val,ci_low=lo,ci_high=hi,bootstrap_repetitions=nb))
   logs.append(dict(cohort=co,n=len(d),events=int(y.sum()),fold_sizes=np.bincount(fold).tolist(),fold_events=[int(y[fold==f].sum()) for f in range(5)],households=len(groups),no_household_leakage=True,max_iterations=max(mi,int(bs.max_iterations.max())),bootstrap_failures=0,seed=20261031))
   pd.DataFrame(rows).to_csv(out/'prediction_increment.csv',index=False)
 (out/'prediction_audit.json').write_text(json.dumps({'models':'Fixed L2-penalized logistic regression, C=1; age restricted cubic spline knots 65/70/80/90; training-only weighted mode/median imputation and numerical scaling; all covariate missingness indicators; no future-observed followup duration; no tuning','validation':'Household-grouped five-fold CV; bootstrap households within fixed folds, refitting preprocessing and both models per replicate; percentile intervals conditional on split; no external transport validation','cohorts':logs},indent=2),encoding='utf8')
 print('DONE',flush=True)
if __name__=='__main__':main()
