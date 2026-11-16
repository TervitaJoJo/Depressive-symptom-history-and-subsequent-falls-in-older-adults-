"""Link historical falls to frozen CHARLS/ELSA windows; keep all records private."""
from pathlib import Path
import argparse, sys, json, hashlib
import numpy as np
import pandas as pd
import pyreadstat
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"workflow"))
from runtime import external_path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("raw_root",type=Path);p.add_argument("prepared",type=Path);p.add_argument("output",type=Path)
a=p.parse_args();out=external_path(a.output);out.mkdir(parents=True,exist_ok=True)
checks=[];counts=[];sources=[]
def read(path,cols):
    frame,md=pyreadstat.read_dta(path,usecols=cols)
    sources.append({"source":path.relative_to(a.raw_root).as_posix(),"columns":cols,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    return frame
h=read(a.raw_root/"CHARLS/Harmonized_CHARLS/H_CHARLS_D_Data.dta",["ID","ID_w1"]).set_index("ID")
h.index=h.index.astype(str)
fallmaps={}
for wave,year in [(1,2011),(2,2013)]:
    files=[f for f in (a.raw_root/"CHARLS"/str(year)).glob("*.dta") if f.name.lower()=="health_status_and_functioning.dta"]
    assert len(files)==1
    _,md=pyreadstat.read_dta(files[0],metadataonly=True);uid="ID" if "ID" in md.column_names else "id"
    r=read(files[0],[uid,"da023"]).set_index(uid);r.index=r.index.astype(str);assert r.index.is_unique
    v=pd.to_numeric(r.da023,errors="coerce").map({1:1,2:0})
    if wave==1:
        v=v.reindex(h.ID_w1.astype(str));v.index=h.index
    fallmaps[wave]=v
el=read(a.raw_root/"ELSA/Raw_data/UKDA-5050-stata/stata/stata13_se/h_elsa_g3.dta",["idauniq","r4fall1y"]+[f"r{w}fall" for w in range(5,8)]).set_index("idauniq")
assert el.index.is_unique
for co,cut in [("CHARLS",10),("ELSA",4)]:
    original=pd.read_csv(a.prepared/f"{co}_analysis.csv",dtype={"id":str,"cluster":str})
    d=original.copy();assert d.id.is_unique
    d["history_fall"]=np.nan
    for w in sorted(d.wave.unique()):
        ix=d.wave.eq(w);history_wave=int(w)-1
        if co=="CHARLS":v=d.loc[ix,"id"].map(fallmaps[history_wave])
        else:
            ids=pd.to_numeric(d.loc[ix,"id"],errors="raise")
            var="r4fall1y" if history_wave==4 else f"r{history_wave}fall"
            v=ids.map(el[var]).where(lambda x:x.isin([0,1]))
        d.loc[ix,"history_fall"]=v.to_numpy(dtype=float, na_value=np.nan)
        checks.append({"cohort":co,"index_wave":int(w),"history_wave":history_wave,"eligible_rows":int(ix.sum()),"missing_history_fall":int(v.isna().sum())})
    pd.testing.assert_frame_equal(d[original.columns],original)
    main=d[d.cesd_current.lt(cut)&d.fall.notna()]
    for label,z in [("primary",main),("history_fall_available",main[main.history_fall.notna()]),("history_fall_absent",main[main.history_fall.eq(0)]),("history_fall_present",main[main.history_fall.eq(1)]),("history_fall_missing",main[main.history_fall.isna()])]:
        for hist in [0,1]:
            q=z[z.cesd_prior.ge(cut).astype(int).eq(hist)]
            counts.append({"cohort":co,"sample":label,"hist":hist,"n":len(q),"events":int(q.fall.sum()),"clusters":q.cluster.nunique()})
    d.to_csv(out/f"{co}_history_falls.csv",index=False)
pd.DataFrame(checks).to_csv(out/"linkage_summary.csv",index=False)
pd.DataFrame(counts).to_csv(out/"sample_counts.csv",index=False)
(out/"source_provenance.json").write_text(json.dumps(sources,indent=2),encoding="utf8")
print(pd.DataFrame(counts).to_string(index=False))
