from pathlib import Path
import sys,json,hashlib
import pandas as pd, numpy as np
from read_sav_fallback import read_sav_numeric
ROOT=Path(__file__).resolve().parents[1]
sav=ROOT/'data/raw/GSS2024.sav'
cols='ID YEAR HAPPY TRUST CONFINAN CONEDUC AGE SEX DEGREE MARITAL WTSSNRPS WTSSPS WTSSNRPS_AS VPSU VSTRAT SAMPLE BALLOT MODE'.split()
d,metadata=read_sav_numeric(sav,cols)
core=d.loc[d.SAMPLE==13].copy()
assert len(d)==3986 and len(core)==3309
assert int((d.SAMPLE==14).sum())==677
assert core.YEAR.eq(2024).all()
for var,expected in {'HAPPY':3281,'TRUST':946,'CONFINAN':2164,'CONEDUC':2164}.items():
    got=int(core[var].isin([1,2,3]).sum())
    assert got==expected,(var,got,expected)
assert core.WTSSNRPS.gt(0).all()
fields=['HAPPY','TRUST','CONFINAN','CONEDUC']
valid=pd.DataFrame({v:core[v].isin([1,2,3]) for v in fields})
counts=pd.DataFrame({'variable':fields,'valid_responses':[int(valid[v].sum()) for v in fields], 'ineligible_or_reserved':[int((~valid[v]).sum()) for v in fields]})
joint=pd.DataFrame({x:[int((valid[x]&valid[y]).sum()) for y in fields] for x in fields},index=fields)
exps=[]
for v in fields:
    mask=valid[v]
    denominator=core.loc[mask,'WTSSNRPS'].sum()
    for i in [1,2,3]:
        sub=mask&(core[v]==i)
        exps.append({'variable':v,'category':i,'valid_n':int(mask.sum()),'category_n':int(sub.sum()),'unweighted_pct':float(100*sub.sum()/mask.sum()),'weighted_pct_preview':float(100*core.loc[sub,'WTSSNRPS'].sum()/denominator)})
table=pd.DataFrame(exps)
missing_codes=[]
for col in fields + ['AGE','SEX','DEGREE','MARITAL']:
    valid_values=([1,2,3] if col in fields else None)
    x=core[col]
    if valid_values is not None: excluded=x[~x.isin(valid_values)]
    else: excluded=x[x.isna() | (x<0)]
    for code,nn in excluded.fillna('SYSTEM_MISSING').value_counts(dropna=False).items():
        missing_codes.append({'variable':col,'reserved_code':code,'count':int(nn)})
missing_codes=pd.DataFrame(missing_codes)
missing_codes.to_csv(ROOT/'qa/generated/reserved_code_frequencies.csv',index=False)
res=pd.DataFrame({'source':'2024 GSS Release 3a','original_file_bytes':sav.stat().st_size,'SHA256':hashlib.sha256(sav.read_bytes()).hexdigest(),'total_rows':[len(d)],'core_rows':[len(core)],'amerispeak_rows':[int((d.SAMPLE==14).sum())],'main_weight':'WTSSNRPS'})
for nm,x in [('eligibility_counts.csv',counts),('joint_eligibility_matrix.csv',joint),('unweighted_and_weighted_previews.csv',table),('source_manifest.csv',res)]:
 path=ROOT/'qa/generated'/nm;path.parent.mkdir(parents=True,exist_ok=True);x.to_csv(path,index=True if nm=='joint_eligibility_matrix.csv' else False)
print('Rows:',len(d),'Core:',len(core),'Amerispeak:',len(d)-len(core))
print('Eligible\n',counts.to_string(index=False))
print('Overlaps\n',joint.to_string())
print('Distributions:\n',table.to_string(index=False,float_format=lambda a:f'{a:.1f}'))
print('Validation: 4 counts match official codebook')
