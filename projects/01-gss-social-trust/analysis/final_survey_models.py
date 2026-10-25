"""Reproducible GSS 2024 core survey analysis: weighted proportional-odds, diagnostics and comparisons.

Prerequisite: user-downloaded official NORC Release 3a GSS2024.sav in data/raw/.
The native SPSS file is not distributed in this repository.
Use a narrow-purpose SPSS reader when pyreadstat is unavailable (verified against GSS codebook
and independent likelihood replication with statsmodels on unweighted data).
Survey variance: PSU-centred stratified with-replacement linearization (67 df),
using weighted observation scores; a sensitivity assumption, not an exact replication of R survey.
"""
from pathlib import Path
import hashlib, json, warnings
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.optimize import minimize
from scipy import stats
from statsmodels.miscmodels.ordinal_model import OrderedModel
import statsmodels.api as sm
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from read_sav_fallback import read_sav_numeric

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/'deliverables'/'generated'; DIR.mkdir(parents=True,exist_ok=True)
Q=ROOT/'qa'/'generated';Q.mkdir(parents=True,exist_ok=True)
SAV=ROOT/'data/raw/GSS2024.sav'
selected='ID SAMPLE HAPPY TRUST CONFINAN CONEDUC AGE SEX DEGREE MARITAL WTSSNRPS WTSSPS VSTRAT VPSU MODE YEAR'.split()
try:
    import pyreadstat
    raw,_=pyreadstat.read_sav(str(SAV),usecols=selected,apply_value_formats=False,user_missing=True)
    raw.columns=raw.columns.str.upper()
    parser='pyreadstat'
except ImportError:
    raw,meta=read_sav_numeric(SAV, selected)
    parser='narrow-purpose bytecode SAV reader, see documented limitations'
assert len(raw)==3986 and (raw['SAMPLE']==13).sum()==3309
assert (raw['YEAR']==2024).all()
d=raw[raw['SAMPLE']==13].copy()
for col,allowed in {'HAPPY':[1,2,3],'TRUST':[1,2,3],'CONFINAN':[1,2,3],'CONEDUC':[1,2,3],
                    'SEX':[1,2],'AGE':list(range(18,90)),'DEGREE':[0,1,2,3,4],'MARITAL':[1,2,3,4,5]}.items():
    d[col]=d[col].where(d[col].isin(allowed),np.nan)
assert [d[x].notna().sum() for x in ['HAPPY','TRUST','CONFINAN','CONEDUC']]==[3281,946,2164,2164]
assert ((d['WTSSNRPS']>0)&d['VSTRAT'].notna()&d['VPSU'].notna()).all()
strata=d[['VSTRAT','VPSU']].drop_duplicates().groupby('VSTRAT').size()
assert (strata==2).all() and len(strata)==67
DF=int(strata.sum()-len(strata))

model=d.dropna(subset=['HAPPY','TRUST','AGE','SEX','DEGREE','MARITAL','WTSSNRPS','VSTRAT','VPSU']).copy()
assert len(model)==911
model['female']=(model.SEX==2).astype(float)
model['trust_careful']=(model.TRUST==2).astype(float)
model['trust_depends']=(model.TRUST==3).astype(float)
for j in [1,2,3,4]:model[f'degree_{j}']=(model.DEGREE==j).astype(float)
for j in [2,3,4,5]:model[f'marital_{j}']=(model.MARITAL==j).astype(float)
for j in [2,3,4]:model[f'mode_{j}']=(model.MODE==j).astype(float)
cols=['trust_careful','trust_depends','AGE','female','degree_1','degree_2','degree_3','degree_4',
      'marital_2','marital_3','marital_4','marital_5']
y=model.HAPPY.to_numpy(dtype=int);X=model[cols].to_numpy(dtype=float)
w=model.WTSSNRPS.to_numpy(float)

# Ordered response: 1=very happy, 2=pretty happy, 3=not too happy.
# P(Y<=j|x) = logistic(cut_j - xB). Positive coefficients indicate less happiness.
def ordinal_probability_scores(param,xx, yy):
    a,delta=param[:2];t1=a;t2=a+np.exp(np.clip(delta,-15,15))
    z=xx @ param[2:]
    f1=expit(np.clip(t1-z,-35,35));f2=expit(np.clip(t2-z,-35,35))
    g1=f1*(1-f1);g2=f2*(1-f2)
    p=np.where(yy==1,f1,np.where(yy==2,f2-f1,1-f2))
    p=np.maximum(p,1e-12)
    dt1=np.where(yy==1,g1/p,np.where(yy==2,-g1/p,0))
    dt2=np.where(yy==3,-g2/p,np.where(yy==2,g2/p,0))
    dz=np.where(yy==1,-g1/p,np.where(yy==2,(g1-g2)/p,g2/p))
    sc=np.column_stack([dt1+dt2,dt2*np.exp(np.clip(delta,-15,15)),xx*dz[:,None]])
    return p,sc

def fit_ordinal(xx,yy,ww,slower=False):
    p1=np.sum(ww*(yy==1))/sum(ww);p2=np.sum(ww*(yy<=2))/sum(ww)
    th1=logit(p1);th2=logit(p2)
    initial=np.r_[th1,np.log(th2-th1),np.zeros(xx.shape[1])]
    def fun(par):
        p,s=ordinal_probability_scores(par,xx,yy)
        return -np.dot(ww,np.log(p)),-np.sum(s*ww[:,None],axis=0)
    r=minimize(fun,initial,method='L-BFGS-B',jac=True,options={'maxiter':1800,'ftol':2e-12,'gtol':1e-7})
    if not r.success and not (np.linalg.norm(r.jac)<1e-3):
        raise RuntimeError('Ordinal fit failed: '+r.message)
    if not np.all(np.isfinite(r.x)):raise RuntimeError('Ordinal parameters nonfinite')
    return r

fit=fit_ordinal(X,y,w)
param=fit.x
p,s=ordinal_probability_scores(param,X,y)
score_per_case=s*w[:,None]
def local_jacobian(fun, at, step=1e-4):
    """Central-difference Jacobian of weighted score; avoids private SciPy APIs."""
    cols=[]
    for j in range(len(at)):
        plus=at.copy(); minus=at.copy()
        plus[j]+=step;minus[j]-=step
        cols.append((fun(plus)-fun(minus))/(2*step))
    return np.column_stack(cols)
H=-local_jacobian(lambda t:(ordinal_probability_scores(t,X,y)[1]*w[:,None]).sum(axis=0),param)
H=(H+H.T)/2
if np.linalg.eigvalsh(H).min()<=0:raise RuntimeError('Non-positive fitted Hessian')
Hinv=np.linalg.inv(H)

def design_sandwich(contribs, sampling):
    frame=pd.DataFrame(contribs);frame['stratum']=sampling['VSTRAT'].to_numpy();frame['psu']=sampling['VPSU'].to_numpy()
    group=frame.groupby(['stratum','psu'],sort=False).sum(numeric_only=True).reset_index()
    J=contribs.shape[1];meat=np.zeros((J,J));lonely=[]
    for strat,g in group.groupby('stratum',sort=False):
        mat=g.iloc[:,2:].to_numpy(float)
        if len(mat)<2:lonely.append(strat);continue
        mat-=mat.mean(axis=0)
        meat+=len(mat)/(len(mat)-1)*(mat.T@mat)
    if lonely:raise RuntimeError(f'Lonely PSUs in model domain: {lonely}')
    return meat

vcov=Hinv@design_sandwich(score_per_case,model)@Hinv
vcov=(vcov+vcov.T)/2
stderr=np.sqrt(np.diag(vcov))
tcrit=stats.t.ppf(.975,DF)
rows=[]
for i,name in enumerate(cols,2):
    coef=float(param[i]);se=float(stderr[i]);z=coef/se
    rows.append(dict(term=name,coefficient=coef,design_se=se,odds_ratio=np.exp(coef),
                     or_ci_low=np.exp(coef-tcrit*se),or_ci_high=np.exp(coef+tcrit*se),
                     p_t67=2*stats.t.sf(abs(z),DF)))
result=pd.DataFrame(rows)
result.to_csv(DIR/'table_03_weighted_ordinal_model.csv',index=False,float_format='%.6f')
trust_c=result[result.term.isin(['trust_careful','trust_depends'])]
W=param[2:4]@np.linalg.inv(vcov[2:4,2:4])@param[2:4]
p_joint=float(stats.f.sf(W/2,2,DF))

unw=fit_ordinal(X,y,np.ones(len(y)))
sm_fit=OrderedModel(y, model[cols].reset_index(drop=True),distr='logit').fit(method='bfgs',disp=False,maxiter=1500)
check=float(np.max(np.abs(sm_fit.params.loc[cols].values-unw.x[2:])))
if check>0.025:raise AssertionError(f'Unweighted ordinal regression did not replicate statsmodels: {check:.3g}')

comparison=[]
for name,cols2,weight in [
    ('Primary WTSSNRPS',cols,'WTSSNRPS'),
    ('Alternative weight WTSSPS',cols,'WTSSPS'),
    ('Interview mode adjusted',cols+[f'mode_{j}' for j in [2,3,4]],'WTSSNRPS'),
    ('Unweighted benchmark',cols,None),
]:
    ww=model[weight].to_numpy() if weight else np.ones(len(model))
    f=fit_ordinal(model[cols2].to_numpy(float),y,ww)
    comparison.append({'specification':name,'n':len(model),'or_trust_careful':np.exp(f.x[2]),
                       'or_trust_depends':np.exp(f.x[3]),'converged':f.success or np.linalg.norm(f.jac)<1e-3})
pd.DataFrame(comparison).to_csv(DIR/'table_04_ordinal_sensitivity.csv',index=False,float_format='%.6f')

thresholds=[]
for cutpoint in [2,3]:
    Y=(y>=cutpoint).astype(int)
    fitb=sm.GLM(Y,sm.add_constant(X),family=sm.families.Binomial(),freq_weights=w).fit()
    thresholds.append({'threshold':f'HAPPY >= {cutpoint}',
                       'coef_trust_careful':float(fitb.params[1]),'coef_trust_depends':float(fitb.params[2]),
                       'or_trust_careful':np.exp(fitb.params[1]),'or_trust_depends':np.exp(fitb.params[2])})
pd.DataFrame(thresholds).to_csv(Q/'proportional_odds_threshold_diagnostic.csv',index=False,float_format='%.6f')

inst=d[d.CONFINAN.notna()&d.CONEDUC.notna()].copy()
assert len(inst)==2142
wt=inst.WTSSNRPS.to_numpy()
instrows=[]
for category,name in [(1,'Great deal'),(3,'Hardly any')]:
    a=(inst.CONFINAN==category).to_numpy(int); b=(inst.CONEDUC==category).to_numpy(int)
    rate_a=np.dot(wt,a)/wt.sum();rate_b=np.dot(wt,b)/wt.sum();difference=rate_a-rate_b
    qq=d[['VSTRAT','VPSU']].copy()
    e=(d.CONFINAN.notna()&d.CONEDUC.notna()).to_numpy()
    a_all=(d.CONFINAN==category).to_numpy(int);b_all=(d.CONEDUC==category).to_numpy(int)
    influence=np.where(e,d.WTSSNRPS.to_numpy()*((a_all-b_all)-difference)/wt.sum(),0.)
    F=pd.DataFrame({'stratum':d.VSTRAT,'psu':d.VPSU,'influence':influence})
    sums=F.groupby(['stratum','psu']).influence.sum().reset_index()
    var=0
    for _,g in sums.groupby('stratum'):
        u=g.influence.to_numpy();var+=len(u)/(len(u)-1)*np.sum((u-u.mean())**2)
    se=np.sqrt(var)
    instrows.append(dict(category=name,n_pair=len(inst),finance_pct=100*rate_a,education_pct=100*rate_b,
                         finance_minus_education_pp=100*difference,ci95_low_pp=100*(difference-tcrit*se),
                         ci95_high_pp=100*(difference+tcrit*se),p_t67=2*stats.t.sf(abs(difference/se),DF)))
insttab=pd.DataFrame(instrows);insttab.to_csv(DIR/'table_05_institutional_confidence_paired.csv',index=False,float_format='%.6f')

manifest={'source':'NORC GSS 2024 SPSS Release 3a','source_sha256':hashlib.sha256(SAV.read_bytes()).hexdigest(),
          'parser':parser,'n_full_source':len(raw),'n_core':len(d),'n_trust_valid':int(d.TRUST.notna().sum()),
          'n_happy_valid':int(d.HAPPY.notna().sum()),'n_trust_happy_joint':int((d.TRUST.notna()&d.HAPPY.notna()).sum()),
          'n_primary_model':len(model),'n_paired_institutional_comparison':len(inst),
          'n_strata':int(len(strata)),'n_psus':int(strata.sum()),'design_df':DF,
          'weight':'WTSSNRPS','ordinal_ordering':'HAPPY 1 very, 2 pretty, 3 not too',
          'ordinal_log_likelihood':-float(fit.fun),
          'ordinal_p_value_joint_trust_f_approx':p_joint,
          'unweighted_optimizer_vs_statsmodels_max_absolute_beta_diff':check,
          'flags':['R survey independent replication unavailable in present runtime',
                   'IBM SPSS executable and .spv output unavailable in present runtime',
                   'Mode sensitivity and questionnaire routing constrain population interpretation',
                   'Custom compressed SAV numeric parser verified by NORC codebook counts but not pyreadstat/haven parity in this runtime']}
(Q/'final_model_quality_summary.json').write_text(json.dumps(manifest,indent=2))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(9,4))
terms=['trust_careful','trust_depends'];names=['Cannot be too careful','Depends / conditional trust']
yy=[1,0]
subset=result.set_index('term').loc[terms]
for loc,(idx,row) in enumerate(subset.iterrows()):
    ax.plot([row.or_ci_low,row.or_ci_high],[yy[loc]]*2,color='#21456b',lw=2)
    ax.scatter(row.odds_ratio,yy[loc],color='#7b52e7',s=90,zorder=3)
ax.axvline(1,color='#69788f',ls='--',lw=1)
ax.set_yticks(yy,names);ax.set_xscale('log');ax.set_xlim(.15,4.5)
ax.set_xlabel('Adjusted proportional-odds ratio (95% design-based CI)')
ax.set_title('Association between interpersonal trust and lower happiness',loc='left',fontweight='bold',pad=15)
ax.text(.01,-.32,'Core GSS 2024; weighted ordinal model; n=911 | No causal interpretation',transform=ax.transAxes,fontsize=8)
fig.tight_layout();fig.savefig(DIR/'figure_03_ordinal_adjusted_odds.png',dpi=230,bbox_inches='tight');plt.close(fig)

print(json.dumps(manifest,indent=2))
print(result.loc[result.term.str.contains('trust')].to_string(index=False))
print(insttab.to_string(index=False))
print('Thresholds',thresholds)
print('Sensitivities',comparison)
