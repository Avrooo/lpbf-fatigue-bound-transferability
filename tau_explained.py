import sys, re, math, numpy as np, pandas as pd
X_,J_=sys.argv[1],sys.argv[2]
sys.argv=['x',X_,J_,'primary']; src=open('reproduce.py').read(); exec(src.split("# ---------------- pair records")[0])
from scipy.optimize import minimize
def num(v):
    m=re.findall(r'[-+]?\d*\.?\d+',str(v)); return float(m[0]) if m and str(v)!='nan' else np.nan
P=par.set_index('id')
att=pd.DataFrame({'dir':P['direciton of specimen \n(°)'].map(num),'freq':P['frequency \n(Hz)'].map(num),'load':P['load control'].astype(str).str.lower(),'shape':P['specimens description'].astype(str).str.lower().str.extract(r'^(circular|rectangular|square)')[0],'amm':P['AM machine'].astype(str).str.lower().str.strip()})
att['fband']=pd.cut(att.freq,[0,30,100,1e9],labels=['<=30','30-100','>100']).astype(object)
D2=D.merge(att,left_on='id',right_index=True,how='left')
GHx,GHw=np.polynomial.hermite_e.hermegauss(21); GHw=GHw/GHw.sum()
def lpf(y,mu,s,c):
    z=(y-mu)/s; return np.where(c,norm.logsf(z),norm.logpdf(z)-np.log(s))
def nll(p,y,X,c,u,K):
    k=X.shape[1]; b=p[:k]; s=math.exp(p[k]); t=math.exp(p[k+1]); l=lpf(y[:,None],(X@b)[:,None]+t*GHx[None,:],s,c[:,None]); L=np.zeros((K,len(GHx))); np.add.at(L,u,l); mx=L.max(1,keepdims=True)
    return -(mx[:,0]+np.log((np.exp(L-mx)*GHw).sum(1))).sum()
def fit_re(y,X,c,g):
    u=np.unique(g,return_inverse=True)[1]; K=u.max()+1; b0,s0=censored_lognormal_em(X,y,c)
    r=minimize(nll,np.r_[b0,math.log(s0*0.9),math.log(0.5)],args=(y,X,c,u,K),method='L-BFGS-B'); k=X.shape[1]; return math.exp(r.x[k]),math.exp(r.x[k+1])
def dummies(col,G,minn=0):
    v=G[col].astype(object).where(G[col].notna(),'unk'); v=v.where(v!='nan','unk')
    if minn: vc=v.value_counts(); v=v.where(v.map(vc)>=minn,'other')
    d=pd.get_dummies(v,drop_first=True).astype(float); return d.loc[:,d.std()>0]
SETS=[('none',[]),('+ surface/HIP/heat treatment',['proc']),('+ build direction',['proc','dir']),('+ specimen shape, load control',['proc','dir','shape','load']),('+ frequency band',['proc','dir','shape','load','fband']),('+ AM machine',['proc','dir','shape','load','fband','amm'])]
rows=[]
for mat in ALLOYS:
    for R in (-1,0.1):
        G=D2[(D2.mat==mat)&(D2.R==R)].reset_index(drop=True)
        if G.DOI.nunique()<6: continue
        y=np.log10(G.N.values); x=np.log10(G.S.values/np.median(G.S.values)); cens=G.ro.values.astype(bool); g=G.DOI.values
        for name,cols in SETS:
            parts=[np.ones(len(G)),x]
            for c_ in cols:
                if c_=='proc':
                    C=G[['asb','hip','hthi']].values.astype(float); C=C[:,C.std(0)>0]; parts+= [C[:,i] for i in range(C.shape[1])]
                else:
                    d_=dummies('dir' if c_=='dir' else c_,G if c_!='dir' else G.assign(dir=G['dir'].astype(object)),minn=(40 if c_=='amm' else 0)); parts+=[d_.values[:,i] for i in range(d_.shape[1])]
            X=np.column_stack(parts); s,t=fit_re(y,X,cens,g); rows.append((mat,R,name,X.shape[1]-2,len(G),s,t))
T=pd.DataFrame(rows,columns=['mat','R','covariates','n_cov','n','s_within','tau_paper']); T['ICC']=T.tau_paper**2/(T.tau_paper**2+T.s_within**2)
w=T.groupby('covariates',sort=False).apply(lambda q: pd.Series({'mean_covariate_cols':q.n_cov.mean(),'specimen-weighted tau':np.average(q.tau_paper,weights=q.n),'specimen-weighted within-paper s':np.average(q.s_within,weights=q.n),'ICC (weighted)':np.average(q.ICC,weights=q.n)}))
print('\nBETWEEN-PAPER SD (tau, log10 life) after adding recorded covariates (random-intercept censored model, 8 alloy-R groups):'); print(w.round(3).to_string())
t0=w.iloc[0]['specimen-weighted tau']; print('\nShare of between-paper variance removed vs no covariates:'); print(((1-(w['specimen-weighted tau']/t0)**2)*100).round(1).to_string())
print('\nPer group tau, none -> all recorded covariates:'); a=T[T.covariates=='none'].set_index(['mat','R']).tau_paper; b=T[T.covariates=='+ AM machine'].set_index(['mat','R']).tau_paper; print(pd.DataFrame({'none':a,'all':b,'n_cov_all':T[T.covariates=='+ AM machine'].set_index(['mat','R']).n_cov}).round(2).to_string())
