import sys, re, math, numpy as np, pandas as pd
X_,J_=sys.argv[1],sys.argv[2]
sys.argv=['x',X_,J_,'primary']; src=open('reproduce.py').read(); exec(src.split("# ---------------- pair records")[0])
from scipy.optimize import minimize
# ---- curve attributes from the parameter sheet
def num(v):
    m=re.findall(r'[-+]?\d*\.?\d+',str(v)); return float(m[0]) if m and str(v)!='nan' else np.nan
P=par.set_index('id')
attr=pd.DataFrame({'dir':P['direciton of specimen \n(°)'].map(num),'freq':P['frequency \n(Hz)'].map(num),'load':P['load control'].astype(str).str.lower().replace({'nan':np.nan,'':np.nan}),
 'shape':P['specimens description'].astype(str).str.lower().str.extract(r'^(circular|rectangular|square)')[0],'amm':P['AM machine'].astype(str).str.lower().str.strip().replace({'nan':np.nan,'':np.nan})})
attr['fband']=pd.cut(attr.freq,[0,30,100,1e9],labels=['<=30','30-100','>100'])
print('attribute availability over curves in the primary set (non-missing %):'); cur=D.drop_duplicates('id').id
for c in ['dir','freq','load','shape','amm']: print(f'  {c:6s} {100*attr.loc[cur,c].notna().mean():.0f}%')
def eq(a,b): return (pd.notna(a) and pd.notna(b) and a==b)
LEVELS=['L0 none (same alloy, R, L-PBF, uniaxial)','L1 + same surface class, HIP, high-T heat treatment','L2 + same build direction','L3 + same specimen shape and load control','L4 + same frequency band','L5 + same AM machine']
rec=[]
for mat in ALLOYS:
    for R in (-1,0.1):
        G=D[(D.mat==mat)&(D.R==R)].reset_index(drop=True)
        if G.DOI.nunique()<6: continue
        y=np.log10(G.N.values); x=np.log10(G.S.values/np.median(G.S.values)); C=G[['asb','hip','hthi']].values.astype(float); Cm=C.copy(); C=C[:,C.std(0)>0]
        X=np.column_stack([np.ones(len(G)),x,C]); cens=G.ro.values.astype(bool); doi=G.DOI.values; cid=G.id.values
        curves={}
        for sid in np.unique(cid):
            m=(cid==sid)&~cens
            if m.sum()<6 or len(np.unique(np.round(x[m],3)))<3: continue
            curves[sid]=(doi[cid==sid][0],np.column_stack([np.ones(m.sum()),x[m]]),y[m],x[m].min(),x[m].max())
        for B in np.unique(doi):
            tr=doi!=B; te=doi==B; bc,sc=censored_lognormal_em(X[tr],y[tr],cens[tr])
            for sid,(A,Xs,ys,lo,hi) in curves.items():
                if A==B: continue
                n=len(ys); bs,sd=ols_fail(Xs,ys); XtXi=np.linalg.inv(Xs.T@Xs); asrc=attr.loc[sid]; csrc=Cm[cid==sid][0]
                for tid in np.unique(cid[te]):
                    sel=(cid==tid)&(x>=lo)&(x<=hi)
                    if not sel.any(): continue
                    x0=np.column_stack([np.ones(sel.sum()),x[sel]]); q=np.einsum('ij,jk,ik->i',x0,XtXi,x0); mu=x0@bs
                    pred=mu-tdist.ppf(0.95,n-2)*sd*np.sqrt(1+q); pool=X[sel]@bc+Z05*sc
                    at=attr.loc[tid]; ct=Cm[cid==tid][0]
                    l1=bool((csrc==ct).all()); l2=l1 and eq(asrc['dir'],at['dir']); l3=l2 and eq(asrc['shape'],at['shape']) and eq(asrc['load'],at['load'])
                    l4=l3 and eq(asrc['fband'],at['fband']); l5=l4 and eq(asrc['amm'],at['amm'])
                    lvl=0+l1+l2+l3+l4+l5
                    rec.append((mat,R,A,sid,B,tid,int(sel.sum()),int(np.sum(~cens[sel]&(y[sel]<pred))),int(np.sum(~cens[sel]&(y[sel]<pool))),lvl))
R_=pd.DataFrame(rec,columns=['mat','R','A','sid','B','tid','n','b_single','b_pool','lvl']); R_.to_pickle('comparability_records.pkl')
def weighted(pr,wA,wB):
    perB={}
    for (A,B),r in pr.items():
        w=wA.get(A,0)
        if w: s=perB.setdefault(B,[0,0]); s[0]+=w*r; s[1]+=w
    num_=den=0
    for B,(s,w) in perB.items():
        wb=wB.get(B,0)
        if wb and w: num_+=wb*s/w; den+=wb
    return num_/den if den else np.nan
def level_table(minlvl,col,rs):
    q=R_[R_.lvl>=minlvl]; g=q.groupby(['mat','R','A','sid','B'])[[col,'n']].sum().reset_index(); g['r']=g[col]/g.n
    p1=g.groupby(['mat','R','A','B']).r.mean().reset_index().groupby(['A','B']).r.mean(); pr=p1.to_dict(); As=sorted({a for a,b in pr}); Bs=sorted({b for a,b in pr})
    est=weighted(pr,{a:1 for a in As},{b:1 for b in Bs}); bs=[]
    for _ in range(1000):
        ca=np.bincount(rs.integers(0,len(As),len(As)),minlength=len(As)); cb=np.bincount(rs.integers(0,len(Bs),len(Bs)),minlength=len(Bs)); bs.append(weighted(pr,dict(zip(As,ca)),dict(zip(Bs,cb))))
    return est,*np.percentile(bs,[2.5,97.5]),len(q),len(pr),len(As),len(Bs),int(q.n.sum())
rs=np.random.default_rng(5)
print('\nSINGLE-STUDY 5% prediction bound breach rate under progressively strict source-target matching (primary set; equal-weight per paper pair)')
print(f'{"matching level":55s} {"breach":>7s} {"95% CI":>14s} | curve pairs  paper pairs  src  tgt  specimens | pooled bound on same specimens')
for L,name in enumerate(LEVELS):
    e=level_table(L,'b_single',rs); p=level_table(L,'b_pool',rs)
    print(f'{name:55s} {100*e[0]:6.1f}% [{100*e[1]:5.1f},{100*e[2]:5.1f}] | {e[3]:6d} {e[4]:11d} {e[5]:5d} {e[6]:4d} {e[7]:9d} | {100*p[0]:5.1f}% [{100*p[1]:4.1f},{100*p[2]:5.1f}]')
