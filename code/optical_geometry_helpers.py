"""Original helper functions retained for interpreting raw optical exports.
Not a complete raw-to-mask pipeline and not called by reproduce.py.
"""
from pathlib import Path
import csv
import numpy as np
import pandas as pd
def read_pose_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.reader(f))
    st=next(i for i,r in enumerate(rows) if len(r)>1 and r[1].strip()=='Timestamp')
    cols=[x.strip() for x in rows[st]]
    d=pd.DataFrame(rows[st+1:],columns=cols).apply(pd.to_numeric,errors='coerce')
    pos=d[['XToGlobal1','YToGlobal1','ZToGlobal1']].to_numpy(float)
    q=d[['QxToGlobal1','QyToGlobal1','QzToGlobal1','QwToGlobal1']].to_numpy(float)
    return d.iloc[:,0].to_numpy(int),d.Timestamp.to_numpy(float),pos,q

def rigid_fit_subset(template,points):
    """Batched proper Kabsch fit: template-local mm -> optical-world mm."""
    xc=template-template.mean(axis=1,keepdims=True); yc=points-points.mean(axis=1,keepdims=True)
    h=np.einsum('nki,nkj->nij',xc,yc)
    u,_,vh=np.linalg.svd(h)
    r=vh.transpose(0,2,1)@u.transpose(0,2,1)
    det=np.linalg.det(r); bad=det<0
    if bad.any():
        vh[bad,-1,:]*=-1
        r=vh.transpose(0,2,1)@u.transpose(0,2,1)
    t=points.mean(axis=1)-np.einsum('nij,nj->ni',r,template.mean(axis=1))
    pred=np.einsum('nij,nkj->nki',r,template)+t[:,None,:]
    rms=np.sqrt(np.mean(np.sum((pred-points)**2,axis=2),axis=1))
    return r,t,rms

def fit_markers(raw,template):
    n=len(raw); obs=np.isfinite(raw).all(axis=2)&np.any(raw!=0,axis=2); count=obs.sum(axis=1)
    outR=np.full((n,3,3),np.nan); outp=np.full((n,3),np.nan); rms=np.full(n,np.nan); fitok=np.zeros(n,bool)
    for miss in range(5):
        if miss==4: ids=np.flatnonzero(count==4); keep=np.arange(4)
        else: ids=np.flatnonzero((count==3)&(~obs[:,miss])); keep=np.array([k for k in range(4) if k!=miss])
        if len(ids)==0: continue
        r,t,e=rigid_fit_subset(np.broadcast_to(template[keep],(len(ids),len(keep),3)),raw[ids][:,keep])
        outR[ids]=r;outp[ids]=t;rms[ids]=e;fitok[ids]=True
    # Triangle non-degeneracy and per-marker edge-pattern residual.
    tri=np.zeros(n)
    for a,b,c in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
        e1=raw[:,b]-raw[:,a];e2=raw[:,c]-raw[:,a]
        sn=np.linalg.norm(np.cross(e1,e2),axis=1)/np.maximum(np.linalg.norm(e1,axis=1)*np.linalg.norm(e2,axis=1),1e-12)
        ok=obs[:,[a,b,c]].all(axis=1);tri=np.maximum(tri,np.where(ok,sn,0))
    edges=list(__import__('itertools').combinations(range(4),2))
    edge_err=np.full(n,np.nan)
    for i in np.flatnonzero(count>=3):
        ids=[(a,b) for a,b in edges if obs[i,a] and obs[i,b]]
        if ids:
            dd=np.array([np.linalg.norm(raw[i,a]-raw[i,b])-np.linalg.norm(template[a]-template[b]) for a,b in ids])
            edge_err[i]=np.max(np.abs(dd))
    return dict(observed=obs,count=count,R=outR,p=outp,fit_rms=rms,fit_valid=fitok,triangle_sine=tri,edge_error=edge_err)
