"""Original optical analysis kernels; coordinate conventions are unchanged."""
import itertools
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation,Slerp
P=['relative_x_mm','relative_y_mm','relative_z_mm']
Q=['q_x','q_y','q_z','q_w']
def runs(g,m):
 ix=np.flatnonzero(m);t=g.optical_s.to_numpy();fr=g.frame.to_numpy()
 if not len(ix):return []
 cut=np.flatnonzero((np.diff(ix)>1)|(np.diff(fr[ix])!=1)|(np.diff(t[ix])<=0)|(np.diff(t[ix])>1.5/90))+1
 return [z for z in np.split(ix,cut) if len(z)>1]

def meanR(q):return Rotation.from_quat(q).mean()

def derivatives(g,m):
 t=g.optical_s.to_numpy();p=g[P].to_numpy();q=g[Q].to_numpy();v=np.full((len(g),3),np.nan);w=np.full_like(v,np.nan)
 for ii in runs(g,m):
  r=Rotation.from_quat(q[ii]);tt=t[ii];v[ii]=np.gradient(p[ii],tt,axis=0)
  if len(ii)>=3:
   w[ii[1:-1]]=(r[2:]*r[:-2].inv()).as_rotvec(degrees=True)/(tt[2:]-tt[:-2])[:,None]
  w[ii[0]]=(r[1]*r[0].inv()).as_rotvec(degrees=True)/(tt[1]-tt[0]);w[ii[-1]]=(r[-1]*r[-2].inv()).as_rotvec(degrees=True)/(tt[-1]-tt[-2])
 return v,w

def profile(g,m,ref,start,end,omega):
 t=g.optical_s.to_numpy();ix=np.flatnonzero(m&(t>=start)&(t<=end));p=g[P].to_numpy();q=g[Q].to_numpy()
 pp=np.full((101,3),np.nan);qq=np.full((101,4),np.nan);ww=np.full((101,3),np.nan);count=np.zeros(101,int);p0,r0=ref
 if len(ix)<5:return dict(p=pp,q=qq,w=ww,n=count)
 bid=np.minimum(((t[ix]-start)/(end-start)*101).astype(int),100)
 for k in range(101):
  zi=ix[bid==k]
  if len(zi)<5:continue
  pp[k]=np.mean(p[zi]-p0,axis=0);qq[k]=(r0.inv()*meanR(q[zi])).as_quat();ww[k]=r0.inv().apply(np.nanmean(omega[zi],axis=0));count[k]=len(zi)
 return dict(p=pp,q=qq,w=ww,n=count)

def cycle_profiles(g,b,m,p0,r0,v,w,cycles):
 rid=g.record_id.iloc[0];t=g.optical_s.to_numpy();p=g[P].to_numpy();q=g[Q].to_numpy();fr=g.frame.to_numpy();curves=[];qr=[];out=[];window=[]
 for c in cycles[cycles.record_id==rid].itertuples():
  ii=np.flatnonzero((fr>=c.start_frame)&(fr<=c.end_frame));reason='PASS'
  if len(ii)<10:reason='INSUFFICIENT_NATIVE_SAMPLES'
  elif t[ii[0]]<b.start_optical_s or t[ii[-1]]>b.end_optical_s:reason='NON_TASK_BOUNDARY'
  elif not g.valid.to_numpy(bool)[ii].all():reason='INVALID_TRACKING'
  elif not m[ii].all():reason='PAUSE_OR_NON_TASK'
  elif np.any(np.diff(fr[ii])!=1) or np.any(np.diff(t[ii])>1.5/90):reason='NATIVE_GAP'
  out.append(dict(record_id=rid,cycle_index=c.cycle_index,start_frame=c.start_frame,end_frame=c.end_frame,start_s=t[ii[0]] if len(ii) else np.nan,end_s=t[ii[-1]] if len(ii) else np.nan,status=reason))
  if reason!='PASS':continue
  grid=np.linspace(t[ii[0]],t[ii[-1]],101);pp=np.column_stack([np.interp(grid,t[ii],p[ii,k]-p[ii[0],k]) for k in range(3)]);vv=np.column_stack([np.interp(grid,t[ii],v[ii,k]) for k in range(3)]);ww=np.column_stack([np.interp(grid,t[ii],w[ii,k]) for k in range(3)]);rr=Rotation.from_quat(q[ii]);rc=rr[0].inv()*Slerp(t[ii],rr)(grid)
  curves.append(np.column_stack([pp,vv,ww]));qr.append(rc.as_quat());window.append((t[ii[0]],t[ii[-1]],t[ii[-1]]-t[ii[0]]))
 if not curves:return None,pd.DataFrame(out),[]
 a=np.stack(curves);rq=np.stack(qr);rows=[];meanq=[]
 for j in range(101):
  rm=meanR(rq[:,j]);meanq.append(rm.as_quat());rv=rm.as_rotvec(degrees=True);dev=(rm.inv()*Rotation.from_quat(rq[:,j])).as_rotvec(degrees=True)
  row=dict(record_id=rid,phase=j/100,cycles=len(a))
  for k,name in enumerate(['dx_mm','dy_mm','dz_mm','vx_mm_s','vy_mm_s','vz_mm_s','wx_spatial_deg_s','wy_spatial_deg_s','wz_spatial_deg_s']):row[name+'_mean']=np.mean(a[:,j,k]);row[name+'_q25']=np.percentile(a[:,j,k],25);row[name+'_q75']=np.percentile(a[:,j,k],75)
  for k,c in enumerate('xyz'):row['rotvec_'+c+'_mean_deg']=rv[k];row['rotvec_'+c+'_q25_deg']=rv[k]+np.percentile(dev[:,k],25);row['rotvec_'+c+'_q75_deg']=rv[k]+np.percentile(dev[:,k],75)
  for k,c in enumerate('xyzw'):row['cycle_mean_q_'+c]=rm.as_quat()[k]
  rows.append(row)
 return dict(array=a,quat=rq,meanq=np.array(meanq),times=np.array(window),task_start=b.start_optical_s,task_end=b.end_optical_s),pd.DataFrame(out),rows

def formal_native_events(g,m,w):
 """Reapply native signed-PC1 positive-peak events to formal data, including t<15.

 Reuses historical distance/prominence/duration limits; no 1/f forcing. PCA is
 fitted on formal Tracker0-spatial signed angular velocity only.
 """
 from scipy.signal import find_peaks
 t=g.optical_s.to_numpy();fr=g.frame.to_numpy();finite=m&np.isfinite(w).all(axis=1);x=w[finite]
 if len(x)<50:return pd.DataFrame(columns=['record_id','cycle_index','start_frame','end_frame']),{}
 center=np.mean(x,axis=0);cov=np.cov(x-center,rowvar=False);values,axes=np.linalg.eigh(cov);axis=axes[:,-1]
 # PCA axes are unoriented. Fix the sign to a single physical reference (+Y
 # of Tracker0), never to whichever loading happens to be largest in a record.
 # The latter can reverse phase when |Y| and |Z| swap rank across repetitions.
 polarity_support=float(abs(axis[1]))
 if axis[1]<0:axis=-axis
 proj=(w-center)@axis;sd=np.std(proj[finite]);prom=max(.35*sd,1.);events=[];rejected=0
 for ii in runs(g,finite):
  peaks,_=find_peaks(proj[ii],distance=20,prominence=prom)
  for a,b in zip(peaks[:-1],peaks[1:]):
   ia,ib=ii[a],ii[b];duration=t[ib]-t[ia]
   if duration<.2 or duration>.7:rejected+=1;continue
   events.append(dict(record_id=g.record_id.iloc[0],cycle_index=len(events)+1,start_frame=int(fr[ia]),end_frame=int(fr[ib]),start_s=t[ia],end_s=t[ib],duration_s=duration,frequency_Hz=1/duration))
 diag=dict(record_id=g.record_id.iloc[0],event_version='formal_native_cycle_v2',axis_polarity='Tracker0 +Y dot PC1 positive',axis_polarity_condition_abs_Y=polarity_support,axis_polarity_supported=bool(polarity_support>=.2),PCA_axis_Tracker0_x=axis[0],PCA_axis_Tracker0_y=axis[1],PCA_axis_Tracker0_z=axis[2],PC1_variance_fraction=float(values[-1]/values.sum()),prominence_deg_s=prom,candidate_cycles=len(events),rejected_duration_intervals=rejected)
 start,end=t[np.flatnonzero(m)[0]],t[np.flatnonzero(m)[-1]];edges=np.linspace(start,end,4)
 for k,name in enumerate(['early','middle','late']):
  z=finite&(t>=edges[k])&(t<=edges[k+1])
  if z.sum()<50:continue
  _,u=np.linalg.eigh(np.cov(w[z],rowvar=False));cos=np.clip(abs(axis@u[:,-1]),0,1);diag[name+'_axis_deviation_deg']=float(np.rad2deg(np.arccos(cos)))
 return pd.DataFrame(events),diag

def compare_pairs(profiles,meta,scenario):
 rows=[]
 for a,b in itertools.combinations(profiles,2):
  A,B=meta.loc[a],meta.loc[b];kind=None
  if A.expert==B.expert and A.material==B.material and A.block==B.block and A['repeat']!=B['repeat']:kind='repeat'
  elif A.expert!=B.expert and A.material==B.material and A.block==B.block and A['repeat']==B['repeat']:kind='expert'
  elif A.expert==B.expert and A.block==B.block and A['repeat']==B['repeat'] and A.material!=B.material:kind='material'
  if kind is None:continue
  pa,pb=profiles[a],profiles[b];m=(pa['n']>=5)&(pb['n']>=5)&np.isfinite(pa['q']).all(axis=1)&np.isfinite(pb['q']).all(axis=1)
  if m.sum()<50:continue
  theta=np.rad2deg((Rotation.from_quat(pa['q'][m]).inv()*Rotation.from_quat(pb['q'][m])).magnitude());dp=pa['p'][m]-pb['p'][m]
  row=dict(comparison=kind,record_a=a,record_b=b,expert_a=A.expert,expert_b=B.expert,material_a=A.material,material_b=B.material,block=int(A.block),repeat_a=int(A['repeat']),repeat_b=int(B['repeat']),reference_scenario=scenario,common_bins=int(m.sum()),position_shape_RMS_mm=float(np.sqrt(np.mean(np.sum(dp*dp,axis=1)))),SO3_RMS_deg=float(np.sqrt(np.mean(theta**2))),SO3_P95_deg=float(np.percentile(theta,95)))
  for k,c in enumerate('xyz'):row[f'position_{c}_RMS_mm']=float(np.sqrt(np.mean(dp[:,k]**2)))
  rows.append(row)
 return rows
