"""Current-paper scalar decomposition and stratified response analysis.
Read-only frozen inputs; no signal filtering, event detection, alignment or mask edits.
Run through reproduce.py. Outputs are written to evidence/quantitative in its output directory.
"""
from pathlib import Path
import csv, json, hashlib, math, itertools
from collections import defaultdict, Counter
import numpy as np
import os
R=Path(os.environ['TBME_RUN_ROOT'])
Q=R/'evidence/quantitative'; Q.mkdir(exist_ok=True)
ORDER=['P1','P3','P2','P4']; CHANNELS=['Fx','Fy','Fz','Tx','Ty','Tz','wx','wy','wz']
INPUTS=['evidence/formal_ati_imu_phase_profiles.csv','evidence/formal_ati_imu_pair_distances.csv','evidence/formal_ati_imu_records.csv','evidence/temporal/formal_third_imu_profiles.csv',None,'evidence/optical_formal_record_features.csv',None,None,None]
def read(f):
    with (R/f).open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def write(name, rows):
    assert rows
    with (Q/name).open('w',encoding='utf-8',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def decompose(a,b,object_type='scalar'):
    """a,b map phase->scalar. Uniform weights on finite intersection, recompute m,a.
    Degenerate shape/rho and zero-distance fractions are undefined (empty CSV).
    Numerical near-zero criterion = 64*eps*max(1 unit, max|a|, max|b|).
    This is a numerical conditioning check, not a biological equivalence margin.
    No clipping of correlations, distances or fractions.
    """
    if object_type!='scalar':raise ValueError('Only signed scalar physical channels; not SO(3).')
    shared=sorted(set(a)&set(b)); phases=[p for p in shared if math.isfinite(a[p]) and math.isfinite(b[p])]
    base={'n_common':len(phases),'n_union':len(set(a)|set(b)),'n_shared_nonfinite':len(shared)-len(phases),'phase_ids':';'.join(f'{p:g}' for p in phases),'weight_rule':'uniform finite common phase points'}
    keys=['mean_a','mean_b','amplitude_a','amplitude_b','delta_mean_b_minus_a','delta_amplitude_b_minus_a','rho','shape_RMS_dimensionless','distance','distance_squared','center_MSE','amplitude_MSE','shape_MSE','center_fraction','amplitude_fraction','shape_fraction','identity_error','near_zero_tolerance']
    if len(phases)<2:return dict(base,**{k:'' for k in keys},status='insufficient_common_support')
    x=np.array([a[p] for p in phases]);y=np.array([b[p] for p in phases])
    mx,my=x.mean(),y.mean();cx=x-mx;cy=y-my
    ax,ay=np.sqrt(np.mean(cx**2)),np.sqrt(np.mean(cy**2))
    tol=64*np.finfo(float).eps*max(1.,np.max(np.abs(x)),np.max(np.abs(y)))
    d2=np.mean((x-y)**2); cm=(mx-my)**2; am=(ax-ay)**2
    rho=sd=''
    if ax>tol and ay>tol:
        sx=cx/ax;sy=cy/ay;rho=float(np.mean(sx*sy));sd=float(np.sqrt(np.mean((sx-sy)**2)))
        sm=ax*ay*sd**2
        assert abs(sm-2*ax*ay*(1-rho))<=1e-10*max(1.,d2)
        status='ok'
    else:
        # Algebraic continuation is valid even when standardized shapes are undefined.
        sm=2*(ax*ay-np.mean(cx*cy));status='shape_undefined_near_zero_amplitude'
    err=cm+am+sm-d2
    assert abs(err)<=1e-10*max(1.,d2),(d2,err)
    # Preserve any roundoff residual, never fabricate nonnegative fractions.
    fractions=[float(v/d2) for v in [cm,am,sm]] if d2>tol**2 else ['','','']
    if d2<=tol**2:status+=';zero_or_numerically_zero_distance_fractions_undefined'
    return dict(base,mean_a=float(mx),mean_b=float(my),amplitude_a=float(ax),amplitude_b=float(ay),delta_mean_b_minus_a=float(my-mx),delta_amplitude_b_minus_a=float(ay-ax),rho=rho,shape_RMS_dimensionless=sd,distance=float(math.sqrt(d2)),distance_squared=float(d2),center_MSE=float(cm),amplitude_MSE=float(am),shape_MSE=float(sm),center_fraction=fractions[0],amplitude_fraction=fractions[1],shape_fraction=fractions[2],identity_error=float(err),near_zero_tolerance=float(tol),status=status)

def main():
    # Inputs are frozen; this entry recomputes current-paper quantities only.
    rec={q['record_id']:q for q in read('evidence/formal_ati_imu_records.csv')};assert len(rec)==32
    def meta(rid):return {k:rec[rid][k] for k in ['expert','material','block','repeat']}
    profiles=defaultdict(dict);counts={}; units={}
    for q in read(INPUTS[0]):
        key=q['record_id'],q['channel'];p=float(q['phase']);assert p not in profiles[key]
        profiles[key][p]=float(q['mean']);counts[key,p]=int(q['n_cycles']);units[q['channel']]=q['unit']
    assert len(profiles)==32*9
    assert all(v>0 for v in counts.values())
    features=[]
    for (rid,ch),a in sorted(profiles.items()):
        finite={p:v for p,v in a.items() if math.isfinite(v)};d=decompose(finite,finite)
        for key,val in [('cycle_center_'+ch,d['mean_a']),('cycle_amplitude_'+ch,d['amplitude_a'])]:
            features.append(dict(record_id=rid,**meta(rid),feature=key,value=val,unit=units[ch],n_support=d['n_common'],support_kind='finite full-record mean-cycle phase bins',source=INPUTS[0]))
        if ch=='Fz':features.append(dict(record_id=rid,**meta(rid),feature='cycle_compression_center',value=-d['mean_a'],unit='N',n_support=d['n_common'],support_kind='negative mean signed Fz cycle',source=INPUTS[0]))
    for q in read(INPUTS[5]):
        rid=q['record_id']
        for feature,unit in [('position_x_range_p05_p95_mm','mm'),('position_y_range_p05_p95_mm','mm'),('position_z_range_p05_p95_mm','mm'),('orientation_excursion_p95_deg','deg')]:
            features.append(dict(record_id=rid,**meta(rid),feature=feature,value=float(q[feature]),unit=unit,n_support=int(q['formal_valid_frames']),support_kind='qualified native optical frames',source=INPUTS[5]))
    write('record_features.csv',features)
    pairs=[]
    for q in read(INPUTS[1]):
        if q['channel'] not in CHANNELS:raise ValueError(q)
        a,b,ch=q['record_a'],q['record_b'],q['channel'];d=decompose(profiles[a,ch],profiles[b,ch])
        assert abs(d['distance']-float(q['RMS']))<1e-9
        ma,mb=meta(a),meta(b)
        if q['comparison']=='expert':assert ma['expert']!=mb['expert'] and all(ma[k]==mb[k] for k in ['material','block','repeat'])
        if q['comparison']=='repeat':assert all(ma[k]==mb[k] for k in ['expert','material','block']) and ma['repeat']!=mb['repeat']
        nca=min(counts[(a,ch),p] for p in profiles[a,ch] if math.isfinite(profiles[a,ch][p]));ncb=min(counts[(b,ch),p] for p in profiles[b,ch] if math.isfinite(profiles[b,ch][p]))
        pairs.append(dict(pair_id=a+'|'+b,comparison=q['comparison'],record_a=a,record_b=b,channel=ch,unit=units[ch],squared_unit={'N':'N^2','N m':'(N m)^2','deg/s':'(deg/s)^2'}[units[ch]],material_a=ma['material'],material_b=mb['material'],block_a=ma['block'],block_b=mb['block'],repeat_a=ma['repeat'],repeat_b=mb['repeat'],n_cycles_a=nca,n_cycles_b=ncb,**d))
    assert len(pairs)==720
    write('pair_decomposition.csv',pairs)
    summ=[]
    for comp in ['expert','repeat','material']:
        for ch in CHANNELS:
            ss=[q for q in pairs if q['comparison']==comp and q['channel']==ch]
            out=dict(comparison=comp,channel=ch,unit=units[ch],squared_unit=ss[0]['squared_unit'],n_pairs=len(ss))
            for key in ['center_MSE','amplitude_MSE','shape_MSE','center_fraction','amplitude_fraction','shape_fraction','rho','shape_RMS_dimensionless']:
                vals=[q[key] for q in ss if q[key]!=''];out[key+'_median']=float(np.median(vals));out[key+'_min']=float(min(vals));out[key+'_max']=float(max(vals))
            out['largest_center_count']=sum(q['center_MSE']>max(q['amplitude_MSE'],q['shape_MSE']) for q in ss)
            out['largest_amplitude_count']=sum(q['amplitude_MSE']>max(q['center_MSE'],q['shape_MSE']) for q in ss)
            out['largest_shape_count']=sum(q['shape_MSE']>max(q['center_MSE'],q['amplitude_MSE']) for q in ss)
            summ.append(out)
    write('decomposition_summary.csv',summ)
    # Material features: average the two repeats arithmetically, retain both rows and trajectories.
    fmap={(q['record_id'],q['feature']):q for q in features}; fns=sorted(set(q['feature'] for q in features))
    cells=[];step=[];mono=[];endpoints=[];did=[]
    def sequence(vals):
        dif=np.diff(vals);tol=64*np.finfo(float).eps*max(1.,max(abs(x) for x in vals));sg=np.where(dif>tol,1,np.where(dif < -tol,-1,0))
        return ('increasing' if np.all(sg>=0) and np.any(sg>0) else 'decreasing' if np.all(sg<=0) and np.any(sg<0) else 'flat' if np.all(sg==0) else 'mixed'),sg,tol
    def choose(e,b,m,k):return next(rid for rid,q in rec.items() if (q['expert'],q['block'],q['material'],q['repeat'])==(e,b,m,k))
    for f in fns:
        for e,b in itertools.product(['ExpertA','ExpertB'],['1','2']):
            for k in ['1','2','mean_of_two']:
                vals=[];ids=[]
                for m in ORDER:
                    rr=[choose(e,b,m,kk) for kk in (['1','2'] if k=='mean_of_two' else [k])]
                    vv=[float(fmap[rid,f]['value']) for rid in rr];v=float(np.mean(vv));vals.append(v);ids.append('|'.join(rr));u=fmap[rr[0],f]['unit']
                    if k=='mean_of_two':cells.append(dict(feature=f,expert=e,block=b,material=m,value=v,unit=u,record_ids='|'.join(rr),repeat_values=';'.join(str(x) for x in vv),aggregation='arithmetic mean of the two record features; equal record weights'))
                category,signs,tol=sequence(vals)
                mono.append(dict(feature=f,expert=e,block=b,repeat=k,unit=u,material_order='P1<P3<P2<P4',P1=vals[0],P3=vals[1],P2=vals[2],P4=vals[3],adjacent_signs=';'.join(str(int(x)) for x in signs),monotonicity=category,numerical_sign_tolerance=tol,record_ids=';'.join(ids)))
                for i in range(3):step.append(dict(feature=f,expert=e,block=b,repeat=k,material_from=ORDER[i],material_to=ORDER[i+1],delta=vals[i+1]-vals[i],unit=u,direction=int(signs[i]),record_ids_from=ids[i],record_ids_to=ids[i+1]))
                endpoints.append(dict(feature=f,expert=e,block=b,repeat=k,delta_P4_minus_P1=vals[-1]-vals[0],unit=u,P1_record_ids=ids[0],P4_record_ids=ids[-1]))
    for n,rows in [('material_cell_means.csv',cells),('material_adjacent_changes.csv',step),('material_monotonicity.csv',mono),('material_endpoint_changes.csv',endpoints)]:write(n,rows)
    # Time curves are actual third-level means, not representative single cycles.
    thirds=defaultdict(dict);tcount={}
    for q in read(INPUTS[3]):
        key=(q['record_id'],['early','middle','late'][int(q['third'])-1],q['channel']);p=float(q['phase']);assert p not in thirds[key]
        thirds[key][p]=float(q['mean_deg_s']);tcount[key,p]=int(q['n_cycles'])
    assert len(thirds)==32*3*3
    contrasts=[('early','middle'),('middle','late'),('early','late')]
    temporal=[];repeat_temporal=[];tfeatures=[]
    for rid in sorted(rec):
        for seg,ch in itertools.product(['early','middle','late'],['wx','wy','wz']):
            d=decompose(thirds[rid,seg,ch],thirds[rid,seg,ch]);tfeatures.append(dict(record_id=rid,**meta(rid),segment=seg,channel=ch,mean=d['mean_a'],amplitude=d['amplitude_a'],unit='deg/s',n_common=d['n_common'],n_cycles=min(tcount[(rid,seg,ch),p] for p in thirds[rid,seg,ch]),phase_ids=d['phase_ids']))
        for (sa,sb),ch in itertools.product(contrasts,['wx','wy','wz']):
            temporal.append(dict(record_id=rid,**meta(rid),segment_a=sa,segment_b=sb,contrast=sb+'-'+sa,channel=ch,unit='deg/s',squared_unit='(deg/s)^2',n_cycles_a=min(tcount[(rid,sa,ch),p] for p in thirds[rid,sa,ch]),n_cycles_b=min(tcount[(rid,sb,ch),p] for p in thirds[rid,sb,ch]),**decompose(thirds[rid,sa,ch],thirds[rid,sb,ch])))
    for q in [q for q in pairs if q['comparison']=='repeat' and q['channel'] in ['wx','wy','wz']]:
        a,b,ch=q['record_a'],q['record_b'],q['channel']
        for seg in ['early','middle','late']:
            repeat_temporal.append(dict(pair_id=a+'|'+b,record_a=a,record_b=b,**meta(a),segment=seg,channel=ch,unit='deg/s',squared_unit='(deg/s)^2',n_cycles_a=min(tcount[(a,seg,ch),p] for p in thirds[a,seg,ch]),n_cycles_b=min(tcount[(b,seg,ch),p] for p in thirds[b,seg,ch]),**decompose(thirds[a,seg,ch],thirds[b,seg,ch])))
    write('temporal_decomposition.csv',temporal);write('temporal_repeat_decomposition.csv',repeat_temporal);write('temporal_record_features.csv',tfeatures)
    temporal_summ=[]
    for contrast,ch in itertools.product(['middle-early','late-middle','late-early'],['wx','wy','wz']):
        ss=[q for q in temporal if q['contrast']==contrast and q['channel']==ch];out=dict(contrast=contrast,channel=ch,mean_amplitude_unit='deg/s',squared_unit='(deg/s)^2',shape_unit='dimensionless',n_records=len(ss))
        for k in ['delta_mean_b_minus_a','delta_amplitude_b_minus_a','rho','shape_RMS_dimensionless','center_fraction','amplitude_fraction','shape_fraction']:
            vals=[q[k] for q in ss if q[k]!=''];out[k+'_median']=float(np.median(vals));out[k+'_min']=float(min(vals));out[k+'_max']=float(max(vals))
        out['amplitude_increase_count']=sum(q['delta_amplitude_b_minus_a']>0 for q in ss);out['amplitude_decrease_count']=sum(q['delta_amplitude_b_minus_a']<0 for q in ss)
        rr=[q for q in repeat_temporal if q['channel']==ch];out['same_third_repeat_shape_distance_median']=float(np.median([q['shape_RMS_dimensionless'] for q in rr]));temporal_summ.append(out)
    write('temporal_summary.csv',temporal_summ)
    assert len(temporal)==288 and len(repeat_temporal)==144
    assert all(q['n_common']==101 and q['status']=='ok' for q in pairs+temporal+repeat_temporal)
if __name__=='__main__':main()
