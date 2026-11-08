#!/usr/bin/env python3
"""Portable release entry. Reuse frozen inputs and original scientific kernels."""
from pathlib import Path
import argparse,csv,json,os,shutil,subprocess,sys,itertools
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
import optical_core as optical
from phantom_force_signal import measured_force_at

def compare_numeric(generated, saved):
    a=pd.read_csv(generated);b=pd.read_csv(saved)
    assert len(a)==len(b),(generated,len(a),len(b))
    # Public record identifiers sort differently from the internal acquisition IDs.
    # Compare rows by their scientific identity, not the old filename order.
    candidates=['record_id','comparison','record_a','record_b','channel','feature','expert','block','repeat','material','material_from','material_to','segment','segment_a','segment_b','contrast','cycle_index','start_frame','phase']
    keys=[k for k in candidates if k in a and k in b]
    if keys:
        a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True)
    errors={}
    for k in a.columns.intersection(b.columns):
        if pd.api.types.is_numeric_dtype(a[k]) and pd.api.types.is_numeric_dtype(b[k]):
            av=a[k].to_numpy(float);bv=b[k].to_numpy(float)
            assert np.allclose(av,bv,rtol=1e-9,atol=1e-8,equal_nan=True),(generated,k)
            errors[k]=float(np.nanmax(np.abs(av-bv))) if np.isfinite(av-bv).any() else 0.
        elif k not in ['source']:
            assert a[k].fillna('').astype(str).tolist()==b[k].fillna('').astype(str).tolist(),(generated,k)
    return dict(rows=len(a),numeric_columns=len(errors),max_absolute_difference=max(errors.values(),default=0.))

def optical_reproduction(data,out):
    bounds=pd.read_csv(data/'processed/config/task_boundaries.csv').set_index('record_id')
    meta=pd.read_csv(data/'processed/record_index.csv').set_index('record_id')
    dest=out/'optical';dest.mkdir();(dest/'cycles').mkdir();(dest/'task_profiles').mkdir()
    profiles=[];cycles={};tasks={};events=[];diagnostics=[];qc=[]
    for rid,b in bounds.iterrows():
        g=pd.read_csv(data/f'processed/records/{rid}/optical_qualified_poses.csv.gz')
        t=g.optical_s.to_numpy();m=(t>=b.start_optical_s)&(t<=b.end_optical_s)&g.valid.to_numpy(bool)
        v,w=optical.derivatives(g,m)
        ev,diag=optical.formal_native_events(g,m,w)
        ii=np.flatnonzero(m&(t<b.start_optical_s+1))
        ref=(g[optical.P].to_numpy()[ii].mean(0),Rotation.from_quat(g[optical.Q].to_numpy()[ii]).mean())
        cy,cq,pr=optical.cycle_profiles(g,b,m,*ref,v,w,ev)
        assert cy is not None
        np.savez_compressed(dest/'cycles'/f'{rid}.npz',**cy)
        tp=optical.profile(g,m,ref,b.start_optical_s,b.end_optical_s,w)
        np.savez_compressed(dest/'task_profiles'/f'{rid}.npz',**tp)
        cycles[rid]=cy;tasks[rid]=tp;profiles+=pr;events+=ev.to_dict('records');diagnostics.append(diag);qc+=cq.to_dict('records')
        print('Optical reconstruction',rid,'accepted',len(cy['times']),flush=True)
    pd.DataFrame(profiles).to_csv(dest/'cycle_mean_profiles.csv',index=False)
    pd.DataFrame(events).to_csv(dest/'cycle_candidates.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(dest/'axis_diagnostics.csv',index=False)
    pd.DataFrame(qc).to_csv(dest/'cycle_quality.csv',index=False)
    checks={}
    for new,old in [('cycle_mean_profiles.csv','processed/figure_sources/formal_optical_cycle_profiles.csv'),('cycle_candidates.csv','processed/optical/formal_native_cycle_candidates.csv'),('axis_diagnostics.csv','processed/optical/formal_cycle_event_axis_QC.csv'),('cycle_quality.csv','processed/optical/formal_cycle_QC.csv')]:
        checks[new]=compare_numeric(dest/new,data/old)
    # Apply the paper's pair identities to the rebuilt mean paths and SO(3) means.
    saved=pd.read_csv(data/'processed/evidence/optical_formal_cycle_pairwise_distances.csv')
    rows=[]
    for r in saved.itertuples():
        ca,cb=cycles[r.record_a],cycles[r.record_b]
        dp=ca['array'][:,:,:3].mean(0)-cb['array'][:,:,:3].mean(0)
        theta=(Rotation.from_quat(ca['meanq']).inv()*Rotation.from_quat(cb['meanq'])).magnitude()*180/np.pi
        rows.append(dict(comparison=r.comparison,record_a=r.record_a,record_b=r.record_b,cycles_a=len(ca['times']),cycles_b=len(cb['times']),cycle_mean_position_shape_RMS_mm=np.sqrt(np.mean(np.sum(dp**2,axis=1))),cycle_mean_SO3_shape_RMS_deg=np.sqrt(np.mean(theta**2))))
    pd.DataFrame(rows).to_csv(dest/'cycle_pair_distances.csv',index=False)
    checks['cycle_pair_distances']=compare_numeric(dest/'cycle_pair_distances.csv',data/'processed/evidence/optical_formal_cycle_pairwise_distances.csv')
    taskrows=optical.compare_pairs(tasks,meta,'primary3_4__formal_mean_1s')
    pd.DataFrame(taskrows).to_csv(dest/'task_pair_distances.csv',index=False)
    checks['task_pair_distances']=compare_numeric(dest/'task_pair_distances.csv',data/'processed/evidence/optical_formal_pairwise_distances.csv')
    return checks

def material_reproduction(data,out,repo):
    rows=[]
    for r in pd.read_csv(repo/'material_manifest.csv').itertuples():
        rr=pd.read_csv(data.parent/r.replay_path);rr=rr[rr.state=='TRACKING']
        mm=pd.read_csv(data.parent/r.material_path);mm=mm[(mm.replay_state=='TRACKING')&(mm.fit_input_valid==1)]
        with (data.parent/r.replay_path).open(newline='') as f:records=[q for q in csv.DictReader(f) if q['state']=='TRACKING']
        mt=mm.timestamp.to_numpy();tm=np.interp(mt,rr.timestamp,rr.force_time,left=np.nan,right=np.nan)
        force=measured_force_at(records,mt,r.max_sync_gap_sec);depth=mm.d_rel_m.to_numpy()*1000
        ok=np.isfinite(tm)&np.isfinite(force)&np.isfinite(depth);grid=np.arange(.25,57.9,.04)
        for t,ff,dd in zip(grid,np.interp(grid,tm[ok],force[ok]),np.interp(grid,tm[ok],depth[ok])):
            rows.append(dict(specimen=r.specimen,repeat=r.repetition,time_s=t,force_N=ff,relative_depth_mm=dd))
    dest=out/'material';dest.mkdir();pd.DataFrame(rows).to_csv(dest/'mechanical_curves.csv',index=False)
    # Reorder both tables by specimen/repeat/time: source file uses its historical run order.
    a=pd.DataFrame(rows);b=pd.read_csv(data/'processed/figure_sources/mechanical_curves.csv')
    keys=['specimen','repeat','time_s'];a=a.sort_values(keys).reset_index(drop=True);b=b.sort_values(keys).reset_index(drop=True)
    errors={}
    for k in ['time_s','force_N','relative_depth_mm']:
        errors[k]=float(np.max(np.abs(a[k].to_numpy()-b[k].to_numpy())))
        assert np.allclose(a[k],b[k],rtol=1e-7,atol=1e-6),(k,errors[k])
    return dict(rows=len(a),max_absolute_errors=errors,slopes='saved trial fits supplied; slope fitting not rerun')

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=Path('data'));p.add_argument('--output',type=Path,default=Path('reproduced'));p.add_argument('--skip-plots',action='store_true');args=p.parse_args()
    data=args.data.resolve();out=args.output.resolve();repo=Path(__file__).resolve().parents[1]
    if out.exists():raise SystemExit('Output must be a new directory, to protect prior results.')
    out.mkdir(parents=True);shutil.copytree(data/'processed/evidence',out/'evidence');shutil.copytree(data/'processed/figure_sources',out/'figure_sources')
    for name in ['figures','verification']:(out/name).mkdir()
    env=dict(os.environ,TBME_RUN_ROOT=str(out));code=Path(__file__).resolve().parent
    subprocess.run([sys.executable,'-s',str(code/'analyze_quantitative.py')],env=env,check=True)
    checks={}
    for f in sorted((data/'results/quantitative').glob('*.csv')):checks[f.name]=compare_numeric(out/'evidence/quantitative'/f.name,f)
    checks['optical']=optical_reproduction(data,out)
    checks['material_curves']=material_reproduction(data,out,repo)
    if not args.skip_plots:subprocess.run([sys.executable,'-s',str(code/'plot_figures.py')],env=env,check=True)
    s=pd.read_csv(out/'evidence/quantitative/decomposition_summary.csv');t=pd.read_csv(out/'evidence/quantitative/temporal_decomposition.csv');cy=pd.read_csv(out/'optical/cycle_pair_distances.csv')
    key={}
    for ch,term,expected in [('Fz','center',91.4),('wx','shape',92.2),('wy','shape',95.9)]:
        row=s[(s.comparison=='expert')&(s.channel==ch)].iloc[0];value=float(row[term+'_fraction_median']*100);count=int(row['largest_'+term+'_count'])
        assert round(value,1)==expected and count==16
        key[ch]=dict(median_squared_distance_contribution_percent=value,dominant_pairs=count,total_pairs=16)
    for comp,expected in [('expert',19.37),('repeat',4.34)]:
        value=float(cy.loc[cy.comparison==comp,'cycle_mean_position_shape_RMS_mm'].median());assert round(value,2)==expected;key[comp+'_spatial_cycle_median_mm']=value
    for e,expected in [('ExpertA',10),('ExpertB',2)]:
        z=t[(t.expert==e)&(t.channel=='wx')&(t.contrast=='late-early')];value=int((z.delta_amplitude_b_minus_a>0).sum());assert value==expected and len(z)==16;key[e+'_wx_amplitude_increases']=dict(increases=value,records=16)
    report=dict(status='PASS',comparison_checks=checks,key_results=key,starting_point='frozen recording-mean ATI/IMU waveforms and qualified native optical poses; material recorder columns',not_rerun=['raw optical marker reconstruction and primary tracking masks','raw ATI/IMU event detection and mean-waveform generation','clock estimation or physical synchronization validation','material slope fitting'])
    (out/'verification/REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(key,indent=2),flush=True)
if __name__=='__main__':main()
