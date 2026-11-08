#!/usr/bin/env python3
"""Validate extracted release files without changing them."""
from pathlib import Path
import argparse,csv,gzip,hashlib,itertools,json
import numpy as np
import pandas as pd
from optical_geometry_helpers import read_pose_csv

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path,default=Path('data'));p.add_argument('--report',type=Path);args=p.parse_args()
    data=args.data.resolve();repo=Path(__file__).resolve().parents[1]
    manifest=pd.read_csv(repo/'record_manifest.csv');assert len(manifest)==32 and manifest.record_id.nunique()==32
    expected=set(itertools.product(['ExpertA','ExpertB'],['P1','P2','P3','P4'],[1,2],[1,2]))
    assert set(zip(manifest.practitioner,manifest.specimen,manifest['round'],manifest.repetition))==expected
    paths=[c for c in manifest if c.endswith('_path')]
    for r in manifest.to_dict('records'):
        for k in paths:assert (data.parent/r[k]).exists(),(r['record_id'],k)
    with (repo/'FILE_MANIFEST.csv').open(newline='') as f:files=list(csv.DictReader(f))
    for r in files:
        path=data.parent/r['path'];assert path.is_file() and not path.is_symlink(),path
        assert path.stat().st_size==int(r['size_bytes']) and sha(path)==r['sha256'],path
    npz_count=csv_count=json_count=0;arrays=0
    for path in data.rglob('*'):
        if not path.is_file():continue
        assert not path.is_symlink()
        if path.suffix=='.npz':
            with np.load(path,allow_pickle=False) as z:
                for k in z.files:a=z[k];assert a.dtype.kind!='O';arrays+=1
            npz_count+=1
        elif path.suffix=='.json':json.loads(path.read_text());json_count+=1
        elif path.suffix=='.csv' or path.name.endswith('.csv.gz'):
            op=gzip.open if path.suffix=='.gz' else open
            with op(path,'rt',newline='',encoding='utf-8-sig') as f:
                reader=csv.reader(f);header=next(reader);n=sum(1 for _ in reader);assert n>0
            csv_count+=1
    records=[]
    for r in manifest.itertuples():
        rid=r.record_id;proc=data/'processed/records'/rid
        force=pd.read_csv(data.parent/r.ati_path);imu=pd.read_csv(data.parent/r.imu_path)
        rawdir=data/'raw/records'/rid
        assert (rawdir/'ati/raw.rdt').stat().st_size==36*len(force)
        assert (rawdir/'imu/raw.bin').stat().st_size==28*len(imu)
        for modality,ledger,binary in [('ati','packets.csv','raw.rdt'),('imu','blocks.csv','raw.bin')]:
            blocks=pd.read_csv(rawdir/modality/ledger)
            assert int(blocks.offset.iloc[0])==0
            assert np.array_equal((blocks.offset+blocks.length).to_numpy()[:-1],blocks.offset.to_numpy()[1:])
            assert int((blocks.offset+blocks.length).iloc[-1])==(rawdir/modality/binary).stat().st_size
        with np.load(proc/'ati_imu_native.npz') as z:
            ix=z['force_source_row'].astype(int);iy=z['imu_source_row'].astype(int)
            wrench=force[['zeroed_fx_N','zeroed_fy_N','zeroed_fz_N','zeroed_tx_Nm','zeroed_ty_Nm','zeroed_tz_Nm']].to_numpy()[ix]
            assert z['force_sensor_zeroed'].shape==(len(ix),6)
            assert np.allclose(wrench,z['force_sensor_zeroed'],rtol=0,atol=1e-12,equal_nan=True)
            assert np.allclose(-wrench[:,2],z['force_down_source_N'],rtol=0,atol=1e-12,equal_nan=True)
            assert np.allclose(imu[['wx_deg_s','wy_deg_s','wz_deg_s']].to_numpy()[iy],z['imu_gyro_deg_s'],rtol=0,atol=1e-12,equal_nan=True)
        config=json.loads((proc/'ati_calibration.json').read_text());bias=np.array(config['record_bias']);fields=config['sensor_fields']
        raw=force[['raw_fx_count','raw_fy_count','raw_fz_count','raw_tx_count','raw_ty_count','raw_tz_count']].to_numpy(float)
        si=force[['fx_N','fy_N','fz_N','tx_Nm','ty_Nm','tz_Nm']].to_numpy()
        assert np.allclose(raw/np.array([float(fields['cfgcpf'])]*3+[float(fields['cfgcpt'])]*3),si,atol=1e-9,equal_nan=True)
        assert np.allclose(si-bias,force[['zeroed_fx_N','zeroed_fy_N','zeroed_fz_N','zeroed_tx_Nm','zeroed_ty_Nm','zeroed_tz_Nm']].to_numpy(),atol=1e-9,equal_nan=True)
        optical=[]
        for tracker in ['Tracker0','Tracker1']:
            fr,ts,pos,q=read_pose_csv(data.parent/r.nokov_path/(tracker+'.csv'));assert len(fr)==len(ts)==len(pos)==len(q);assert pos.shape[1]==3 and q.shape[1]==4;optical.append(len(fr))
        records.append(dict(record_id=rid,ati_rows=len(force),imu_rows=len(imu),optical_frames=optical,ati_output_rate_hz=float(fields['comrdtrate']),six_axis_verified=True))
        print('Validated',rid,flush=True)
    clocks=json.loads((data/'processed/config/optical_clock_models.json').read_text());assert set(clocks)==set(manifest.record_id)
    assert all(c['physical_sample_latency_verified'] is False and c['whole_record_motion_validation_passed'] is False for c in clocks.values())
    align=pd.read_csv(data/'processed/config/ati_imu_alignment.csv');assert len(align)==32 and not align.sample_epoch_phase_calibrated.any()
    material=pd.read_csv(repo/'material_manifest.csv');assert len(material)==12 and material.trial_id.nunique()==12
    for k in ['replay_path','material_path']:
        assert all((data.parent/v).is_file() for v in material[k])
    result=dict(status='PASS',main_records=32,material_trials=12,manifest_files_verified=len(files),readable_npz_files=npz_count,readable_arrays=arrays,readable_csv_files=csv_count,readable_json_files=json_count,records=records,physical_synchronization_verified=False)
    if args.report:
        if args.report.exists():raise SystemExit('Refusing to overwrite a report')
        args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()
