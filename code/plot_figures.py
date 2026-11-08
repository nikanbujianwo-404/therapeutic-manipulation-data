"""Numerical panels for Fig. 1c and Figs. 2–6 from unchanged measured/frozen tables."""
from pathlib import Path
import csv,json,hashlib,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
import os
R=Path(os.environ['TBME_RUN_ROOT']);F=R/'figures';Q=R/'evidence/quantitative'
def read(p):
 with (R/p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def export(n,rows):
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (R/'figure_sources'/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def save(fig,n):
 for ext in ['pdf','svg','png']:fig.savefig(F/(n+'.'+ext),dpi=300,bbox_inches='tight',pad_inches=.06)
 plt.close(fig)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
order=['P1','P3','P2','P4'];colors=['#2c6cad','#2f8755','#ca633a','#864da5']
rids=['ExpertA_P1_B1_R1','ExpertA_P1_B1_R2','ExpertB_P1_B1_R1','ExpertB_P1_B1_R2'];ec=['#246ba2','#bd562d'];styles=['-','--']
handles=[Line2D([],[],color=ec[i//2],ls=styles[i%2],lw=1.5) for i in range(4)];labels=['ExpertA E1','ExpertA E2','ExpertB E1','ExpertB E2']
ONLY={int(v) for v in sys.argv[1:]} or set(range(1,7))
if 1 in ONLY:
 fig,ax=plt.subplots(figsize=(4.2,2.5))
 data=read('figure_sources/mechanical_curves.csv')
 for mat,color in zip(order,colors):
  for k in ['1','2','3']:
   q=[v for v in data if v['specimen']==mat and v['repeat']==k]
   ax.plot([float(v['relative_depth_mm']) for v in q],[float(v['force_N']) for v in q],color=color,lw=.9,alpha=.8,label=mat if k=='1' else None)
 ax.set(title='(c) Independent indentation',xlabel='Indentation from preload (mm)',ylabel='Compression (N)');ax.spines[['top','right']].set_visible(False);ax.legend(ncol=2,frameon=False,loc='upper left');ax.grid(alpha=.12)
 fig.subplots_adjust(left=.02,right=.985,bottom=.25,top=.86,wspace=.35);save(fig,'Fig1c_indentation')


if 2 in ONLY:
 native=read('figure_sources/P1_native_cycle_state_sequences.csv');cloud=read('figure_sources/P1_actual_native_state_display.csv');profiles=read('evidence/formal_ati_imu_phase_profiles.csv')
 fig=plt.figure(figsize=(7.2,5.1));gs=fig.add_gridspec(2,1,height_ratios=[1.35,1],hspace=.53);topgs=gs[0].subgridspec(1,3,wspace=.25);botgs=gs[1].subgridspec(1,4,wspace=.85)
 paths=[]
 for j,(group,chs,unit) in enumerate([('force',['Fx','Fy','Fz'],'N'),('moment',['Tx','Ty','Tz'],'N m'),('gyro',['wx','wy','wz'],'deg/s')]):
  ax=fig.add_subplot(topgs[0,j],projection='3d',computed_zorder=False)
  q=[v for v in cloud if v['group']==group and v['record_id'] in rids];xx=np.array([[float(v[c]) for c in chs] for v in q]);ax.scatter(*xx.T,s=.55,c='#aaaaaa',alpha=.18,depthshade=False,rasterized=True,zorder=1)
  for k,rid in enumerate(rids):
   q=[v for v in native if v['group']==group and v['record_id']==rid and v['segment']=='middle'];q.sort(key=lambda v:float(v['host_time_s']))
   assert q;yy=np.array([[float(v[c]) for c in chs] for v in q]);tt=np.array([float(v['host_time_s']) for v in q]);cuts=np.r_[0,np.flatnonzero(np.diff(tt)>.02)+1,len(q)]
   for a,b in zip(cuts[:-1],cuts[1:]):
    ax.plot(*yy[a:b].T,color=ec[k//2],ls=styles[k%2],lw=1.5,zorder=3)
    for frac in [.25,.60]:
     i=a+int((b-a-1)*frac)
     if i+2<b:ax.quiver(*yy[i],*(1.6*(yy[i+2]-yy[i])),color='#303030',arrow_length_ratio=.55,linewidth=.9,zorder=4)
   paths.append(dict(panel=chr(97+j),record_id=rid,group=group,segment='middle',start_host_s=q[0]['cycle_start_host_s'],end_host_s=q[0]['cycle_end_host_s'],native_samples=len(q),continuous_runs=len(cuts)-1,background_rows=sum(v['group']==group and v['record_id']==rid for v in cloud)))
  ax.set(xlabel=chs[0],ylabel=chs[1],zlabel=chs[2],title=f'({chr(97+j)}) '+('Angular velocity' if group=='gyro' else group.capitalize())+' ('+unit+')');ax.view_init(elev=21,azim=-58);ax.tick_params(labelsize=8,pad=0);ax.xaxis.labelpad=1;ax.yaxis.labelpad=1;ax.zaxis.labelpad=1;ax.locator_params(nbins=3)
 for j,ch in enumerate(['Fz','Fy','Ty','wx']):
  ax=fig.add_subplot(botgs[0,j]);unit='N' if ch.startswith('F') else ('N m' if ch.startswith('T') else 'deg/s')
  for k,rid in enumerate(rids):
   q=sorted([v for v in profiles if v['record_id']==rid and v['channel']==ch],key=lambda v:float(v['phase']));assert len(q)==101
   ax.plot([float(v['phase']) for v in q],[float(v['mean']) for v in q],color=ec[k//2],ls=styles[k%2],lw=1.35)
  ax.set(title=f'({chr(100+j)}) {ch}',xlabel='IMU-event progress',ylabel=f'{ch} ({unit})');ax.set_xticks([0,.5,1]);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
 fig.subplots_adjust(left=.065,right=.965,bottom=.12,top=.86);fig.legend(handles,labels,ncol=4,loc='upper center',frameon=False,bbox_to_anchor=(.5,1.01),columnspacing=1.0,handlelength=2)
 save(fig,'Fig2_states_signed_cycles');export('Fig2_native_path_IDs.csv',paths)
 export('Fig2_mean_cycle_values.csv',[v for v in profiles if v['record_id'] in rids and v['channel'] in ['Fz','Fy','Ty','wx']])


if 3 in ONLY:
 opt=read('figure_sources/formal_optical_cycle_profiles.csv');cy=read('evidence/optical_formal_cycle_pairwise_distances.csv');task=read('evidence/optical_formal_pairwise_distances.csv');meta={q['record_id']:q for q in read('evidence/formal_ati_imu_records.csv')}
 def rk(rid):
  q=meta[rid];return(q['expert'],int(q['block']),order.index(q['material']),int(q['repeat']))
 def pk(q):return(rk(q['record_a']),rk(q['record_b']))
 pp=[]
 for comp in ['repeat','expert','material']:
  for q in sorted([v for v in cy if v['comparison']==comp],key=pk):
   t=next(v for v in task if v['comparison']==comp and v['record_a']==q['record_a'] and v['record_b']==q['record_b'])
   pp.append(dict(pair_index=len(pp)+1,pair_id=comp+'|'+q['record_a']+'|'+q['record_b'],comparison=comp,record_a=q['record_a'],record_b=q['record_b'],cycles_a=q['cycles_a'],cycles_b=q['cycles_b'],task_common_bins=t['common_bins'],cycle_position_mm=q['cycle_mean_position_shape_RMS_mm'],cycle_orientation_deg=q['cycle_mean_SO3_shape_RMS_deg'],task_position_mm=t['position_shape_RMS_mm'],task_orientation_deg=t['SO3_RMS_deg']))
 assert len(pp)==80
 fig=plt.figure(figsize=(7.2,5.1));gs=fig.add_gridspec(2,1,height_ratios=[1.65,1],hspace=.55);topgs=gs[0].subgridspec(1,2,wspace=.35);botgs=gs[1].subgridspec(1,4,wspace=.85)
 ax=fig.add_subplot(topgs[0,0],projection='3d');rgs=topgs[0,1].subgridspec(3,1,hspace=.30);rots=[fig.add_subplot(rgs[j,0]) for j in range(3)];rotrows=[]
 for k,rid in enumerate(rids):
  q=sorted([v for v in opt if v['record_id']==rid],key=lambda v:float(v['phase']));assert len(q)==101
  xyz=np.array([[float(v[c+'_mm_mean']) for c in ['dx','dy','dz']] for v in q]);ax.plot(*xyz.T,color=ec[k//2],ls=styles[k%2],lw=1.5)
  quat=np.array([[float(v['cycle_mean_q_'+c]) for c in ['x','y','z','w']] for v in q]);quat/=np.linalg.norm(quat,axis=1)[:,None];ang=np.degrees(2*np.arccos(np.clip(abs(quat[:,3]),0,1)))
  quat[quat[:,3]<0]*=-1
  si=np.linalg.norm(quat[:,:3],axis=1);theta=2*np.arctan2(si,quat[:,3]);vv=np.zeros((len(q),3));good=si>1e-12;vv[good]=quat[good,:3]/si[good,None]*theta[good,None]
  recon=np.column_stack([np.zeros_like(vv),np.ones(len(q))]);recon[good,:3]=vv[good]/theta[good,None]*np.sin(theta[good,None]/2);recon[good,3]=np.cos(theta[good]/2)
  err=np.max(np.linalg.norm(recon-quat,axis=1));assert err<1e-10 and np.max(theta)<np.pi
  for j,rr in enumerate(rots):rr.plot([float(v['phase']) for v in q],np.degrees(vv[:,j]),color=ec[k//2],ls=styles[k%2],lw=1.3)
  rotrows.extend(dict(record_id=rid,phase=v['phase'],log_mean_cycle_x_deg=z[0],log_mean_cycle_y_deg=z[1],log_mean_cycle_z_deg=z[2],reference='cycle-start moving-frame axes',quaternion_reconstruction_max_error=err) for v,z in zip(q,np.degrees(vv)))
 ax.set(title='(a) Mean spatial cycles',xlabel='dx cycle (mm)',ylabel='dy cycle (mm)',zlabel='dz cycle (mm)');ax.view_init(elev=22,azim=-58);ax.tick_params(labelsize=8,pad=0);ax.xaxis.labelpad=1;ax.yaxis.labelpad=1;ax.zaxis.labelpad=1;ax.locator_params(nbins=3)
 for j,rr in enumerate(rots):
  rr.set_ylabel(r'$\lambda_'+['x','y','z'][j]+r'$ (deg)');rr.set_xticks([0,.5,1]);rr.spines[['top','right']].set_visible(False);rr.grid(alpha=.12);rr.text(.02,.82,'('+chr(98+j)+')',transform=rr.transAxes,fontsize=9)
  if j<2:rr.tick_params(axis='x',labelbottom=False)
  else:rr.set_xlabel('Optical-cycle progress')
 for j,(key,title,unit) in enumerate([('cycle_position_mm','Cycle position','mm'),('cycle_orientation_deg','Cycle orientation','deg'),('task_position_mm','Task position','mm'),('task_orientation_deg','Task orientation','deg')]):
  ax=fig.add_subplot(botgs[0,j]);
  for comp,color in [('repeat','#687f91'),('expert','#b34f32'),('material','#368fa2')]:
   q=[v for v in pp if v['comparison']==comp];ax.scatter([v['pair_index'] for v in q],[float(v[key]) for v in q],s=8,color=color,alpha=.8,edgecolors='none')
  base=np.median([float(v[key]) for v in pp if v['comparison']=='repeat']);ax.axhline(base,color='black',ls='--',lw=.7);ax.set(title=f'({chr(101+j)}) '+title,ylabel='RMS distance ('+unit+')');ax.set_xticks([1,16,32,80]);ax.set_xlabel('Pair index');ax.tick_params(axis='x',length=2);ax.axvline(16.5,color='#bbbbbb',lw=.5);ax.axvline(32.5,color='#bbbbbb',lw=.5);ax.set_xlim(0,81);ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False)
 fig.subplots_adjust(left=.065,right=.99,bottom=.13,top=.86);fig.legend(handles,labels,ncol=4,loc='upper center',frameon=False,bbox_to_anchor=(.5,1.01),columnspacing=1.0)
 save(fig,'Fig3_space_physical_distances');export('Fig3_mean_cycle_rotation_logs.csv',rotrows);export('Fig3_all_spatial_pair_IDs.csv',pp);export('Fig3_spatial_cycle_values.csv',[v for v in opt if v['record_id'] in rids])

if 4 in ONLY:
 P=read('evidence/quantitative/pair_decomposition.csv');chans=['Fx','Fy','Fz','Tx','Ty','Tz','wx','wy','wz'];cs=['#617b96','#d49536','#a84a55']
 fig=plt.figure(figsize=(7.2,2.95));gs=fig.add_gridspec(1,3,width_ratios=[2.0,1,1]);ax=fig.add_subplot(gs[0]);xx=np.arange(9)
 for offset,comp,hatch in [(-.18,'expert',None),(.18,'repeat','//')]:
     bottom=np.zeros(9)
     for key,c in zip(['center','amplitude','shape'],cs):
         vals=np.array([np.mean([float(q[key+'_fraction']) for q in P if q['comparison']==comp and q['channel']==ch])*100 for ch in chans])
         ax.bar(xx+offset,vals,width=.34,bottom=bottom,color=c,hatch=hatch,edgecolor='white',linewidth=.25);bottom+=vals
     assert np.max(abs(bottom-100))<1e-10
 dominance=[]
 for i,ch in enumerate(chans):
  ss=[v for v in P if v['comparison']=='expert' and v['channel']==ch];counts={k:0 for k in ['center','amplitude','shape']};ties=0
  for v in ss:
   values=np.array([float(v[k+'_MSE']) for k in counts]);winner=np.flatnonzero(values==values.max())
   if len(winner)>1:ties+=1
   else:counts[list(counts)[int(winner[0])]]+=1
  top=max(counts,key=counts.get)
  dominance.append(dict(channel=ch,n_practitioner_pairs=len(ss),dominant_center=counts['center'],dominant_amplitude=counts['amplitude'],dominant_shape=counts['shape'],exact_ties=ties,displayed_most_frequent_term=top,rule='strict unique maximum of the three unrounded MSE terms; exact ties separate'))
 export('Fig4_dominance_counts.csv',dominance)
 dmap={v['channel']:v for v in dominance}
 count_text=[]
 for key,label in [('center','C'),('amplitude','A'),('shape','S')]:
  groups={}
  for ch in chans:
   n=dmap[ch]['dominant_'+key]
   if n:groups.setdefault(n,[]).append(ch)
  count_text.append(label+': '+'; '.join('/'.join(groups[n])+' '+str(n)+'/16' for n in sorted(groups,reverse=True)))
 for text,y in zip(count_text,[.98,.895,.81]):ax.text(.01,y,text,transform=ax.transAxes,ha='left',va='top',fontsize=8)
 ax.set_title('(a) Mean shares and practitioner counts',loc='left');ax.set_ylabel('Mean pair contribution (%)');ax.set_xticks(xx,chans,fontsize=8);ax.set_ylim(0,151);ax.set_yticks([0,25,50,75,100]);ax.spines[['top','right']].set_visible(False)
 for j,ch in enumerate(['Fz','wx'],1):
     ax=fig.add_subplot(gs[j]);ss=sorted([q for q in P if q['comparison']=='expert' and q['channel']==ch],key=lambda q:(order.index(q['material_a']),q['block_a'],q['repeat_a']));bottom=np.zeros(16)
     for key,c in zip(['center','amplitude','shape'],cs):
         vals=np.array([float(q[key+'_MSE']) for q in ss]);ax.bar(range(16),vals,bottom=bottom,width=.85,color=c);bottom+=vals
     base=np.median([float(q['distance_squared']) for q in P if q['comparison']=='repeat' and q['channel']==ch]);ax.axhline(base,color='black',ls='--',lw=.8)
     ax.set_title('('+['','b','c'][j]+') '+ch+' expert pairs',loc='left');ax.set_ylabel('Distance squared ('+('N$^2$' if ch=='Fz' else '(deg/s)$^2$')+')',fontsize=8);ax.set_xticks([0,7,15],['1','8','16']);ax.spines[['top','right']].set_visible(False)
 fig.subplots_adjust(left=.075,right=.99,bottom=.26,top=.88,wspace=.6)
 fig.legend([Patch(color=c) for c in cs]+[Patch(facecolor='white',edgecolor='gray',hatch='//'),Line2D([],[],color='black',ls='--')],['Center','Amplitude','Shape','Repeat (hatched)','Repeat median in (b,c)'],ncol=3,loc='lower center',frameon=False,fontsize=8)
 save(fig,'Fig4_distance_composition')

if 5 in ONLY:
 features=['cycle_compression_center','cycle_amplitude_wx','position_x_range_p05_p95_mm','orientation_excursion_p95_deg'];names=['Compression center\n(N)','wx cycle amplitude\n(deg/s)','X activity range\n(mm)','Orientation excursion\nP95 (deg)'];ff=read('evidence/quantitative/record_features.csv');cells=read('evidence/quantitative/material_cell_means.csv')
 fig,axs=plt.subplots(2,4,figsize=(7.2,3.7))
 for row,e in enumerate(['ExpertA','ExpertB']):
     for col,(feat,name) in enumerate(zip(features,names)):
         ax=axs[row,col]
         for b,c in [('1','#246ba2'),('2','#bd562d')]:
             for k,marker in [('1','o'),('2','s')]:
                 vals=[next(float(q['value']) for q in ff if (q['feature'],q['expert'],q['block'],q['repeat'],q['material'])==(feat,e,b,k,m)) for m in order]
                 ax.plot(range(4),vals,color=c,ls='--',lw=.65,alpha=.42,marker=marker,ms=3,mfc='white')
             vals=[next(float(q['value']) for q in cells if (q['feature'],q['expert'],q['block'],q['material'])==(feat,e,b,m)) for m in order]
             ax.plot(range(4),vals,color=c,lw=1.5,marker='.',ms=4)
         ax.set_title(('ExpertA' if e=='ExpertA' else 'ExpertB')+' / '+name,fontsize=8,loc='left');ax.set_xticks(range(4),order);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.1)
 fig.subplots_adjust(left=.065,right=.99,bottom=.2,top=.88,wspace=.48,hspace=.75)
 fig.legend([Line2D([],[],color=c) for c in ['#246ba2','#bd562d']]+[Line2D([],[],color='gray',ls='--',marker=m,mfc='white') for m in ['o','s']]+[Line2D([],[],color='gray',lw=1.5)],['Round 1','Round 2','Execution 1','Execution 2','Mean of two'],ncol=5,loc='lower center',frameon=False,fontsize=8)
 save(fig,'Fig5_material_responses')


if 6 in ONLY:
 tt=read('evidence/quantitative/temporal_decomposition.csv');thirds=read('evidence/temporal/formal_third_imu_profiles.csv')
 fig=plt.figure(figsize=(7.2,4.65));gs=fig.add_gridspec(2,1,height_ratios=[1,1.05],hspace=.85);topgs=gs[0].subgridspec(1,2,wspace=.45);botgs=gs[1].subgridspec(1,3,wspace=.55);pcols=['#1769b0','#228b45','#d62728'];wave=[]
 for j,rid in enumerate([rids[0],rids[2]]):
  ax=fig.add_subplot(topgs[0,j]);cs=[]
  for third,col in zip(['1','2','3'],pcols):
   q=sorted([v for v in thirds if v['record_id']==rid and v['channel']=='wx' and v['third']==third],key=lambda v:float(v['phase']));assert len(q)==101;cs.append(q[0]['n_cycles']);wave+=q
   ax.plot([float(v['phase']) for v in q],[float(v['mean_deg_s']) for v in q],color=col,lw=1.5)
  ax.set(title=f'({chr(97+j)}) '+('ExpertA' if j==0 else 'ExpertB')+' / Selected execution',xlabel='IMU-event progress',ylabel='Mean wx (deg/s)');ax.set_xticks([0,.5,1]);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
 for j,ch in enumerate(['wx','wy','wz']):
  ax=fig.add_subplot(botgs[0,j])
  for e,col in [('ExpertA','#246ba2'),('ExpertB','#bd562d')]:
   for b,mark in [('1','o'),('2','s')]:
    q=[v for v in tt if v['contrast']=='late-early' and v['channel']==ch and v['expert']==e and v['block']==b];assert len(q)==8
    ax.scatter([float(v['delta_amplitude_b_minus_a']) for v in q],[float(v['shape_RMS_dimensionless']) for v in q],s=17,color=col,marker=mark,alpha=.8,edgecolors='white',linewidth=.25)
  ax.set(title=f'({chr(99+j)}) '+ch,xlabel='Amplitude change (deg/s)',ylabel='Shape RMS (1)');ax.axvline(0,color='gray',lw=.6);ax.set_ylim(bottom=0);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.12)
 fig.subplots_adjust(left=.075,right=.99,bottom=.20,top=.87)
 fig.legend([Line2D([],[],color=c,lw=1.5) for c in pcols],['Early period mean','Middle period mean','Late period mean'],ncol=3,loc='upper center',bbox_to_anchor=(.5,1.01),frameon=False)
 fig.legend([Line2D([],[],color=c,marker='o',ls='') for c in ec]+[Line2D([],[],color='gray',marker=m,ls='') for m in ['o','s']],['ExpertA','ExpertB','Round 1','Round 2'],ncol=4,loc='lower center',frameon=False)
 save(fig,'Fig6_temporal_waveforms_components');export('Fig6_period_mean_examples.csv',wave)
 export('Fig6_late_early_all_components.csv',[v for v in tt if v['contrast']=='late-early'])
 # All three contrasts are tabulated, without creating independent-cycle inference.
 summary=[]
 for contrast in ['middle-early','late-middle','late-early']:
  for ch in ['wx','wy','wz']:
   q=[v for v in tt if v['contrast']==contrast and v['channel']==ch];assert len(q)==32
   summary.append(dict(contrast=contrast,channel=ch,unit='deg/s',median_center_change=np.median([float(v['delta_mean_b_minus_a']) for v in q]),median_amplitude_change=np.median([float(v['delta_amplitude_b_minus_a']) for v in q]),median_shape_RMS=np.median([float(v['shape_RMS_dimensionless']) for v in q]),amplitude_increases=sum(float(v['delta_amplitude_b_minus_a'])>0 for v in q),amplitude_decreases=sum(float(v['delta_amplitude_b_minus_a'])<0 for v in q)))
 export('Fig6_all_contrast_summary.csv',summary)
(R/'verification/FIGURE_DATA_CHECK.json').write_text(json.dumps(dict(figures=6,updated_figures=sorted(ONLY),case_rule='P1, round 1, both repeats/both practitioners; chronological middle central valid event for native state paths; P1/round1/repeat1 for period-mean examples',native_state_path_count=12,spatial_pairs=80,all_records_material_time=32,unchanged_queries_masks_coordinates=True,normal_fonts_source_pt=8),indent=2))
print('Six measured-data figures and source identity manifests generated.')
