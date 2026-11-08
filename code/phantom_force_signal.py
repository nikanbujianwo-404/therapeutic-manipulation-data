"""Shared offline force convention: replay controller measured positive compression."""
import numpy as np
FORCE_SOURCE = 'replay.controller_measured_z_positive_compression_v1'
def measured_force_at(rows, timestamps, max_gap=.02):
    key = 'controller_measured_z' if rows and 'controller_measured_z' in rows[0] else 'measured_vertical_down'
    def value(r,k):
        try:return float(r[k])
        except (KeyError,TypeError,ValueError):return float('nan')
    rt=np.array([value(r,'timestamp') for r in rows]);f=np.array([value(r,key) for r in rows])
    if len(rt)<2 or not np.isfinite(rt).all() or np.any(np.diff(rt)<=0):
        raise ValueError('控制器实测力时间轴缺失或不单调')
    t=np.asarray(timestamps);result=np.interp(t,rt,f,left=np.nan,right=np.nan)
    j=np.searchsorted(rt,t,side='left');j=np.clip(j,1,len(rt)-1)
    exact=np.isin(t,rt)
    result[(rt[j]-rt[j-1]>max_gap)&~exact]=np.nan
    return result
