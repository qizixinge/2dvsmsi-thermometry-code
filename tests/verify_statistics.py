"""Independent verification of MATLAB results using least-squares design matrices.

This checks computed values, not text copied from the manuscript. A rebuilt PFA
result is not required to match the pre-fix numeric baseline.

2026-11-05: nine time points (calibration/LOOCV on all nine; Bland-Altman on 3-9),
reference-invariance check, B0 read from the DICOM header, acquisition-1 probe temperatures recomputed
from the 1 Hz log, and a check that every macro in output/results/numbers*.tex equals
the value recomputed here / stored in statistics_results_final.mat.
"""
from pathlib import Path
import sys, os
import json
import re
import numpy as np
import h5py
from scipy.io import loadmat

base = Path(__file__).resolve().parents[1]
data = base/'data'
out = Path(os.environ.get('VSMSI_OUTPUT_DIR', str(base/'output')))/'results'
s = loadmat(out/'statistics_results_final.mat', simplify_cells=True)
checks = []
def close(name, actual, expected, atol=1e-9, kind='independent'):
    # kind = 'independent': the expected value is recomputed here from the inputs (log, DICOM bytes, raw ROI means).
    # kind = 'consistency': a value stored in one file is compared with a copy of the same value in another file.
    err = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
    assert np.allclose(actual, expected, atol=atol, rtol=1e-9), (name,err)
    checks.append({'check':name,'max_absolute_difference':err,'kind':kind})

def regression(x,y,common):
    n,k=x.shape
    # Stack by ROI, with explicit intercept columns. No within-ROI centering.
    intercept=np.repeat(np.eye(k),n,axis=0)
    if common:
        design=np.column_stack((intercept,x.T.ravel()))
    else:
        slopes=np.zeros((n*k,k))
        for j in range(k): slopes[j*n:(j+1)*n,j]=x[:,j]
        design=np.column_stack((intercept,slopes))
    beta=np.linalg.lstsq(design,y.T.ravel(),rcond=None)[0]
    resid=y.T.ravel()-design@beta
    return beta, float(resid@resid)

cfg=s['loocv_cfg']
r=loadmat(data/'ROI_and_tempr_resolved.mat',simplify_cells=True)
NTP=9
assert int(s['n_tp'])==NTP and list(np.atleast_1d(s['BA_TPS']))==list(range(3,10)) and int(s['TP_CMP'])==3, 'time-point configuration changed: expected 9 time points, BA_TPS=3:9, TP_CMP=3'
for k_ in ['temp1','temp3','temprp1']:
    assert np.asarray(r[k_]).shape==(NTP,4),(k_,np.asarray(r[k_]).shape)
for k_ in ['ROI_2DVSMSI_100_series','ROI_PRF_128_series']:
    for q_,row_ in enumerate(r[k_]):
        row_=np.asarray(row_); base_=np.asarray(r['ROI_2DVSMSI_100' if '2DVSMSI' in k_ else 'ROI_PRF_128'])[q_]
        assert row_.shape==(NTP,2) and np.all(row_==base_),(k_,q_)
checks.append({'check':'nine time points in probe temperatures and ROI series (ROI coordinates constant)','passed':True})
# acquisition-1 probe temperatures recomputed from the 1 Hz fluoroptic log (scanner clock - 51 s)
import csv, io
def hms(t):
    h,m,sec=t.split(':'); return int(h)*3600+int(m)*60+float(sec)
rows=[]
with open(data/'probe_log_20250425.csv','rb') as fh:
    raw=fh.read().decode('gbk',errors='replace').splitlines()
for ln in raw[1:]:
    c=ln.split(',')
    try: vals=[float(c[i]) for i in (2,3,4,5)]; ts=hms(c[1])
    except Exception: continue
    if all(np.isfinite(vals)) and any(v>0 for v in vals): rows.append([ts]+vals)
rows=np.array(rows)
def winmean(t0,t1):
    a=np.ceil(hms(t0)-51); b=np.floor(hms(t1)-51)
    m=(rows[:,0]>=a)&(rows[:,0]<=b); return rows[m,1:].mean(axis=0)
# Stated windows (scanner clock): GRE = DICOM AcquisitionTime + 25.4 s (SliceMeasurementDuration); PFA/Full = first to last
# ADC time stamp of the raw .dat file (twix).  All 9 acquisitions x {GRE, PFA, Full} = 27 windows are recomputed from the log.
gre_start=['17:03:10.285','17:13:43.26','17:20:26.28','17:26:39.25','17:34:52.31','17:42:56.28','17:49:18.32','17:56:44.30','18:05:09.25']
pfa_win=[('17:08:12.32','17:08:42.00'),('17:14:16.40','17:14:46.08'),('17:21:00.65','17:21:30.33'),('17:27:12.46','17:27:42.14'),('17:35:27.39','17:35:57.07'),
         ('17:43:29.43','17:43:59.11'),('17:49:53.43','17:50:23.11'),('17:57:20.43','17:57:50.11'),('18:05:42.38','18:06:12.06')]
full_win=[('17:09:55.46','17:10:45.14'),('17:16:03.43','17:16:53.11'),('17:22:51.43','17:23:41.11'),('17:29:02.43','17:29:52.11'),('17:37:13.43','17:38:03.11'),
          ('17:45:16.43','17:46:06.11'),('17:51:45.42','17:52:35.10'),('17:59:07.39','17:59:57.07'),('18:07:26.43','18:08:16.11')]
GRE_WINDOW_S=25.4
for k_,wins in (('temprp1',[(g_,'%.3f'%(hms(g_)+GRE_WINDOW_S)) for g_ in gre_start]),('temp1',pfa_win),('temp3',full_win)):
    for t_,(t0_,t1_) in enumerate(wins):
        t1n=hms(t1_) if isinstance(t1_,str) and ':' in t1_ else float(t1_)
        a_=np.ceil(hms(t0_)-51); b_=np.floor(t1n-51)
        m_=(rows[:,0]>=a_)&(rows[:,0]<=b_); assert m_.sum()>=20,(k_,t_,int(m_.sum()))
        close('Probe temperature from log, acquisition %d, %s'%(t_+1,k_),np.asarray(r[k_])[t_],rows[m_,1:].mean(axis=0),1e-9)
g=loadmat(data/'phaseMAPr_and_magmap.mat',simplify_cells=True)
with h5py.File(data/'T2MAP_and_img4i_native100.mat') as h:
    mip=np.asarray(h['T2MAP']).transpose(3,2,1,0)
    sos=np.asarray(h['T2MAPsos']).transpose(3,2,1,0)
legacy=np.asarray(r['ROI_2DVSMSI_100'],dtype=int)
close('Fixed PFA coordinate transform',s['C2D'],legacy+[0,1],kind='consistency')
close('Full ROI coordinates unchanged',s['C2D_full'],legacy,kind='consistency')
for idx,vol,variant,coords in [(1,sos,0,legacy+[0,1]),(2,sos,2,legacy),(3,mip,0,legacy+[0,1]),(4,mip,2,legacy)]:
    means=np.column_stack([vol[row-2:row+1,col-2:col+1,:,variant].mean(axis=(0,1)) for row,col in coords])
    close(f'Input ROI extraction {cfg[idx][2]}',cfg[idx][0],means-means.mean(axis=0))
    close(f'Probe-temperature arrays (acquisition-window log means) {cfg[idx][2]}',cfg[idx][1],r['temp1' if variant==0 else 'temp3'],kind='consistency')
phase=np.unwrap(g['phaseMAPr'],axis=2)
dphase=phase-phase.mean(axis=2,keepdims=True)
oil=np.logical_or.reduce([r['oil'+str(i)].astype(bool) for i in range(1,5)])
dphase=dphase-dphase[oil,:].mean(axis=0)[None,None,:]
prfcoords=np.asarray(r['ROI_PRF_128'],dtype=int)
# Effective PRF time = echo-time difference of the dual-echo GRE, read here directly from the DICOM header bytes
# (independent of code/gre_params.m): phase series 0015 is EchoNumber 2, EchoTime 7.46 ms (phase difference);
# magnitude series 0016 is echo 1, EchoTime 5 ms.
def _dicom_te(pattern):
    raw=next(data.glob(pattern)).read_bytes(); i=raw.find(bytes([0x18,0x00,0x81,0x00])+b'DS'); n=int.from_bytes(raw[i+6:i+8],'little')
    return float(raw[i+8:i+8+n].decode().strip())*1e-3
dTE=_dicom_te('GRE_series0015_*.IMA')-_dicom_te('GRE_series0016_*.IMA')
assert abs(dTE-2.46e-3)<1e-9,('GRE dTE from DICOM header',dTE)
# B0 is read independently from the DICOM bytes: ImagingFrequency (0018,0084) / proton gamma-bar 42.577478 MHz/T
# (the nominal 3 T of the FieldStrength tag is NOT used).
def _dicom_ds(pattern,tag):
    raw=next(data.glob(pattern)).read_bytes(); i=raw.find(bytes(tag)+b'DS'); n=int.from_bytes(raw[i+6:i+8],'little')
    return float(raw[i+8:i+8+n].decode().strip())
f0_MHz=_dicom_ds('GRE_series0015_*.IMA',[0x18,0x00,0x84,0x00])
assert abs(f0_MHz-123.251054)<1e-6,('ImagingFrequency from DICOM header',f0_MHz)
GAMMA_BAR=42.577478
B0=f0_MHz/GAMMA_BAR
assert 2.89<B0<2.90,('B0 from the DICOM header should be ~2.8947 T, not the nominal 3 T',B0)
E_B0=(B0,f0_MHz,GAMMA_BAR)
checks.append({'check':'B0 = %.4f T from DICOM ImagingFrequency %.6f MHz / %.6f MHz/T (nominal 3 T not used)'%(B0,f0_MHz,GAMMA_BAR),'passed':True})
prf=-np.column_stack([dphase[row-2:row+1,col-2:col+1,:].mean(axis=(0,1)) for row,col in prfcoords])/(B0*2*np.pi*GAMMA_BAR*1e6*dTE)
close('PRF phase unwrap/oil drift/ROI extraction',cfg[0][0],prf,1e-12)
for i,row in enumerate(cfg):
    y,x,name=row
    x=np.asarray(x);y=np.asarray(y)
    pred=np.empty_like(y)
    for j in range(x.shape[1]):
        for t in range(x.shape[0]):
            train=np.arange(len(x))!=t
            d=np.column_stack((np.ones(train.sum()),x[train,j]))
            intercept,slope=np.linalg.lstsq(d,y[train,j],rcond=None)[0]
            pred[t,j]=(y[t,j]-intercept)/slope
    e=pred-x
    close(f'{name}: all held-out errors',s['loocv_errors'][:,:,i],e,1e-8)
    for metric,expected in [('rmse',np.sqrt(np.mean(e**2,axis=0))),('mae',np.mean(abs(e),axis=0)),('bias',np.mean(e,axis=0)),('max',np.max(abs(e),axis=0))]:
        close(f'{name}: {metric}',s['loocv_'+metric][i,:],expected,1e-8)
    close(f'{name}: aggregation identity',np.sqrt(np.mean(s['loocv_rmse'][i,:]**2)),np.sqrt(np.mean(e**2)))

BI=np.asarray(s['bootstrap_indices'],dtype=int).T-1     # stored resampled time points (20000 x 9)
def boot_slopes(x,y):
    xr=x[BI];yr=y[BI]
    xc=xr-xr.mean(axis=1,keepdims=True);yc=yr-yr.mean(axis=1,keepdims=True)
    valid=np.all(xr.std(axis=1,ddof=1)>=1e-6,axis=1)
    return (xc*yc).sum(axis=(1,2))/(xc*xc).sum(axis=(1,2)),valid
for c,idx in enumerate([1,2,3,4]):
    y,x,name=cfg[idx]; pool=s['POOL'][c]
    beta,rssc=regression(x,y,True); bf,rssf=regression(x,y,False)
    close(f'{name}: common slope',pool['pooled'],beta[-1])
    close(f'{name}: separate slopes',pool['slopes'],bf[-4:])
    assert bool(pool['common_reported'])==(c<2),'common slope must be reported for SOS only'
    if pool['common_reported']:
        bs,valid=boot_slopes(x,y)
        # MATLAB quantile uses sample positions (i-0.5)/n, Hyndman-Fan type 5.
        ci=np.quantile(bs[valid],[.025,.975],method='hazen')
        close(f'{name}: 20000 cluster-bootstrap percentile CI',pool['ci'],ci)
# paired PFA - Full difference: same resampled time points for both variants
yP,xP,_=cfg[1]; yF,xF,_=cfg[2]
bP,vP=boot_slopes(xP,yP); bF,vF=boot_slopes(xF,yF)
dd_=(bP-bF)[vP&vF]
close('PFA-SOS minus Full-SOS common slope',s['PAIR_DIFF'],regression(xP,yP,True)[0][-1]-regression(xF,yF,True)[0][-1])
close('paired cluster-bootstrap CI of the PFA-Full difference',s['PAIR_CI'],np.quantile(dd_,[.025,.975],method='hazen'))
assert int(s['PAIR_NVALID'])==dd_.size
checks.append({'check':'paired bootstrap uses the same stored time-point resamples for PFA and Full','passed':True})

for name,ba,d in [('PFA',s['BA1'],s['dT_pfa']),('Full',s['BA2'],s['dT_full'])]:
    # Summary identity independently checks LOA and RMSE consistency even when
    # registered voxel arrays have been generated by the MATLAB display path.
    sd=(ba['loa_hi']-ba['loa_lo'])/(2*1.96)
    close(name+': LOA midpoint',ba['bias'],(ba['loa_hi']+ba['loa_lo'])/2)
    close(name+': BA RMSE identity',ba['rmse']**2,ba['bias']**2+sd**2*(int(s['n_pairs'])-1)/int(s['n_pairs']))

assert int(s['NB'])==20000 and int(s['SEED'])==20260830
assert int(s['n_vox'])==240 and int(s['n_pairs'])==240*7==1680
assert not any('sd_noise' in p for p in s['POOL'])
checks.append({'check':'seed, replicate count, pair count and noise-only removal','passed':True,'kind':'consistency'})

# ---- reference invariance: dT2/ds referenced to time points 3-9 instead of all nine
for i,row in enumerate(cfg):
    y,x,name=row; y=np.asarray(y); x=np.asarray(x)
    y_alt=y-y[2:9].mean(axis=0)
    ba,_=regression(x,y,False); bb,_=regression(x,y_alt,False)
    close(f'{name}: separate slopes unchanged by the 3-9 reference',ba[-4:],bb[-4:],1e-10)
    if i>0:
        ca,_=regression(x,y,True); cb,_=regression(x,y_alt,True)
        close(f'{name}: common slope unchanged by the 3-9 reference',ca[-1],cb[-1],1e-10)
    pred=np.empty_like(y); pred_alt=np.empty_like(y)
    for j in range(4):
        for t in range(NTP):
            tr=np.arange(NTP)!=t
            d=np.column_stack((np.ones(tr.sum()),x[tr,j]))
            for yy,pp in ((y,pred),(y_alt,pred_alt)):
                a0,b0=np.linalg.lstsq(d,yy[tr,j],rcond=None)[0]; pp[t,j]=(yy[t,j]-a0)/b0
    close(f'{name}: LOOCV predictions unchanged by the 3-9 reference',pred,pred_alt,1e-9)
rc=s['ref_check']
assert max(np.max(rc['max_slope_diff']),np.max(rc['max_loocv_diff']))<1e-7
checks.append({'check':'MATLAB reference-invariance diagnostics (slopes/LOOCV, all-nine vs 3-9 reference) below tolerance','passed':True,'kind':'consistency'})

# ---- every LaTeX macro equals the value recomputed here / stored in the .mat
RN=['I','II','III','IV']; K5=['Gre','PfaSos','FullSos','PfaMip','FullMip']; K4=K5[1:]
E={}
def put(name,fmt,val): E[name]=fmt%val
TT=np.concatenate([np.asarray(r[k_]).ravel() for k_ in ('temprp1','temp1','temp3')])
put('GreWinSec','%.1f',GRE_WINDOW_S); put('GreScanSec','%.0f',27); put('BZero','%.2f',E_B0[0]); put('FZeroMHz','%.2f',E_B0[1]); put('GammaBar','%.3f',E_B0[2])
put('TempMin','%.1f',TT.min()); put('TempMax','%.1f',TT.max()); put('TempSpan','%.1f',TT.max()-TT.min())
put('CoolMin','%.0f',(np.max(r['scan_end_s'])-np.min(r['scan_start_s']))/60)
put('TempStep','%.2f',np.mean(-np.diff(np.hstack([r['temprp1'],r['temp1'],r['temp3']]),axis=0)))
words=['one','two','three','four','five','six','seven','eight','nine','ten']
put('NumTp','%d',NTP); E['NumTpWord']=words[NTP-1]; E['NumTrainWord']=words[NTP-2]
put('NumPredRoi','%d',NTP*4)
put('NumBaTp','%d',7); E['NumBaTpWord']='seven'; put('BaFirstTp','%d',3); put('BaLastTp','%d',9)
put('NumPairs','%d',int(s['n_pairs'])); put('NumBoxVox','%d',int(s['n_vox']))
put('TpCmp','%d',3); E['TpCmpWord']='three'
gap=np.asarray(r['gre_to_pfa_gap_s']); gap2=r['scan_start_s'][:,1]-r['scan_start_s'][:,0]; close('GRE->PFA gap from scan times',gap,gap2,kind='consistency')
put('GapOneSec','%.0f',gap[0]); put('GapOneMin','%.1f',gap[0]/60); put('GapOtherMinSec','%.0f',gap[1:].min()); put('GapOtherMaxSec','%.0f',gap[1:].max())
put('WinRangeMax','%.2f',r['win_range_max']); put('WinDevMax','%.2f',r['win_dev_max'])
for k_ in range(4): put('Dist'+RN[k_],'%.0f',s['dist_mm'][k_])
with h5py.File(data/'T2MAP_and_img4i_native100.mat') as h:
    if 'fit_seconds' in h:
        fs=np.asarray(h['fit_seconds']).transpose(2,1,0).sum(axis=2)   # MATLAB (tp,variant,type) -> per file MIP+SOS
        put('FitMedian','%.1f',np.median(fs)); put('FitMin','%.1f',fs.min()); put('FitMax','%.1f',fs.max())
m2=[s['m2_p_sos'],s['m2_f_sos'],s['m2_p_mip'],s['m2_f_mip']]
IND={}
for m_,row in enumerate(cfg):
    y,x,name=row; y=np.asarray(y); x=np.asarray(x)
    pred=np.empty_like(y)
    for j in range(4):
        for t in range(NTP):
            tr=np.arange(NTP)!=t; d=np.column_stack((np.ones(tr.sum()),x[tr,j]))
            a0,b0=np.linalg.lstsq(d,y[tr,j],rcond=None)[0]; pred[t,j]=(y[t,j]-a0)/b0
    e=pred-x; IND[m_]=e
    for k_ in range(4):
        put('Rmse'+K5[m_]+RN[k_],'%.2f',np.sqrt(np.mean(e[:,k_]**2))); put('Mae'+K5[m_]+RN[k_],'%.2f',np.mean(abs(e[:,k_])))
        put('Bias'+K5[m_]+RN[k_],'%+.2f',np.mean(e[:,k_])); put('Maxe'+K5[m_]+RN[k_],'%.2f',np.max(abs(e[:,k_])))
    put('RmsOverall'+K5[m_],'%.2f',np.sqrt(np.mean(e**2)))
rm=lambda ms:np.array([[np.sqrt(np.mean(IND[m_][:,k_]**2)) for k_ in range(4)] for m_ in ms])
put('RmseMinSos','%.2f',rm([1,2]).min()); put('RmseMaxSos','%.2f',rm([1,2]).max())
put('RmseMinMip','%.2f',rm([3,4]).min()); put('RmseMaxMip','%.2f',rm([3,4]).max())
g_=rm([0])[0]
put('RmseGreRemoteMin','%.2f',g_[[0,3]].min()); put('RmseGreRemoteMax','%.2f',g_[[0,3]].max())
put('RmseGreNearMin','%.2f',g_[[1,2]].min()); put('RmseGreNearMax','%.2f',g_[[1,2]].max())
alpha=np.asarray(s['ALPHA_PAPER'])*1e8
for k_ in range(4):
    put('AlphaGre'+RN[k_],'%.3f',alpha[k_])
    for c_ in range(4): put('Slope'+K4[c_]+RN[k_],'%.3f',m2[c_][k_])
put('AlphaGreRemote','%.3f',alpha[[0,3]].mean()); put('AlphaGreRemoteTwo','%.2f',alpha[[0,3]].mean())
for c_,idx in enumerate([1,2,3,4]):
    y,x,name=cfg[idx]; y=np.asarray(y); x=np.asarray(x); pool=s['POOL'][c_]
    beta,rssc=regression(x,y,True); bf,rssf=regression(x,y,False); sl=bf[-4:]
    put('Cv'+K4[c_],'%.1f',100*np.std(sl,ddof=1)/np.mean(sl))
    if c_<2:
        ci=np.asarray(pool['ci']); rel=(ci[1]-ci[0])/2/beta[-1]
        put('Pool'+K4[c_],'%.4f',beta[-1]); put('PoolTwo'+K4[c_],'%.2f',beta[-1]); put('PoolLo'+K4[c_],'%.4f',ci[0]); put('PoolHi'+K4[c_],'%.4f',ci[1])
        put('RelHw'+K4[c_],'%.2f',100*rel)
        put('UncOne'+K4[c_],'%.2f',rel*1); put('UncEight'+K4[c_],'%.2f',rel*8)
# per-ROI R^2 (slope fit with its own intercept), independent of the MATLAB roifit
for m_,row in enumerate(cfg):
    y,x,name=row; y=np.asarray(y); x=np.asarray(x)
    for k_ in range(4):
        r_=np.corrcoef(x[:,k_],y[:,k_])[0,1]; put('Rsq'+K5[m_]+RN[k_],'%.3f',r_**2)
# between-ROI residual correlations of the SOS common-slope fits
rc_=[]
for idx in (1,2):
    y,x,name=cfg[idx]; y=np.asarray(y); x=np.asarray(x)
    beta,_=regression(x,y,True)
    res=(y-y.mean(axis=0))-beta[-1]*(x-x.mean(axis=0)); cm=np.corrcoef(res.T); rc_+=list(cm[np.triu_indices(4,1)])
put('ResCorrMax','%.2f',max(rc_))
put('PairDiff','%.3f',s['PAIR_DIFF']); put('PairDiffLo','%.3f',s['PAIR_CI'][0]); put('PairDiffHi','%.3f',s['PAIR_CI'][1])
BAt=np.asarray(s['BA_tp']); BS=['Pfa','Full']; BAs=[s['BA1'],s['BA2']]
tw=['One','Two','Three','Four','Five','Six','Seven','Eight','Nine']
for q_ in range(2):
    b=BAs[q_]
    put('BaLoa'+BS[q_],'%.2f',(b['loa_hi']-b['loa_lo'])/2); put('BaBiasPool'+BS[q_],'%.2f',round(float(b['bias']),2)+0.0)
    put('BaLoPool'+BS[q_],'%+.2f',b['loa_lo']); put('BaHiPool'+BS[q_],'%+.2f',b['loa_hi']); put('BaRms'+BS[q_],'%.3f',b['rmse'])
    for it,tp in enumerate(range(3,10)):
        for nm_,col in (('BaBias',0),('BaLo',1),('BaHi',2)): put(nm_+BS[q_]+'Tp'+tw[tp-1],'%+.2f',BAt[it,3*q_+col])
def _deming(xx,pp):
    sxx=np.var(pp);syy=np.var(xx);sxy=np.cov(pp,xx,bias=True)[0,1]
    return (syy-sxx+np.sqrt((syy-sxx)**2+4*sxy**2))/(2*sxy)
aa_,bb_,pp_=[np.asarray(s[k_],dtype=float).ravel() for k_ in ('a_alt','b_alt','p_alt')]
assert list(np.atleast_1d(s['BA_ALT']))==list(range(2,10))
put('BaAltFirstTp','%d',2); put('NumBaAltPairs','%d',aa_.size)
for nm_,xx_,dem_,ci_ in (('Pfa',aa_,s['DEMalt_pfa'],s['DEMalt_pfa_ci']),('Full',bb_,s['DEMalt_full'],s['DEMalt_full_ci'])):
    d_=xx_-pp_; put('BaAltLoa'+nm_,'%.2f',1.96*np.std(d_,ddof=1))
    close('Deming slope with time point 2 included, '+nm_,dem_,_deming(xx_,pp_),1e-9)
    put('BaAltDem'+nm_,'%.2f',dem_); put('BaAltDemLo'+nm_,'%.2f',ci_[0]); put('BaAltDemHi'+nm_,'%.2f',ci_[1])
put('BaBiasMin','%+.2f',BAt[:,[0,3]].min()); put('BaBiasMax','%+.2f',BAt[:,[0,3]].max())
for nm_,v_,ci_ in (('Pfa',s['DEM_pfa'],s['DEM_pfa_ci']),('Full',s['DEM_full'],s['DEM_full_ci'])):
    put('Dem'+nm_,'%.2f',v_); put('DemLo'+nm_,'%.2f',ci_[0]); put('DemHi'+nm_,'%.2f',ci_[1])
TA={'Gre':r['temprp1'],'Pfa':r['temp1'],'Full':r['temp3']}; TR={'Gre':s['Tref_prf'],'Pfa':s['Tref_pfa'],'Full':s['Tref_full']}
for nm_,A in TA.items():
    for k_ in range(4): put('Tcmp'+nm_+RN[k_],'%.1f',A[2,k_])
    put('TcmpMin'+nm_,'%.1f',A[2].min()); put('TcmpMax'+nm_,'%.1f',A[2].max()); put('Tref'+nm_,'%.1f',TR[nm_])
    close('Tref '+nm_+' = mean over time points 3-9 of all ROIs',TR[nm_],A[2:9].mean())
for nm_,v_ in zip(['Pfa','Full','Gre'],np.asarray(s['thr_pct'])): put('ThrPct'+nm_,'%.0f',v_)
# ---- secondary metrics added 2026-11-10: common-coefficient LOOCV, bootstrap intervals, in-sample SD, cooling rate
CIDX=[1,2,3,4]
common_e={}
for q_,idx in enumerate(CIDX):
    y,x,name=cfg[idx]; y=np.asarray(y); x=np.asarray(x); e=np.empty_like(y)
    for t in range(NTP):
        tr=np.arange(NTP)!=t
        beta,_=regression(x[tr],y[tr],True)          # ROI intercepts + one common slope, refitted without time point t
        e[t]=(y[t]-beta[:4])/beta[-1]-x[t]
    common_e[q_]=e
    close('%s: common-coefficient LOOCV errors'%name,np.asarray(s['common_errors'])[:,:,q_],e,1e-8)
    for k_ in range(4): put('RmseCommon'+K4[q_]+RN[k_],'%.2f',np.sqrt(np.mean(e[:,k_]**2)))
    put('RmsOverallCommon'+K4[q_],'%.2f',np.sqrt(np.mean(e**2)))
sets={'GreRemote':IND[0][:,[0,3]],'GreNear':IND[0][:,[1,2]],'PfaSos':IND[1],'FullSos':IND[2],'PfaMip':IND[3],'FullMip':IND[4],
      'CommonPfaSos':common_e[0],'CommonFullSos':common_e[1],'CommonPfaMip':common_e[2],'CommonFullMip':common_e[3]}
assert list(s['boot_names'])==list(sets),'bootstrap set order changed'
rng_=np.random.default_rng(20260830)
for q_,(nm_,E_) in enumerate(sets.items()):
    ix_=rng_.integers(0,NTP,size=(20000,NTP)); bs_=np.sqrt(np.mean(E_[ix_]**2,axis=(1,2)))
    ci_=np.quantile(bs_,[.025,.975],method='hazen')
    # numpy's random stream differs from MATLAB's: agreement is approximate (Monte Carlo error of a 20000-draw percentile)
    close('Overall LOOCV RMSE interval, %s (independent numpy bootstrap, tolerance 0.03 degC)'%nm_,np.asarray(s['boot_ci'])[q_],ci_,0.03)
    close('Overall LOOCV RMSE point value, %s'%nm_,np.asarray(s['boot_rmse'])[q_],np.sqrt(np.mean(E_**2)),1e-8)
    put('RmsBoot'+nm_,'%.2f',np.sqrt(np.mean(E_**2)))
    put('RmsBootLo'+nm_,'%.2f',np.asarray(s['boot_ci'])[q_][0]); put('RmsBootHi'+nm_,'%.2f',np.asarray(s['boot_ci'])[q_][1])
sdfit=np.empty((5,4))
for m_,row in enumerate(cfg):
    y,x,name=row; y=np.asarray(y); x=np.asarray(x)
    for k_ in range(4):
        sl_,ic_=np.polyfit(x[:,k_],y[:,k_],1); sdfit[m_,k_]=np.std((y[:,k_]-ic_)/sl_-x[:,k_],ddof=1)   # in-sample SD, N-1
close('In-sample SD of (fit-derived minus probe temperature), per-ROI fits',np.asarray(s['sd_fit_roi']),sdfit,1e-8)
put('SdFitPfaSosMin','%.2f',sdfit[1].min()); put('SdFitPfaSosMax','%.2f',sdfit[1].max())
put('SdFitFullSosMin','%.2f',sdfit[2].min()); put('SdFitFullSosMax','%.2f',sdfit[2].max())
put('SdFitSosMin','%.2f',sdfit[1:3].min()); put('SdFitSosMax','%.2f',sdfit[1:3].max())
put('SdFitGreRemoteMin','%.2f',sdfit[0,[0,3]].min()); put('SdFitGreRemoteMax','%.2f',sdfit[0,[0,3]].max())
cool_min=(np.max(r['scan_end_s'])-np.min(r['scan_start_s']))/60
put('CoolRate','%.2f',(TT.max()-TT.min())/cool_min); put('CoolFiveMin','%.1f',5*(TT.max()-TT.min())/cool_min)

# satellite files (written by supporting scripts) against their own result files
nd=loadmat(data/'noise_dependence_summary.mat',simplify_cells=True)
if int(nd.get('tp',0))==3:
    for nm_,v_ in zip(K4,np.asarray(nd['snr_replica'])): put('Snr'+nm_,'%.0f',v_)
    put('ReplicasPfa','%d',nd['n_replicas'][0]); put('ReplicasFull','%d',nd['n_replicas'][1]); put('SnrTp','%d',int(nd['tp']))
gr=out/'gre_timepoint_residual.mat'
if gr.exists():
    q=loadmat(gr,simplify_cells=True); sd_=np.asarray(q['resid_sd']); sd2=np.asarray(q['resid_sd_ex2'])
    for t_ in range(9): put('GreResSdTp'+tw[t_],'%.1f',sd_[t_])
    oth=[0,2,3,4,5,6,7,8]
    put('GreResSdOthersMin','%.1f',sd_[oth].min()); put('GreResSdOthersMax','%.1f',sd_[oth].max())
    mg=np.asarray(g['magmap'],dtype=float); w_=10
    edge_=np.concatenate([mg[:w_].ravel(),mg[-w_:].ravel(),mg[:,:w_].ravel(),mg[:,-w_:].ravel()])
    ya_=loadmat(data/'YX.mat',simplify_cells=True)['yangaib'].astype(bool)
    put('GreResNumVox','%d',int(np.sum((mg>=3*edge_.max())&~ya_)))
    put('GreResSdTpOneExTwo','%.1f',sd2[0]); put('GreResSdExTwoMin','%.1f',sd2[oth].min()); put('GreResSdExTwoMax','%.1f',sd2[oth].max())
gs=out/'gre_calibration_sensitivity.mat'
if gs.exists():
    q=loadmat(gs,simplify_cells=True); al=np.asarray(q['alpha'])[1]*1e8; rs=np.asarray(q['rmse'])[1]
    for k_ in range(4): put('GreSensAlpha'+RN[k_],'%.3f',al[k_]); put('GreSensRmse'+RN[k_],'%.2f',rs[k_])
    put('GreSensRmseNearMin','%.2f',rs[[1,2]].min()); put('GreSensRmseNearMax','%.2f',rs[[1,2]].max())
    put('GreSensAlphaRemote','%.3f',al[[0,3]].mean()); put('GreSensAlphaRemoteTwo','%.2f',al[[0,3]].mean())
    close('GRE sensitivity primary set equals the main analysis (all nine time points)',np.asarray(q['alpha'])[0]*1e8,alpha,1e-9,kind='consistency')
    close('GRE sensitivity primary LOOCV RMSE equals the main analysis',np.asarray(q['rmse'])[0],s['loocv_rmse'][0],1e-9,kind='consistency')

def read_macros(path,seen=None):
    seen=seen if seen is not None else set(); path=Path(path)
    if not path.exists() and (base/'results'/path.name).exists():
        # numbers_snr.tex is written only by the raw-data pseudo-replica run; fall back to the reference copy
        print('note: %s not regenerated (needs raw data); using reference copy results/%s'%(path.name,path.name))
        path=base/'results'/path.name
    if path in seen or not path.exists(): return {}
    seen.add(path); d={}
    for ln in path.read_text(encoding='utf-8').splitlines():
        if ln.lstrip().startswith('%'): continue
        m=re.match(r'\\newcommand\{\\(n[A-Za-z]+)\}\{([^{}]*)\}',ln)
        if m:
            assert m.group(1) not in d,('duplicate macro',m.group(1)); d[m.group(1)]=m.group(2)
        m=re.match(r'\\InputIfFileExists\{([^{}]+)\}',ln)
        if m:
            for k_,v_ in read_macros(out.parent/m.group(1),seen).items():
                assert k_ not in d,('duplicate macro across files',k_); d[k_]=v_
    return d
mac=read_macros(out/'numbers.tex')
assert len(mac)>200,('numbers.tex looks incomplete',len(mac))
assert all(re.fullmatch(r'n[A-Za-z]+',k_) for k_ in mac), 'macro names must be letters only'
wrong=[(k_,mac[k_],E['%s'%k_[1:]]) for k_ in mac if k_[1:] in E and mac[k_]!=E[k_[1:]]]
assert not wrong,('macro value differs from recomputed value',wrong[:10])
unchecked=[k_ for k_ in mac if k_[1:] not in E]
assert not unchecked,('macros without a recomputation in verify_statistics.py. Fix: add them to the macro block',unchecked[:20])
missing=[k_ for k_ in E if 'n'+k_ not in mac]
assert not missing,('expected macros missing from numbers*.tex. Fix: rerun statistics_analysis_final.m and the supporting scripts',missing[:20])
checks.append({'check':'%d LaTeX macros equal recomputed / stored values'%len(mac),'passed':True,'kind':'consistency'})
report={'passed':len(checks),'checks':checks,'scope':'Independent least-squares LOOCV, per-ROI slopes and R2, stored cluster resamples (common-slope CI and paired PFA-Full difference), B0 from the DICOM header and BA identities.'}
(out/'statistics_independent_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
n_ind=sum(1 for c_ in checks if c_.get('kind','independent')=='independent'); n_con=len(checks)-n_ind
report['independently_recomputed']=n_ind; report['consistency']=n_con
(out/'statistics_independent_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(f'PASS: {n_ind} independently recomputed checks + {n_con} consistency checks (value copied between files)')
