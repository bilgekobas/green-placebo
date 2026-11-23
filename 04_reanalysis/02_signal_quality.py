"""Quality audit of reconstructed raw Empatica recordings.

Run AFTER rebuilding the raw data with the scripts supplied in the
raw-rebuild archive:
    python 04_reanalysis/02_signal_quality.py

Inputs: 04_reanalysis/outputs/minute_raw_rebuilt.csv or /mnt/data/rebuilt_raw/
and repository 00_data/raw/EMPATICA.
No replacement of original processed data.
Pulse IBI results are not ECG-validated HRV and must pass adequate coverage.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import welch

ROOT=Path(__file__).resolve().parents[1]
DEST=Path(__file__).resolve().parent/'outputs'
DEST.mkdir(exist_ok=True)
INPUT=DEST/'minute_raw_rebuilt.csv'
if not INPUT.exists():
    INPUT=Path('/mnt/data/rebuilt_raw/minute_raw_rebuilt.csv')
M=pd.read_csv(INPUT)
M['eda_usable']=(M.EDA_n_samples>=228)&M.EDA_raw_mean.ge(0)&M.EDA_raw_mean.notna()
M['temp_usable']=(M.TEMP_n_samples>=228)&M.TEMP_raw_mean.between(15,45)
M['hr_usable']=(M.HR_n_samples>=57)&M.HR_raw_mean.between(40,180)
M.to_csv(DEST/'eda_temp_hr_qc_minute.csv',index=False)
P=M.groupby('SessionID',as_index=False).agg(
    Participant=('Participant','first'),
    Condition=('Condition','first'),
    SessionPosition=('SessionPosition','first'),
    eda_valid_minutes=('eda_usable','sum'),
    temperature_valid_minutes=('temp_usable','sum'),
    hr_valid_minutes=('hr_usable','sum'))
folders={}
for f in (ROOT/'00_data/raw/EMPATICA').iterdir():
    if not f.is_dir(): continue
    name=M.loc[M.source_folder==f.name,'Participant'].dropna().unique()
    if len(name)==1:folders[name[0]]=f
out=[]
for participant,folder in folders.items():
    if not (folder/'IBI.csv').exists():continue
    with (folder/'IBI.csv').open() as fp:
        epoch=float(fp.readline().split(',')[0])
    raw=np.genfromtxt(folder/'IBI.csv',delimiter=',',skip_header=1)
    if raw.ndim!=2 or raw.shape[1]!=2:continue
    times=epoch+raw[:,0]
    ibi=raw[:,1]
    for sid,g in M[M.Participant==participant].groupby('SessionID'):
        start=g.minute_unix.min()+900
        end=g.minute_unix.max()+60
        for ws in np.arange(start,end-239,240):
            mask=(times>=ws)&(times<ws+240)
            t,r=times[mask],ibi[mask]
            item=dict(Participant=participant,SessionID=sid,Condition=int(g.Condition.iloc[0]),
                      SessionPosition=int(g.SessionPosition.iloc[0]),window_start_unix=float(ws),
                      n_intervals=int(len(r)),quality='reject',reason='',HF_ms2=np.nan,
                      LF_ms2=np.nan,RMSSD_ms=np.nan,LF_over_LFplusHF=np.nan)
            if len(r)<160:
                item['reason']='insufficient_intervals'
            else:
                ok=np.isfinite(r)&(r>=.33)&(r<=1.5)
                if ok.mean()<.95:
                    item['reason']='invalid_ibi_fraction'
                else:
                    t,r=t[ok],r[ok]
                    if len(r)<160 or t[-1]-t[0]<225 or (np.diff(t)>3).any():
                        item['reason']='large_gaps_or_short_coverage'
                    else:
                        bad=np.r_[False,np.abs(np.diff(r))>np.maximum(.2,.2*r[:-1])]
                        t,r=t[~bad],r[~bad]
                        if bad.mean()>.05 or len(r)<160 or (np.diff(t)>3).any():
                            item['reason']='artifact_or_post_artifact_gap'
                        else:
                            it=np.arange(t[0],t[-1],.25)
                            if len(it)<800:
                                item['reason']='short_interpolated_record'
                            else:
                                y=np.interp(it,t,r)*1000
                                fq,px=welch(y,fs=4,nperseg=min(1024,len(y)),detrend='linear')
                                lm=(fq>=.04)&(fq<.15);hm=(fq>=.15)&(fq<.40)
                                lf=float(np.trapezoid(px[lm],fq[lm]))
                                hf=float(np.trapezoid(px[hm],fq[hm]))
                                if not (np.isfinite(lf+hf) and lf>=0 and hf>0):
                                    item['reason']='invalid_spectrum'
                                else:
                                    item.update(quality='pass_preliminary_ibi_qc',
                                      HF_ms2=hf,LF_ms2=lf,
                                      LF_over_LFplusHF=lf/(lf+hf),
                                      RMSSD_ms=float(np.sqrt(np.mean(np.diff(r)**2))*1000))
            out.append(item)
W=pd.DataFrame(out)
W.to_csv(DEST/'ibi_hrv_strict_qc_windows.csv',index=False)
accepted=W[W.quality=='pass_preliminary_ibi_qc']
if len(accepted):
    H=accepted.groupby('SessionID',as_index=False).agg(
        hrv_valid_windows=('HF_ms2','size'),HF_median_ms2=('HF_ms2','median'),
        LF_median_ms2=('LF_ms2','median'),RMSSD_median_ms=('RMSSD_ms','median'))
    P=P.merge(H,on='SessionID',how='left')
P.to_csv(DEST/'signal_quality_session_summary.csv',index=False)
print('Raw EDA usable minutes:',int(M.eda_usable.sum()))
print('Raw EDA sessions:',M.loc[M.eda_usable,'SessionID'].nunique())
print('HRV windows passing initial QC:',len(accepted),'of',len(W))
print('HRV sessions passing initial QC:',accepted.SessionID.nunique())
