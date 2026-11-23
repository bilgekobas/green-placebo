"""Rebuild EDA and PPG-derived HRV from original Empatica E4 exports.

Requires: python 3.11+, numpy, pandas, scipy, neurokit2, heartpy, openpyxl.
Run: python 04_reanalysis/03_reprocess_eda_bvp.py
All outputs go to 04_reanalysis/outputs and never overwrite original inputs.

EDA is NeuroKit2-derived at the correct 4-Hz sample rate. BVP-based HRV is
a research-grade candidate requiring visual beat/signal QC against raw BVP;
it is NOT automatically validated clinical ECG HRV.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import neurokit2 as nk
import heartpy as hp

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'00_data/raw/EMPATICA'
OUT=ROOT/'04_reanalysis/outputs'
OUT.mkdir(parents=True,exist_ok=True)
sessions=pd.read_excel(ROOT/'00_data/00_Sessions.xlsx')
required={'Participant','SessionID','Test','Session','Start','End'}
if not required.issubset(sessions.columns):
    raise ValueError(f'Missing session columns: {required-set(sessions.columns)}')
# Historical recording dates were stored as local Munich (Europe/Berlin) times.
for col in ('Start','End'):
    v=pd.to_datetime(sessions[col],errors='coerce')
    if v.isna().any():raise ValueError(f'Invalid session {col}')
    sessions[col+'_unix']=v.dt.tz_localize('Europe/Berlin',ambiguous='raise',nonexistent='raise').astype('int64')/1e9
sessions['local_date']=sessions['Start'].dt.date
records=[]; windows=[]; inventory=[]
for folder in sorted(RAW.iterdir()):
    if not folder.is_dir():continue
    try:
        recording_epoch=float(folder.name.split('_')[0])
    except ValueError:continue
    date=pd.Timestamp(recording_epoch,unit='s',tz='UTC').tz_convert('Europe/Berlin').date()
    current=sessions[sessions.local_date==date]
    if current.Participant.nunique()!=1:
        inventory.append(dict(folder=folder.name,reason='unmatched_or_ambiguous_date'))
        continue
    part=current.Participant.iloc[0]
    item=dict(folder=folder.name,Participant=part,eda_available=(folder/'EDA.csv').exists(),bvp_available=(folder/'BVP.csv').exists())
    inventory.append(item)
    eda=folder/'EDA.csv'
    if eda.exists():
        try:
            arr=np.loadtxt(eda,delimiter=',')
            t0,fs=float(arr[0]),float(arr[1]); y=arr[2:]
            if fs!=4: raise ValueError('EDA sample rate not 4Hz')
            if len(y)>5 and np.isfinite(y).all():
                signals,_=nk.eda_process(y,sampling_rate=4,method='neurokit')
                t=t0+np.arange(len(y))/4
                s=pd.DataFrame({'unix':t,'EDA_Raw':y,'EDA_Clean':signals.EDA_Clean.to_numpy(),
                                'EDA_Tonic':signals.EDA_Tonic.to_numpy(),'EDA_Phasic':signals.EDA_Phasic.to_numpy()})
                s['minute_unix']=np.floor(s.unix/60)*60
                m=s.groupby('minute_unix',as_index=False).agg(
                    EDA_Raw=('EDA_Raw','mean'),EDA_Clean=('EDA_Clean','mean'),
                    EDA_Tonic=('EDA_Tonic','mean'),EDA_Phasic=('EDA_Phasic','mean'),
                    n_eda_samples=('unix','size'))
                for _,ss in current.iterrows():
                    clip=m[((m.minute_unix+30)>=ss.Start_unix)&((m.minute_unix+30)<ss.End_unix)].copy()
                    clip['Participant']=part;clip['SessionID']=ss.SessionID
                    clip['Condition']=int(ss.Test);clip['SessionPosition']=int(ss.Session)
                    clip['EDA_quality']='coverage_only_not_visual_qc'
                    records.append(clip)
        except Exception as exc:
            inventory[-1]['eda_processing_error']=str(exc)[:150]
    bvp=folder/'BVP.csv'
    if bvp.exists():
        try:
            arr=np.loadtxt(bvp,delimiter=',')
            t0,fs=float(arr[0]),float(arr[1]); y=arr[2:]
            if fs!=64:raise ValueError('BVP sample rate not 64 Hz')
            # Original published filter specification: order 3, 0.75-3.5 Hz
            filt=sosfiltfilt(butter(3,[.75,3.5],btype='bandpass',fs=64,output='sos'),y)
            # Window separately by session to prevent data from different conditions
            # being included in the same HRV calculation.
            for _,ss in current.iterrows():
                a=max(0,int(np.ceil((ss.Start_unix-t0)*64)))
                b=min(len(filt),int(np.floor((ss.End_unix-t0)*64)))
                if b-a<240*64: continue
                # Avoid first 15 minutes of session (stabilization).
                a+=15*60*64
                if b-a<240*64:continue
                # 240 s windows / 30 s steps as in the historical script
                for start in range(a,b-240*64+1,30*64):
                    seg=filt[start:start+240*64]
                    result=dict(Participant=part,SessionID=ss.SessionID,Condition=int(ss.Test),
                                SessionPosition=int(ss.Session),window_start_unix=t0+start/64,
                                quality='unreviewed',reason='',bpm=np.nan,rmssd_ms=np.nan,
                                hf_ms2=np.nan,lf_ms2=np.nan,lf_over_lf_plus_hf=np.nan,
                                rejected_peak_fraction=np.nan)
                    try:
                        work,meas=hp.process(seg,sample_rate=64,calc_freq=True)
                        bpm=float(meas.get('bpm',np.nan))
                        hf=float(meas.get('hf',np.nan));lf=float(meas.get('lf',np.nan))
                        rmssd=float(meas.get('rmssd',np.nan))
                        # HeartPy raw peak rejection does not account for all motion artifacts.
                        nall=len(work.get('peaklist',[]))
                        bad=len(work.get('removed_beats',[]))
                        frac=bad/nall if nall else np.nan
                        result.update(bpm=bpm,rmssd_ms=rmssd,hf_ms2=hf,lf_ms2=lf,
                                      rejected_peak_fraction=frac)
                        if 40<=bpm<=180 and 0<=frac<=0.05 and np.isfinite(rmssd) and rmssd>0 and np.isfinite(lf) and lf>=0 and np.isfinite(hf) and hf>0:
                            result.update(quality='candidate_pass_requires_signal_review',
                                          lf_over_lf_plus_hf=lf/(lf+hf))
                        else:result.update(quality='reject',reason='physiological_or_peak_qc')
                    except Exception as exc:
                        result.update(quality='reject',reason=('heartpy_error:'+type(exc).__name__))
                    windows.append(result)
        except Exception as exc:
            inventory[-1]['bvp_processing_error']=str(exc)[:150]

E=pd.concat(records,ignore_index=True) if records else pd.DataFrame()
H=pd.DataFrame(windows);I=pd.DataFrame(inventory)
E.to_csv(OUT/'eda_neurokit2_4hz_1min.csv',index=False)
H.to_csv(OUT/'bvp_heartpy_hrv_window_qc.csv',index=False)
I.to_csv(OUT/'raw_reprocessing_inventory.csv',index=False)
summary={'neurokit2_version':getattr(nk,'__version__','unknown'),
         'heartpy_version':getattr(hp,'__version__','unknown'),
         'eda_minutes':len(E),'eda_sessions':E.SessionID.nunique() if len(E) else 0,
         'hrv_windows':len(H),'hrv_candidate_windows':int(H.quality.eq('candidate_pass_requires_signal_review').sum()) if len(H) else 0,
         'hrv_candidate_sessions':H.loc[H.quality.eq('candidate_pass_requires_signal_review'),'SessionID'].nunique() if len(H) else 0}
(OUT/'processing_summary.json').write_text(json.dumps(summary,indent=2))
print(summary)
