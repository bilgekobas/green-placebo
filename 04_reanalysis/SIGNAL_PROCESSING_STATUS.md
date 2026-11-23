# Corrected signal processing: execution and validation status

**As of 2026-11-23, no complete clinically validated EDA/HRV dataset has been produced.** Original data remain untouched.

## Executed locally

- Earlier independent reconstruction: 10,680 raw-derived minute records from 120/132 sessions; 116 sessions had EDA.
- New raw EDA coverage QC: 10,301 minute records passed >=95% raw sample-count and nonnegative-EDA checks, across 116 sessions. These are *signal availability flags*, not a manual sensor-quality assessment.
- Raw HR quality screen: 120 sessions had coverage for further examination.
- Initial conservative Empatica IBI-based HRV screen: **7 / 2,143** candidate 240-second windows passed the completeness/gap checks, across only **3** sessions. **Do not use this HRV subset for condition comparison.** Most rejected windows had insufficient beat intervals, not necessarily bad underlying BVP.

## Original processing problems to resolve

1. The historical EDA script constructed timestamps at 4 Hz but passed `sampling_rate=8` to `nk.eda_process`. For a corrected decomposition use the 4-Hz raw EDA (`03_reprocess_eda_bvp.py`).
2. The original HRV procedure used 64-Hz BVP; Butterworth bandpass 0.75–3.5 Hz (order 3); HeartPy 240-second segments stepped every 30 seconds. Corrected script reproduces these basic settings and records candidate peak-quality flags; **candidate** estimates still require inspection of BVP and beat detections and cannot be considered ECG-validated.
3. Exact analysis duration remains unresolved: 90 minutes minus first 15 and final 30 gives 45 minutes, not 60.
4. Raw EDA files are not available for all study sessions; final exclusions must be disclosed.

## GitHub files

- `01_eda_session_model.py`: previously committed first-pass statistical analysis
- `02_signal_quality.py`: raw minute and IBI quality-control screening
- `03_reprocess_eda_bvp.py`: corrected NeuroKit2 and HeartPy *candidate* processing
- `.github/workflows/reprocess-empatica.yml`: **manual** GitHub Actions workflow, designed to execute in an environment where dependencies can be installed; generates versioned output files in `04_reanalysis/outputs` and pushes them to `reviewer-reanalysis`.

## Why the GitHub Action is not yet a validated analysis

The author must run the workflow from the GitHub Actions tab. This environment could not install NeuroKit2 and HeartPy (package index connectivity unavailable), so the corrected NeuroKit2/HeartPy scripts have **not yet been executed or checked against their output files**. Do not describe those outputs as available until GitHub Actions successfully runs. Manual signal-quality review and completeness audit are required before inferential use.

The independent local reconstruction archive contains existing raw-derived CSVs, scripts, and logs; larger local outputs are **not yet pushed to GitHub**. The workflow will push its *new* outputs after a successful run.
