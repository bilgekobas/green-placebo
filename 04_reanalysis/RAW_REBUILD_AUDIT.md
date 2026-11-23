# Raw-data reconstruction audit — 23 November 2026

Status: **partial independent reconstruction; not yet a validated replacement for the publication dataset**.

The original Empatica archive and supplied `00_Sessions.xlsx` were used to rebuild 1-minute data from raw EDA (4 Hz), wrist temperature (4 Hz), and E4 heart rate (1 Hz). An exploratory 4-Hz EDA filtering + 0.05-Hz tonic estimate was used because the originally specified NeuroKit2 and HeartPy packages were unavailable in the analysis environment. This is **not** a reproduction of the NeuroKit2 tonic decomposition or HeartPy pulse HRV algorithm.

**Coverage:** 33 Empatica folders, 29 with raw EDA, 29 with raw TEMP, 30 with raw HR, 30 with raw IBI, 27 with raw BVP. Rebuilt E4 minute data cover **120/132 sessions, 10,680 minute rows**, with 10,323 raw EDA and wrist-temperature minutes. Missing raw E4 for **P25, P30 and P33**; P28 has HR but lacks raw EDA and TEMP. Independently recovered raw IAQ exports cover **30 Nov 2022–5 Jan 2023 only**, yielding 3,920 minute rows with IAQ readings; other IEQ data may be present in processed datasets but have **not** been recovered from raw files.

**Validation:** 10,656 participant-session minute entries can be aligned to the original clean dataset after accounting for the original midpoint convention (Unix minute start +30 seconds). Raw EDA matches original EDA_Raw exceptionally closely; this agreement alone does not validate tonic decomposition.

**Exploratory pulse HRV:** IBI-based 240-second nonoverlapping windows using interpolation and Welch spectra were produced, but only 580/2,635 windows passed basic preliminary completeness/range checks (84 sessions); these should **not** be substituted for validated HeartPy BVP-derived HF/LF metrics.

**Questionnaires:** The original PANAS, comfort, cognitive spreadsheets retain 132 session-level records; they were not independently re-scored from source forms.

**Next:** recover missing raw E4 and later raw IAQ exports; install and pin NeuroKit2 and HeartPy; recompute EDA at true 4 Hz, audit BVP beat quality, confirm exact session analysis window, then undertake the prespecified mixed-effects models.

A reproducible ZIP of scripts, data extracts and a full QA README was produced in this analysis session; the original data and repository main branch were not changed. To avoid accidentally publishing incomplete analytical substitutes, the derived CSVs have **not** been committed to GitHub.
