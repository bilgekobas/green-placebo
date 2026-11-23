# Reviewer reanalysis — first-pass tonic EDA audit (23 November 2026)

## Source and status

Exploratory, **not submission-ready**. Uses `00_data/10_MergedAllClean.csv`, `00_data/00_Meta.xlsx`, and the existing `01_codes/Vis.py`. Reproducible script: [01_eda_session_model.py](01_eda_session_model.py). No original files or figures were changed.

## Reproduction concerns

1. `Vis.py` computes `EDA_Tonic_Z` with `groupby('SessionID')`. This is **within participant-session** standardization, not across all four sessions per participant as described in the manuscript. The mean of the resulting standardized values within each session is mathematically forced toward zero; hence differences of session means cannot be assessed using this variable.
2. The manuscript's window specification is ambiguous: first 15 minutes plus final 30 minutes excluded from each 90-minute session gives **45**, not 60, analyzed minutes. We therefore checked minutes **15–74** (60 minutes) and **15–59** (45 minutes); the exact protocol window needs confirmation.
3. The data include all **33 participants × 4 sessions**, and five were marked as noticing artificial plants in `00_Meta.xlsx`: P04, P10, P24, P28, P31. The primary all-participant analysis and the exclusion sensitivity are both reported.

## Model

For each participant-session, average positive unnormalized tonic EDA in microSiemens. Model `log(mean_EDA) ~ condition + session_position + (1 | participant)` by Gaussian mixed effects (maximum likelihood). This models the **session average** without treating minutes as independent replicates. Coefficients are exponentiated as ratios of geometric means. Condition contrasts are S2 vs S1, S2 vs S3, S3 vs S1, and S4 vs S1. Holm adjustment across these four *specified exploratory* contrasts, within each sample/window. The family of confirmatory comparisons must be specified before final inference.

## Results: 60-minute window, all participants (N=33)

| Contrast | Geometric-mean ratio | Approx. 95% CI | p (nominal) | p (Holm, four contrasts) |
|---|---:|---:|---:|---:|
| S2 real / S1 no plants | 0.754 | 0.604–0.942 | 0.0128 | 0.0512 |
| S2 real / S3 artificial | 0.758 | 0.605–0.950 | 0.0162 | 0.0512 |
| S3 artificial / S1 no plants | 0.995 | 0.800–1.238 | 0.966 | 0.966 |
| S4 window / S1 no plants | 0.906 | 0.703–1.167 | 0.4442 | 0.8884 |

## Sensitivity: 60-minute window, excluding five participants (N=28)

| Contrast | Ratio | p (nominal) | p (Holm) |
|---|---:|---:|---:|
| S2 / S1 | 0.776 | 0.0285 | 0.0928 |
| S2 / S3 | 0.768 | 0.0232 | 0.0928 |

For the **45-minute window**, 33-participant estimates were S2/S1 = 0.764 (p=0.020; Holm=0.060) and S2/S3 = 0.742 (p=0.0111; Holm=0.0444). In the 28-participant 45-minute sensitivity, S2/S1 = 0.790 (p=0.0496; Holm=0.1488) and S2/S3 = 0.760 (p=0.0224; Holm=0.0896).

## Interpretation and remaining checks

Lower EDA associated with real plants **persists directionally** after adjustment for session position; multiplicity-adjusted evidence is borderline or non-significant across tested specifications. This does **not** prove causal biological plant effects, stress reduction, or a successful placebo: visual morphology and tasks differ, session allocation is not fully balanced, and sample exclusions were based on a post-exposure observation.

Before manuscript revision: confirm actual analyzable window and origin of `StartEnd`, check signal quality and missingness, compare against original figure-generation steps, evaluate residuals and random-effect fit, prespecify a primary contrast and multiplicity family, examine possible time-of-day/sequence confounding, and perform participant-level diagnostics. Models here are preliminary Gaussian random-intercept models; their Wald CIs/p-values should not yet be treated as final.
