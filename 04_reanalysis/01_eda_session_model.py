"""Exploratory first-pass tonic EDA reanalysis.

Run from repository root:
    python 04_reanalysis/01_eda_session_model.py

Requires pandas, numpy, scipy, statsmodels, openpyxl. Outputs a CSV of
condition contrasts by analysis window and sample. This script deliberately
does not alter original data or figures.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.formula.api import mixedlm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
data = pd.read_csv(ROOT / "00_data" / "10_MergedAllClean.csv")
meta = pd.read_excel(ROOT / "00_data" / "00_Meta.xlsx")
excluded = set(meta.loc[meta["Noticed artificial plants"] == 1, "Participant"])
contrasts = [(2, 1), (2, 3), (3, 1), (4, 1)]
rows = []

for end_minute in (74, 59):
    for exclude_unblinded in (False, True):
        subset = data.loc[data["ElapsedMinutes"].between(15, end_minute)].copy()
        if exclude_unblinded:
            subset = subset.loc[~subset["Participant"].isin(excluded)]
        sessions = (
            subset.groupby(["SessionID", "Participant", "Test", "Session"], as_index=False)
            .agg(eda=("EDA_Tonic", "mean"), n_minutes=("EDA_Tonic", "size"))
        )
        if (sessions["eda"] <= 0).any():
            raise ValueError("Nonpositive EDA session mean; log transform not valid")
        sessions["logeda"] = np.log(sessions["eda"])
        formula = "logeda ~ C(Test, Treatment(reference=1)) + C(Session)"
        model = mixedlm(formula, sessions, groups=sessions["Participant"])
        fit = model.fit(reml=False, method="lbfgs")
        if not fit.converged:
            raise RuntimeError("Mixed model failed to converge")
        coef = fit.params
        covariance = fit.cov_params()
        pvalues = []
        current = []
        for condition_a, condition_b in contrasts:
            v = pd.Series(0.0, index=coef.index)
            for condition, sign in ((condition_a, 1), (condition_b, -1)):
                if condition != 1:
                    key = f"C(Test, Treatment(reference=1))[T.{condition}]"
                    v[key] += sign
            log_ratio = float(v @ coef)
            se = float(np.sqrt(v @ covariance @ v))
            p = float(2 * norm.sf(abs(log_ratio / se)))
            current.append(dict(
                last_minute=end_minute,
                exclude_unblinded=exclude_unblinded,
                n_participants=sessions["Participant"].nunique(),
                n_sessions=len(sessions),
                condition_a=condition_a,
                condition_b=condition_b,
                log_ratio=log_ratio,
                ratio=np.exp(log_ratio),
                ci_low=np.exp(log_ratio - 1.96 * se),
                ci_high=np.exp(log_ratio + 1.96 * se),
                p_unadjusted=p,
            ))
            pvalues.append(p)
        p_holm = multipletests(pvalues, method="holm")[1]
        for item, adjusted_p in zip(current, p_holm):
            item["p_holm_four_contrasts"] = float(adjusted_p)
            rows.append(item)

pd.DataFrame(rows).to_csv(OUT / "eda_session_model_contrasts.csv", index=False)
print(pd.DataFrame(rows).round(4).to_string(index=False))
