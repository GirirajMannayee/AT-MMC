"""Descriptive / mixture-level statistics (Section 4.8, 4.13.4, 4.14).

IMPORTANT: all calculations use ONE value per mixture (or per mixture x exposure
cell).  They are descriptive arithmetic.  Their F- and p-values are NOMINAL and
are not evidence of experimental significance (see manuscript Sections 3.9, 4.14).
Replicate-based inference lives in at_mmc_anova.py and needs specimen-level data.
"""
import numpy as np
import pandas as pd
from scipy import stats as st
from . import data as D
from .config import EXPOSURES, MIXES, PROPS

FIBRE_PROPS = ["density_kg_m3", "water_absorption_pct", "apparent_porosity_pct", "sorptivity_mm_min05",
               "compressive_MPa", "flexural_MPa", "splitting_tensile_MPa", "fracture_toughness_MPa_m05",
               "acid_5pct_28d"]


# ---------------------------------------------------------------- Table 18
def fibre_dose_fits() -> pd.DataFrame:
    """OLS of AT1-AT4 means on basalt-fibre dose (0, 0.5, 1.0, 1.5 mass %). F(1,2), nominal p."""
    df = D.load_all().set_index("mix_id")
    x = df.loc[["AT1", "AT2", "AT3", "AT4"], "basalt_fibre_pct"].astype(float).values
    rows = []
    for key in FIBRE_PROPS:
        y = df.loc[["AT1", "AT2", "AT3", "AT4"], key].astype(float).values
        res = st.linregress(x, y)
        r2 = res.rvalue ** 2
        F = r2 / (1 - r2) * 2
        rows.append({"property": key, "slope_per_1pct_BF": res.slope, "intercept_0pct_BF": res.intercept,
                     "R2_descriptive": r2, "F_1_2_nominal": F, "p_nominal": st.f.sf(F, 1, 2)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Table 19 / Kendall W
def two_way_unreplicated() -> pd.DataFrame:
    """Mixture x exposure ANOVA on Table 5 (one value per cell; residual = interaction)."""
    d = D.load_durability().set_index("mix_id")[list(EXPOSURES)]
    Y = d.values.astype(float)
    a, b = Y.shape
    g = Y.mean()
    ss_mix = b * ((Y.mean(axis=1) - g) ** 2).sum()
    ss_exp = a * ((Y.mean(axis=0) - g) ** 2).sum()
    ss_tot = ((Y - g) ** 2).sum()
    ss_res = ss_tot - ss_mix - ss_exp
    df_mix, df_exp, df_res = a - 1, b - 1, (a - 1) * (b - 1)
    ms_mix, ms_exp, ms_res = ss_mix / df_mix, ss_exp / df_exp, ss_res / df_res
    rows = [
        {"source": "Mixture", "SS": ss_mix, "df": df_mix, "MS": ms_mix, "F_nominal": ms_mix / ms_res,
         "p_nominal": st.f.sf(ms_mix / ms_res, df_mix, df_res), "eta2": ss_mix / ss_tot},
        {"source": "Exposure", "SS": ss_exp, "df": df_exp, "MS": ms_exp, "F_nominal": ms_exp / ms_res,
         "p_nominal": st.f.sf(ms_exp / ms_res, df_exp, df_res), "eta2": ss_exp / ss_tot},
        {"source": "Mixture x exposure (residual)", "SS": ss_res, "df": df_res, "MS": ms_res,
         "F_nominal": np.nan, "p_nominal": np.nan, "eta2": ss_res / ss_tot},
        {"source": "Total", "SS": ss_tot, "df": a * b - 1, "MS": ss_tot / (a * b - 1),
         "F_nominal": np.nan, "p_nominal": np.nan, "eta2": 1.0},
    ]
    return pd.DataFrame(rows)


def kendall_w_and_friedman() -> dict:
    """Concordance of the mixture ranking across the five exposure endpoints."""
    d = D.load_durability().set_index("mix_id")[list(EXPOSURES)]
    ranks = d.rank(axis=0)  # ranks of the 7 mixtures within each exposure (ties averaged)
    m, n = ranks.shape[1], ranks.shape[0]
    Rj = ranks.sum(axis=1)
    S = ((Rj - Rj.mean()) ** 2).sum()
    # tie correction
    T = 0.0
    for c in ranks.columns:
        _, counts = np.unique(ranks[c], return_counts=True)
        T += (counts ** 3 - counts).sum()
    W = 12 * S / (m ** 2 * (n ** 3 - n) - m * T)
    chi2 = m * (n - 1) * W
    fr = st.friedmanchisquare(*[d[c].values for c in d.columns])
    return {"kendall_W": W, "chi2": chi2, "df": n - 1, "p_nominal": st.chi2.sf(chi2, n - 1),
            "friedman_stat_scipy": fr.statistic, "friedman_p_scipy": fr.pvalue,
            "n_mixtures": n, "n_exposures": m}


# ---------------------------------------------------------------- Tables 15 / 21
def ryshkewitch_fits() -> pd.DataFrame:
    """ln(y) = a + b*P (P in %).  y0 = exp(a); b_fraction = -100*b  (Eq. 5, P as fraction)."""
    df = D.load_all().set_index("mix_id")
    rows = []
    for key in ["compressive_MPa", "splitting_tensile_MPa", "fracture_toughness_MPa_m05", "flexural_MPa"]:
        for label, mixes in [("C0 + AT1-AT6", MIXES), ("AT1-AT6", MIXES[1:])]:
            P = df.loc[mixes, "apparent_porosity_pct"].astype(float).values
            y = np.log(df.loc[mixes, key].astype(float).values)
            r = st.linregress(P, y)
            rows.append({"response": key, "data_set": label, "n": len(mixes), "intercept_a": r.intercept,
                         "slope_b_per_pct": r.slope, "R2_of_ln_y": r.rvalue ** 2,
                         "y0_fitting_constant": float(np.exp(r.intercept)), "b_Eq5_fraction": -100 * r.slope})
    return pd.DataFrame(rows)


def porosity_acid_fit() -> dict:
    df = D.load_all()
    r = st.linregress(df["apparent_porosity_pct"].astype(float), df["acid_5pct_28d"].astype(float))
    return {"slope": r.slope, "intercept": r.intercept, "R2": r.rvalue ** 2, "n": len(df)}


# ---------------------------------------------------------------- Table 22
def contrasts_2x2() -> pd.DataFrame:
    """AT3-AT6 half-difference contrasts. Rows AT3,AT4 (package A) / AT5,AT6 (package B)."""
    df = D.load_all().set_index("mix_id")
    rows = []
    for key in FIBRE_PROPS:
        v = df[key].astype(float)
        a3, a4, a5, a6 = v["AT3"], v["AT4"], v["AT5"], v["AT6"]
        rows.append({"property": key,
                     "binder_package (AT5,AT6 - AT3,AT4)": ((a5 + a6) - (a3 + a4)) / 2,
                     "fibre 1.5 vs 1.0": ((a4 + a6) - (a3 + a5)) / 2,
                     "interaction": ((a3 + a6) - (a4 + a5)) / 2})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Table 23 / Fig. 14d
def replicate_precision(n_list=(3, 5), cv_f=5.0, cv_p=10.0, alpha=0.05) -> pd.DataFrame:
    """Hypothetical planning scenarios.

    Treat the 7 mixture means as the group means of a balanced one-way design with n
    specimens each and a common within-mixture CV, i.e. SD_i = CV * mean_i, so the pooled
    error variance is CV^2 * mean(mean_i^2):
        F = n * var(means, ddof=1) / (CV^2 * mean(means^2)) ,  df = (6, 7(n-1)).
    Returns the largest CV that still gives p < alpha, plus F at CV 5 % and p at CV 10 % (n=3).
    NOT measured uncertainty and NOT evidence of significance.
    """
    df = D.load_all()
    rows = []
    for key in FIBRE_PROPS:
        m = df[key].astype(float).values
        gm, v = np.sqrt(np.mean(m ** 2)), m.var(ddof=1)  # gm = RMS of means
        row = {"property": key}
        for n in n_list:
            df2 = 7 * (n - 1)
            fcrit = st.f.ppf(1 - alpha, 6, df2)
            row[f"max_CV_pct_n{n}"] = 100 * np.sqrt(n * v / fcrit) / gm
        n = n_list[0]
        df2 = 7 * (n - 1)
        F5 = n * v / (cv_f / 100 * gm) ** 2
        F10 = n * v / (cv_p / 100 * gm) ** 2
        row[f"F_6_{df2}_at_CV{cv_f:g}_n{n}"] = F5
        row[f"p_at_CV{cv_p:g}_n{n}"] = st.f.sf(F10, 6, df2)
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- Desirability / dominance (Eqs. 11-13)
def d_max(y, L, T, s=1.0):
    return np.clip((np.asarray(y, float) - L) / (T - L), 0, 1) ** s


def d_min(y, T, U, s=1.0):
    return np.clip((U - np.asarray(y, float)) / (U - T), 0, 1) ** s


def overall_D(d: np.ndarray, w: np.ndarray) -> float:
    d, w = np.asarray(d, float), np.asarray(w, float)
    if np.any(d == 0):
        return 0.0
    return float(np.exp((w * np.log(d)).sum() / w.sum()))


CRITERIA = [  # key, goal, weight (Table 8)
    ("compressive_MPa", "max", 5),
    ("fracture_toughness_MPa_m05", "max", 4),
    ("water_absorption_pct", "min", 4),
    ("thermal_conductivity_W_mK", "min", 3),
]


def pareto_dominance() -> pd.DataFrame:
    """For each mixture: does any other mixture match or beat it on ALL four criteria (>=, with one >)?"""
    df = D.load_all().set_index("mix_id")
    sign = {k: (1 if g == "max" else -1) for k, g, _ in CRITERIA}
    S = pd.DataFrame({k: sign[k] * df[k].astype(float) for k, _, _ in CRITERIA})
    rows = []
    for m in MIXES:
        dom_by = [o for o in MIXES if o != m and (S.loc[o] >= S.loc[m]).all() and (S.loc[o] > S.loc[m]).any()]
        rows.append({"mix_id": m, "dominated_by": ", ".join(dom_by) if dom_by else "none (non-dominated)"})
    return pd.DataFrame(rows)


def desirability_scenarios() -> pd.DataFrame:
    """Transparent recomputation of D under STATED limits (the original limits are unavailable).

    Scenario A: L/U = worst observed value, T = best observed value (min-max).
    Scenario B: limits widened by 25 % of the observed range beyond the worst value (so no d_i = 0).
    Only the four criteria that exist for all seven mixtures are used (CO2 emission exists for C0, AT3, AT6 only).
    These are NOT a reproduction of the reported D = 0.86.
    """
    df = D.load_all().set_index("mix_id")
    w = np.array([c[2] for c in CRITERIA], float)
    out = []
    for scen, pad in [("A: min-max limits", 0.0), ("B: limits widened by 25% of range", 0.25)]:
        for m in MIXES:
            d = []
            for key, goal, _ in CRITERIA:
                v = df[key].astype(float)
                lo, hi = v.min(), v.max()
                rng = hi - lo
                if goal == "max":
                    d.append(float(d_max(v[m], lo - pad * rng, hi)))
                else:
                    d.append(float(d_min(v[m], lo, hi + pad * rng)))
            out.append({"scenario": scen, "mix_id": m, "D": overall_D(np.array(d), w),
                        **{f"d_{k}": di for (k, _, _), di in zip(CRITERIA, d)}})
    res = pd.DataFrame(out)
    res["rank"] = res.groupby("scenario")["D"].rank(ascending=False, method="min").astype(int)
    return res
