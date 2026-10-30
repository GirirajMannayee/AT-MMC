"""Compare recomputed values against the numbers printed in the manuscript (Tables 10-23)."""
import numpy as np
import pandas as pd
from . import derived as X, stats as S
from .config import MIXES

# ---- values as printed in the manuscript ----------------------------------------------------
T10 = {("compressive_MPa", "splitting_tensile_MPa"): 0.998, ("water_absorption_pct", "apparent_porosity_pct"): 1.000,
       ("apparent_porosity_pct", "thermal_conductivity_W_mK"): 0.998, ("fracture_toughness_MPa_m05", "sorptivity_mm_min05"): -0.975,
       ("density_kg_m3", "compressive_MPa"): -0.255}
T11 = [[97.59, 100.00, 71.60, 7.69, 0.339], [99.37, 100.00, 88.55, 8.66, 0.032], [98.81, 100.00, 95.80, 2.01, 0.009]]
T13 = {"density_kg_m3": [-5.3, -4.8, -4.7, -4.1, -3.4, -2.9], "water_absorption_pct": [-20.2, -25.5, -32.9, -35.8, -38.5, -43.8],
       "apparent_porosity_pct": [-21.3, -27.3, -34.6, -37.3, -41.0, -45.4], "sorptivity_mm_min05": [-31.5, -43.5, -51.1, -55.4, -59.8, -66.3],
       "compressive_MPa": [8.3, 18.4, 26.7, 29.8, 36.3, 43.5], "flexural_MPa": [11.9, 19.6, 30.1, 34.2, 43.1, 51.6],
       "splitting_tensile_MPa": [12.8, 21.8, 33.3, 37.7, 46.7, 54.2], "fracture_toughness_MPa_m05": [29.5, 66.1, 85.7, 97.3, 104.5, 109.8],
       "acid_5pct_28d": [11.1, 15.7, 20.4, 24.4, 27.6, 30.6], "thermal_conductivity_W_mK": [-18.5, -27.2, -32.7, -35.2, -39.5, -42.6],
       "cte_1e6_per_C": [-8.3, -15.6, -20.8, -22.9, -25.0, -27.1]}
T14 = {"density_kg_m3": [1.4, 0.7, 0.5], "water_absorption_pct": [-19.5, -4.4, -8.6], "apparent_porosity_pct": [-20.3, -4.2, -7.5],
       "sorptivity_mm_min05": [-34.9, -8.9, -16.2], "compressive_MPa": [19.9, 2.5, 5.3], "flexural_MPa": [19.9, 3.1, 5.9],
       "splitting_tensile_MPa": [22.1, 3.3, 5.1], "fracture_toughness_MPa_m05": [52.4, 6.3, 2.6], "acid_5pct_28d": [12.0, 3.3, 2.4],
       "thermal_conductivity_W_mK": [-20.5, -3.7, -5.1], "cte_1e6_per_C": [-15.9, -2.6, -2.8]}
T15 = {("compressive_MPa", 7): (83.59, 5.66, 0.939), ("compressive_MPa", 6): (102.56, 7.90, 0.994),
       ("splitting_tensile_MPa", 7): (8.17, 6.79, 0.955), ("splitting_tensile_MPa", 6): (10.11, 9.14, 0.998),
       ("fracture_toughness_MPa_m05", 7): (6.45, 12.38, 0.964), ("fracture_toughness_MPa_m05", 6): (7.23, 13.62, 0.911)}
T16_ALPHA_REP_RHO = [0.753, 0.608, 0.527, 0.480, 0.456, 0.420, 0.394]
T16_ALPHA_IMPL = [0.790, 0.612, 0.539, 0.493, 0.471, 0.444, 0.414]
T16_RHOD = [2229, 2199, 2178, 2172, 2177, 2139, 2166]
T17 = {"C0": (None, None, None, 83.5), "AT3": (-35.5, -32.1, -32.5, 164.1), "AT6": (-46.1, -46.4, -46.0, 222.5)}
T18 = {"density_kg_m3": (18.4, 2215.7, .9342), "water_absorption_pct": (-.694, 5.103, .9772), "apparent_porosity_pct": (-1.584, 11.188, .9721),
       "sorptivity_mm_min05": (-.0146, .0612, .9559), "compressive_MPa": (5.62, 42.41, .9539), "flexural_MPa": (.944, 6.877, .9757),
       "splitting_tensile_MPa": (.554, 3.642, .9735), "fracture_toughness_MPa_m05": (.5, 1.525, .9396), "acid_5pct_28d": (5.56, 69.38, .9988)}
T19 = {"Mixture": (146.81, 30.70, .318), "Exposure": (442.63, 92.57, .640), "Mixture x exposure (residual)": (4.78, None, .041)}
T22 = {"density_kg_m3": (28.5, 13.5, -1.5), "water_absorption_pct": (-.435, -.265, -.075), "apparent_porosity_pct": (-1.03, -.51, -.12),
       "sorptivity_mm_min05": (-.009, -.005, -.001), "compressive_MPa": (4.5, 2, .8), "flexural_MPa": (.935, .385, .135),
       "splitting_tensile_MPa": (.48, .19, .05), "fracture_toughness_MPa_m05": (.175, .095, -.035), "acid_5pct_28d": (4.2, 2.2, -.3)}
T23 = {"density_kg_m3": (1.9, 2.7, .4), "water_absorption_pct": (20.6, 28.7, 48.3), "apparent_porosity_pct": (21.9, 30.5, 54.7),
       "sorptivity_mm_min05": (38.5, 53.7, 169.0), "compressive_MPa": (12.7, 17.7, 18.5), "flexural_MPa": (14.4, 20.0, 23.5),
       "splitting_tensile_MPa": (15.0, 20.9, 25.7), "fracture_toughness_MPa_m05": (24.4, 34.0, 67.7), "acid_5pct_28d": (9.1, 12.7, 9.5)}
TEXT = {"Eq3 Ca/(Si+Al) AT1-AT4": (1.58, 0.005), "Eq3 Ca/(Si+Al) AT5-AT6": (1.31, 0.005), "R2 porosity vs acid residual (Fig 18)": (0.979, 0.0005),
        "Kendall W": (0.989, 0.0005), "AT6 GWP change vs C0 (%)": (-46.1, 0.05), "AT6 embodied energy change (%)": (-46.4, 0.05),
        "AT6 primary energy change (%)": (-46.0, 0.05)}


def compare() -> pd.DataFrame:
    rows = []

    def add(item, ms, mine, tol):
        ok = (ms is None) or (abs(ms - mine) <= tol)
        rows.append({"item": item, "manuscript": ms, "recomputed": mine, "abs_diff": (abs(ms - mine) if ms is not None else np.nan),
                     "tolerance": tol, "status": "match" if ok else "DIFFERS"})

    t10 = X.table10().set_index(["x", "y"])["pearson_r"]
    for k, v in T10.items(): add(f"Table 10 r {k[0]} vs {k[1]}", v, t10[k], 0.0006)
    t11 = X.oxide_consistency().iloc[:, 1:].values
    cols = ["sum_ox", "sum_ox+LOI", "SiO2+Al2O3+Fe2O3", "SiO2/Al2O3", "CaO/SiO2"]
    for i, m in enumerate(["ATM", "GWP", "MK"]):
        for j, c in enumerate(cols): add(f"Table 11 {m} {c}", T11[i][j], t11[i][j], 0.006)
    ch = X.pct_change_vs_c0()
    for k, vals in T13.items():
        for j, v in enumerate(vals): add(f"Table 13 {k} AT{j + 1}", v, ch.loc[k, f"AT{j + 1}"], 0.051)
    fc = X.fibre_series_changes()
    for k, vals in T14.items():
        for j, v in enumerate(vals): add(f"Table 14 {k} col{j + 1}", v, fc.loc[k].iloc[j], 0.051)
    ry = S.ryshkewitch_fits().set_index(["response", "n"])
    for (k, n), (y0, b, r2) in T15.items():
        r = ry.loc[(k, n)]
        add(f"Table 15 {k} n={n} y0", y0, r.y0_fitting_constant, 0.006); add(f"Table 15 {k} n={n} b", b, r.b_Eq5_fraction, 0.006)
        add(f"Table 15 {k} n={n} R2", r2, r.R2_of_ln_y, 0.0006)
    pt = X.physical_thermal_consistency()
    for i, m in enumerate(MIXES):
        add(f"Table 16 {m} rho_d", T16_RHOD[i], pt.iloc[i]["rho_d_implied_kg_m3"], 0.6)
        add(f"Table 16 {m} alpha calc (reported rho)", T16_ALPHA_REP_RHO[i], pt.iloc[i]["alpha_calc_reported_rho_mm2_s"], 0.0006)
        add(f"Table 16 {m} alpha calc (implied rho_d)", T16_ALPHA_IMPL[i], pt.iloc[i]["alpha_calc_implied_rho_mm2_s"], 0.0006)
    lc = X.lca_relative()
    for m, (g, e, p, eco) in T17.items():
        for nm, ms, col in [("dGWP", g, "dGWP_vs_C0_pct"), ("dEE", e, "dEmbodiedEnergy_vs_C0_pct"), ("dPE", p, "dPrimaryEnergy_vs_C0_pct")]:
            if ms is not None: add(f"Table 17 {m} {nm}", ms, lc.loc[m, col], 0.051)
        add(f"Table 17 {m} E_eco", eco, lc.loc[m, "E_eco_MPa_per_kgCO2eq_per_kg"], 0.051)
    f18 = S.fibre_dose_fits().set_index("property")
    for k, (sl, ic, r2) in T18.items():
        add(f"Table 18 {k} slope", sl, f18.loc[k, "slope_per_1pct_BF"], abs(sl) * 0.006 + 1e-4)
        add(f"Table 18 {k} intercept", ic, f18.loc[k, "intercept_0pct_BF"], abs(ic) * 0.0006 + 1e-4)
        add(f"Table 18 {k} R2", r2, f18.loc[k, "R2_descriptive"], 0.00006)
    tw = S.two_way_unreplicated().set_index("source")
    for k, (ms_, F, eta) in T19.items():
        add(f"Table 19 {k} MS", ms_, tw.loc[k, "MS"], 0.006)
        if F is not None: add(f"Table 19 {k} F", F, tw.loc[k, "F_nominal"], 0.006)
        add(f"Table 19 {k} eta2", eta, tw.loc[k, "eta2"], 0.0006)
    ct = S.contrasts_2x2().set_index("property")
    for k, v in T22.items():
        for j, c in enumerate(ct.columns): add(f"Table 22 {k} {c}", v[j], ct.loc[k, c], 0.0006)
    rp = S.replicate_precision().set_index("property")
    for k, (c3, c5, f5) in T23.items():
        add(f"Table 23 {k} maxCV n=3", c3, rp.loc[k, "max_CV_pct_n3"], 0.051); add(f"Table 23 {k} maxCV n=5", c5, rp.loc[k, "max_CV_pct_n5"], 0.051)
        add(f"Table 23 {k} F(6,14) at CV5", f5, rp.loc[k, "F_6_14_at_CV5_n3"], 0.051)
    bi = X.binder_indices().set_index("Mix")
    add("Eq3 Ca/(Si+Al) AT1-AT4", 1.58, bi.loc["AT1", "Ca/(Si+Al) molar, Eq.(3)"], 0.005)
    add("Eq3 Ca/(Si+Al) AT5-AT6", 1.31, bi.loc["AT5", "Ca/(Si+Al) molar, Eq.(3)"], 0.005)
    add("R2 porosity vs acid residual (Fig 18)", 0.979, S.porosity_acid_fit()["R2"], 0.0005)
    add("Kendall W", 0.989, S.kendall_w_and_friedman()["kendall_W"], 0.0005)
    return pd.DataFrame(rows)
