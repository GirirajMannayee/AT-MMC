"""Derived, arithmetic quantities (Tables 11-14, 16, 17, 20; Eqs. 3, 4, 6, 14, 15, 16).

Everything here is computed from the reported mixture-level summaries only.
"""
import numpy as np
import pandas as pd
from . import data as D
from .config import (M_AL2O3, M_CAOH2, M_SIO2, RHO_WATER, MIXES, PROPS)


# ------------------------------------------------------------------ Table 11
def oxide_consistency() -> pd.DataFrame:
    ox = D.load_oxides().set_index("oxide")
    rows = []
    for mat, label in [("ATM_mortar", "Ancient temple mortar (ATM)"),
                       ("granite_waste_powder", "Granite waste powder (GWP)"),
                       ("metakaolin", "Metakaolin (MK)")]:
        c = ox[mat]
        s_ox = c.drop("LOI").sum()
        rows.append({
            "Material": label,
            "Sum oxides excl. LOI": s_ox,
            "Sum oxides + LOI": s_ox + c["LOI"],
            "SiO2+Al2O3+Fe2O3": c["SiO2"] + c["Al2O3"] + c["Fe2O3"],
            "SiO2/Al2O3 (molar)": (c["SiO2"] / M_SIO2) / (c["Al2O3"] / M_AL2O3),
            "CaO/SiO2 (mass)": c["CaO"] / c["SiO2"],
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Table 12 / Eq. (3)
def binder_indices() -> pd.DataFrame:
    mix = D.load_mix().set_index("mix_id")
    mk = D.load_oxides().set_index("oxide")["metakaolin"]
    bal = D.mass_balance(D.load_mix()).set_index("mix_id")
    rows = []
    for m in MIXES:
        r = mix.loc[m]
        if m == "C0" or pd.isna(r.lime_pct):
            rows.append({"Mix": m, "Lime+MK (%)": np.nan, "Lime:MK": np.nan, "GW:(Lime+MK)": np.nan,
                         "Ca/(Si+Al) molar, Eq.(3)": np.nan,
                         "Sum listed solids (%)": bal.loc[m, "sum_listed_solids_pct"]})
            continue
        mL, mMK = r.lime_pct, r.metakaolin_pct
        denom = mMK * mk["SiO2"] / 100 / M_SIO2 + 2 * mMK * mk["Al2O3"] / 100 / M_AL2O3
        rows.append({"Mix": m, "Lime+MK (%)": mL + mMK, "Lime:MK": mL / mMK,
                     "GW:(Lime+MK)": r.granite_waste_pct / (mL + mMK),
                     "Ca/(Si+Al) molar, Eq.(3)": (mL / M_CAOH2) / denom,
                     "Sum listed solids (%)": bal.loc[m, "sum_listed_solids_pct"]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ Eq. (14), Tables 13 / 14 / 20
def pct_change_vs_c0() -> pd.DataFrame:
    df = D.load_all().set_index("mix_id")
    out = {}
    for key, (label, unit) in PROPS.items():
        base = df.loc["C0", key]
        out[key] = 100 * (df[key] - base) / base
    res = pd.DataFrame(out).T
    res.index.name = "property"
    return res[MIXES]


def fibre_series_changes() -> pd.DataFrame:
    """Table 14: AT1->AT4, AT3->AT4, AT5->AT6."""
    df = D.load_all().set_index("mix_id")
    pairs = {"AT1->AT4 (0->1.5% BF; w/b 0.42)": ("AT1", "AT4"),
             "AT3->AT4 (1.0->1.5% BF; w/b 0.42)": ("AT3", "AT4"),
             "AT5->AT6 (1.0->1.5% BF; w/b 0.40)": ("AT5", "AT6")}
    out = {}
    for key in PROPS:
        out[key] = {name: 100 * (df.loc[b, key] - df.loc[a, key]) / df.loc[a, key]
                    for name, (a, b) in pairs.items()}
    res = pd.DataFrame(out).T
    res.index.name = "property"
    return res


def mechanical_change_table() -> pd.DataFrame:
    """Table 20."""
    ch = pct_change_vs_c0()
    keys = ["compressive_MPa", "flexural_MPa", "splitting_tensile_MPa", "fracture_toughness_MPa_m05"]
    t = ch.loc[keys].T
    t["mean_of_four"] = t.mean(axis=1)
    t.index.name = "mix_id"
    return t


# ------------------------------------------------------------------ Eqs. (4), (6), Table 16
def physical_thermal_consistency() -> pd.DataFrame:
    df = D.load_all().set_index("mix_id")
    rho_d = RHO_WATER * df["apparent_porosity_pct"] / df["water_absorption_pct"]  # kg/m3, Eq. (4)
    k, cp = df["thermal_conductivity_W_mK"], df["specific_heat_J_kgK"]
    out = pd.DataFrame({
        "rho_reported_kg_m3": df["density_kg_m3"],
        "rho_d_implied_kg_m3": rho_d,
        "rho_d_vs_reported_pct": 100 * (rho_d - df["density_kg_m3"]) / df["density_kg_m3"],
        "alpha_reported_mm2_s": df["alpha_reported_mm2_s"],
        "alpha_calc_reported_rho_mm2_s": k / (df["density_kg_m3"] * cp) * 1e6,
        "alpha_calc_implied_rho_mm2_s": k / (rho_d * cp) * 1e6,
    })
    out["alpha_reported_minus_calc_reported_rho"] = out["alpha_reported_mm2_s"] - out["alpha_calc_reported_rho_mm2_s"]
    out["alpha_excess_pct_over_calc_reported_rho"] = 100 * out["alpha_reported_minus_calc_reported_rho"] / out["alpha_calc_reported_rho_mm2_s"]
    return out


# ------------------------------------------------------------------ Table 17 / Eq. (16)
def lca_relative() -> pd.DataFrame:
    lca = D.load_lca().set_index("mix_id")
    mech = D.load_mechanical().set_index("mix_id")
    base = lca.loc["C0"]
    out = pd.DataFrame(index=lca.index)
    out["GWP_kgCO2eq_per_kg"] = lca["gwp_kgCO2eq_per_kg"]
    for col, name in [("gwp_kgCO2eq_per_kg", "dGWP_vs_C0_pct"), ("embodied_energy_MJ_per_kg", "dEmbodiedEnergy_vs_C0_pct"),
                      ("primary_energy_MJ_per_kg", "dPrimaryEnergy_vs_C0_pct")]:
        out[name] = 100 * (lca[col] - base[col]) / base[col]
    out["fc_MPa"] = mech.loc[lca.index, "compressive_MPa"]
    out["E_eco_MPa_per_kgCO2eq_per_kg"] = out["fc_MPa"] / out["GWP_kgCO2eq_per_kg"]
    out["E_eco_ratio_to_C0"] = out["E_eco_MPa_per_kgCO2eq_per_kg"] / out.loc["C0", "E_eco_MPa_per_kgCO2eq_per_kg"]
    return out


# ------------------------------------------------------------------ Eq. (15), Table 10, Fig. 19
def correlation_matrix(method: str = "pearson") -> pd.DataFrame:
    df = D.load_all()
    cols = list(PROPS.keys())
    return df[cols].corr(method=method)


TABLE10_PAIRS = [
    ("compressive_MPa", "splitting_tensile_MPa"),
    ("water_absorption_pct", "apparent_porosity_pct"),
    ("apparent_porosity_pct", "thermal_conductivity_W_mK"),
    ("fracture_toughness_MPa_m05", "sorptivity_mm_min05"),
    ("density_kg_m3", "compressive_MPa"),
]


def table10() -> pd.DataFrame:
    c = correlation_matrix()
    return pd.DataFrame([{"x": a, "y": b, "pearson_r": c.loc[a, b],
                          "spearman_rho": correlation_matrix("spearman").loc[a, b]} for a, b in TABLE10_PAIRS])
