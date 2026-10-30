"""Loading and validation of the reported (mixture-level) data."""
import pandas as pd
from .config import DATA, MIXES


def _read(name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA / name)
    df["mix_id"] = pd.Categorical(df["mix_id"], categories=MIXES, ordered=True)
    return df.sort_values("mix_id").reset_index(drop=True)


def load_mix():
    return _read("mix_design.csv")


def load_physical():
    return _read("physical_properties.csv")


def load_mechanical():
    return _read("mechanical_properties.csv")


def load_durability():
    return _read("durability_residual_strength.csv")


def load_thermal():
    return _read("thermal_properties.csv")


def load_oxides():
    return pd.read_csv(DATA / "xrf_oxide_composition.csv")


def load_lca():
    return pd.read_csv(DATA / "lca_summary.csv")


def load_all() -> pd.DataFrame:
    """One wide table, one row per mixture, all reported 28-day summaries."""
    df = load_mix()
    for part in (load_physical(), load_mechanical(), load_durability(), load_thermal()):
        df = df.merge(part.rename(columns={"thermal_diffusivity_mm2_s": "alpha_reported_mm2_s"}),
                      on="mix_id", suffixes=("", "_dup"))
    df["mix_id"] = pd.Categorical(df["mix_id"], categories=MIXES, ordered=True)
    return df.sort_values("mix_id").reset_index(drop=True)


def mass_balance(mix: pd.DataFrame) -> pd.DataFrame:
    """Sum of listed solid fractions per mixture (Section 3.2 mass-balance statement)."""
    cols = ["lime_pct", "metakaolin_pct", "granite_waste_pct", "sand_pct", "basalt_fibre_pct"]
    out = mix[["mix_id"]].copy()
    out["sum_listed_solids_pct"] = mix[cols].fillna(0).sum(axis=1)
    out["closes_to_100"] = (out["sum_listed_solids_pct"] - 100).abs() < 1e-9
    out["note"] = ""
    out.loc[out["mix_id"] == "C0", "note"] = "OPC mass fraction not stated; C0 listed as sand = 100 % (unresolved, see docs/DATA_ISSUES.md)"
    return out
