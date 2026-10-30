"""Paths, labels and constants shared by all modules."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
TEMPLATES = DATA / "templates"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"

MIXES = ["C0", "AT1", "AT2", "AT3", "AT4", "AT5", "AT6"]
AT_MIXES = MIXES[1:]

# Molar masses (g/mol) used in Eq. (3), Table 11 and Table 12
M_CAOH2 = 74.093
M_SIO2 = 60.084
M_AL2O3 = 101.961
M_CAO = 56.077
RHO_WATER = 1000.0  # kg/m3, Eq. (4)

# Short column keys -> (label, unit)
PROPS = {
    "density_kg_m3": ("Density", "kg/m³"),
    "water_absorption_pct": ("Water absorption", "%"),
    "apparent_porosity_pct": ("Apparent porosity", "%"),
    "sorptivity_mm_min05": ("Sorptivity", "mm/min$^{0.5}$"),
    "compressive_MPa": ("Compressive strength", "MPa"),
    "flexural_MPa": ("Flexural strength", "MPa"),
    "splitting_tensile_MPa": ("Splitting tensile strength", "MPa"),
    "fracture_toughness_MPa_m05": ("Fracture toughness", "MPa·m$^{0.5}$"),
    "acid_5pct_28d": ("Residual strength, 5% H$_2$SO$_4$", "%"),
    "thermal_conductivity_W_mK": ("Thermal conductivity", "W/m·K"),
    "cte_1e6_per_C": ("CTE", "10$^{-6}$/°C"),
}
EXPOSURES = {
    "water_immersion_56d": "Water immersion (56 d)",
    "wet_dry_100cyc": "Wet–dry (100 cyc)",
    "sulfate_5pct_56d": "Sulfate 5% (56 d)",
    "acid_5pct_28d": "Acid 5% (28 d)",
    "freeze_thaw_50cyc": "Freeze–thaw (50 cyc)",
}
