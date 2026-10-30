# Data dictionary (`data/`)

| File | Columns |
|------|---------|
| mix_design.csv | mix_id; binder_system; lime_pct, metakaolin_pct, granite_waste_pct, sand_pct, basalt_fibre_pct (mass % of solids); w_b_ratio |
| physical_properties.csv | density_kg_m3; water_absorption_pct; apparent_porosity_pct; sorptivity_mm_min05 (28 d) |
| mechanical_properties.csv | compressive_MPa; flexural_MPa; splitting_tensile_MPa; fracture_toughness_MPa_m05 (28 d) |
| durability_residual_strength.csv | residual compressive strength (%), 28-d baseline: water 56 d, wet-dry 100 cyc, 5 % Na2SO4 56 d, 5 % H2SO4 28 d, freeze-thaw 50 cyc |
| thermal_properties.csv | thermal_conductivity_W_mK; thermal_diffusivity_mm2_s (as reported); specific_heat_J_kgK; cte_1e6_per_C |
| xrf_oxide_composition.csv | oxide mass % for ATM mortar, granite waste powder, metakaolin (Table 2) |
| lca_summary.csv | reported per-kg GWP, embodied and primary energy; waste utilisation (unverified) |
| templates/*.csv | empty record sheets (specimen, durability, lab log, LCA inventory, curves, XRD/XRF, PSD) |
| source/ | the original workbook and manuscript supplied |

All reported values are mixture-level summaries copied from the manuscript tables / workbook and verified identical between the two.
