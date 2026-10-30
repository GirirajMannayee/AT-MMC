# Data issues carried over from the manuscript and workbook

These are items the manuscript itself flags, plus what the recomputation in this repository shows.
Nothing here has been "fixed" by inventing data.

| # | Issue | Evidence from this package | Needed to close it |
|---|-------|----------------------------|--------------------|
| 1 | No specimen-level data, replicate counts or SDs | Only one value per mixture / per mixture x exposure cell | Fill `data/templates/specimen_raw_data_template.csv`, run `python at_mmc_anova.py` |
| 2 | Non-orthogonal mixture design | AT1-AT4 vary fibre and sand together; AT5-AT6 also change binder, filler and w/b (`T12`, `T22`) | Designed mixture experiment |
| 3 | C0 composition | C0 listed as sand = 100 %, OPC mass fraction not stated (`T03_mass_balance`) | Original batching sheet |
| 4 | Thermal diffusivity | Reported alpha exceeds k/(rho cp) by 3.7-9.5 % with reported density (`T16`) | Instrument outputs, density/moisture state |
| 5 | Water absorption vs porosity | Implied dry density (1000 P/WA) is 0.7-5.3 % below listed density (`T16`) | Original mass records |
| 6 | Porosity vs thermal conductivity r = +0.998 | Opposite to usual pore-volume effect (`T10`) | Moisture state, units, composition-controlled specimens |
| 7 | Desirability D = 0.86 not reproducible | Limits, targets, exponents unavailable; CO2 data exist only for C0, AT3, AT6. AT6 is non-dominated on all four criteria with complete data (`S1`), so its rank is limit-independent, but D depends on limits (`S2`) | Original limits/targets/model file |
| 8 | 62 % waste utilisation vs 25 % granite waste | Cannot be reconciled from Table 1 | Numerator/denominator definition |
| 9 | LCA | Inventory sheet is blank; per-kg values only; no functional unit | Populated inventory (`lca_inventory_template.csv`) |
| 10 | Curve/pattern data absent | Figs 4, 8, 9, 11, 13 (XRD, flexural curves, load-CMOD, wet-dry mass, high-temperature residual strength) have no numeric data | Fill the curve templates in `data/templates/`, then add plotting functions |
| 11 | Image duplication | In the .docx the same embedded image (`image5.jpeg`) is referenced twice in the Fig. 5 area | Check SEM panel (a)/(b) source files |
| 12 | Provenance | Sample ID, date, location, permission, SEM/XRD acquisition metadata not recorded | `lab_record_log_template.csv` |
| 13 | Rounding | Table 14, AT3->AT4 fracture toughness is exactly +6.25 %; manuscript prints 6.3, Python rounds to 6.2 | none (cosmetic) |
