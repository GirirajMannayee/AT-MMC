# AT-MMC: data, code, tables, figures and analysis package

Supporting repository for *"Experimental Development and Exploratory Evaluation of Ancient Temple-Inspired Mineral Matrix Composites for Sustainable Engineering Applications"* (G. Mannayee, VIT Vellore).
Intended location: https://github.com/GirirajMannayee/AT-MMC

## What is (and is not) in the data
The supplied records contain **one summary value per mixture** (C0, AT1-AT6) and per mixture x exposure cell. There are **no specimen-level data, replicate counts or standard deviations**. Consequently:

* Everything computed here from the reported tables is **descriptive**. F, p and R2 values are nominal and are *not* evidence of experimental significance (manuscript Sections 3.9 and 4.14).
* Replicate-based inference (one-/two-way ANOVA, Welch, Kruskal-Wallis, Tukey, Dunnett, effect sizes, power) is implemented in `at_mmc_anova.py` and runs as soon as specimen-level data are entered. On the blank templates it deliberately computes nothing.
* No data were invented. Figures that need curve/pattern data (XRD, flexural and load-CMOD curves, wet-dry mass, 200-800 C strength) are not regenerated; templates are provided (`data/templates/`).

## Quick start
```bash
pip install -r requirements.txt
python scripts/run_all.py        # all tables + figures + audit
python -m pytest -q tests        # 8 tests
python at_mmc_anova.py           # replicate ANOVA (needs specimen-level data)
```

## Layout
```
data/                 reported data (CSV), source workbook + manuscript, blank templates
src/at_mmc/           data.py  derived.py  stats.py  plotting.py  reproduction.py
scripts/run_all.py    one-command reproduction
at_mmc_anova.py       replicate-based ANOVA suite (name as cited in the manuscript)
tests/                reproduction + ANOVA self-tests (synthetic data used only inside tests)
results/tables/       20 CSV tables;  results/AT_MMC_results_tables.xlsx (same, one workbook)
results/figures/      18 figures, PNG 300 dpi + PDF, no watermark
results/original_manuscript_images/   images embedded in the submitted .docx
docs/                 DATA_ISSUES.md  FIGURE_MAP.md  DATA_DICTIONARY.md
```

## Main results (reported summaries, descriptive)
* **Audit:** 260 of 260 numeric values printed in manuscript Tables 10-23 and the quoted text values are reproduced from the data (`results/tables/S3_reproduction_check.csv`).
* AT6 vs C0: compressive +43.5 %, flexural +51.6 %, splitting tensile +54.2 %, fracture toughness +109.8 %; water absorption -43.8 %, porosity -45.4 %, sorptivity -66.3 %.
* AT6 is Pareto non-dominated and every other mixture is dominated on compressive strength, fracture toughness, water absorption and thermal conductivity, so its rank-1 position does not depend on desirability limits. The numerical D depends on limits (D = 1.00 under min-max limits, 1.00 with widened limits; the reported 0.86 cannot be reproduced without the original limits).
* Mixture x exposure decomposition (Table 5, one value per cell): mixture 31.8 %, exposure 64.0 %, residual 4.1 % of total SS; Kendall W = 0.989.
* AT1-AT4 means vs fibre dose: R2 = 0.934-0.999 (descriptive; sand changes with fibre dose).
* Thermal check: reported diffusivity exceeds k/(rho cp) by 3.7-9.5 %; implied dry density is 0.7-5.3 % below listed density.
* Eco-indicator fc/GWP: 83.5 (C0), 164.1 (AT3), 222.5 (AT6), per-kg basis only.

Open items and provenance gaps: see `docs/DATA_ISSUES.md`. Figure-by-figure status: `docs/FIGURE_MAP.md`.

## License
Code MIT; data/figures suggested CC BY 4.0 (author to confirm). See `LICENSE`, `CITATION.cff`.
