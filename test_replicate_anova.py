"""Self-test of at_mmc_anova.py on SYNTHETIC data generated here. Synthetic values are never written to data/ or results/."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import at_mmc_anova as A


def synth(tmp_path):
    rng = np.random.default_rng(1)
    rows = []
    for mix, mu in zip(["C0", "AT1", "AT2", "AT3"], [38, 42, 46, 49]):
        for i in range(5):
            rows.append({"Batch ID": f"B{i}", "Mix ID": mix, "Specimen ID": f"{mix}-{i}", "Test": "Compressive",
                         "Calculated result": rng.normal(mu, 1.0), "Valid? (Y/N)": "Y"})
    p = tmp_path / "s.csv"; pd.DataFrame(rows).to_csv(p, index=False); return p


def test_runs_and_detects_difference(tmp_path):
    p = synth(tmp_path)
    A.run(p, None, tmp_path / "out")
    s = pd.read_csv(tmp_path / "out/oneway_anova_summary.csv")
    assert s.loc[0, "p"] < 1e-6 and s.loc[0, "df_within"] == 16
    assert (tmp_path / "out/tukey_hsd.csv").exists() and (tmp_path / "out/dunnett_vs_C0.csv").exists()


def test_blank_template_gives_no_inference(tmp_path):
    A.run(ROOT / "data/templates/specimen_raw_data_template.csv", None, tmp_path / "o")
    assert (tmp_path / "o/NO_DATA.txt").exists()


def test_one_way_matches_scipy():
    from scipy import stats
    g = {"a": [1, 2, 3, 4.], "b": [2, 3, 4, 6.], "c": [5, 6, 7, 9.]}
    r = A.one_way(g)
    f = stats.f_oneway(*g.values())
    assert abs(r["F"] - f.statistic) < 1e-9 and abs(r["p"] - f.pvalue) < 1e-9
