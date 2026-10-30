import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from at_mmc import data as D, reproduction as R, derived as X, stats as S


def test_every_manuscript_value_reproduced():
    c = R.compare()
    bad = c[c.status != "match"]
    assert bad.empty, bad.to_string()


def test_mass_balance_all_100():
    mb = D.mass_balance(D.load_mix())
    assert mb["closes_to_100"].all()


def test_oxide_totals_close_with_loi():
    t = X.oxide_consistency()
    assert (abs(t["Sum oxides + LOI"] - 100) < 0.01).all()


def test_at6_pareto_non_dominated():
    p = S.pareto_dominance().set_index("mix_id")["dominated_by"]
    assert p["AT6"].startswith("none")
    assert all(not p[m].startswith("none") for m in ["C0", "AT1", "AT2", "AT3", "AT4", "AT5"])


def test_no_watermark_text_in_plot_code():
    src = (Path(__file__).resolve().parents[1] / "src/at_mmc/plotting.py").read_text().lower()
    assert "watermark" not in src.replace("no watermark", "")
