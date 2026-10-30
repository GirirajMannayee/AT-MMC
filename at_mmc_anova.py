#!/usr/bin/env python3
"""Replicate-based ANOVA for AT-MMC (the script referred to in manuscript Section 3.9).

Input : specimen-level CSV in the layout of data/templates/specimen_raw_data_template.csv
        (needs columns 'Mix ID', 'Batch ID', 'Test', 'Calculated result', 'Valid? (Y/N)')
        optional durability CSV in the layout of durability_raw_data_template.csv
        (needs 'Mix ID', 'Exposure type', 'Retention (%) formula' holding the numeric retention %).
Output: tables in --outdir: for every Test (and every exposure) one-way ANOVA, Welch ANOVA,
        Kruskal-Wallis, Shapiro-Wilk (residuals), Levene, effect sizes (eta2, omega2, Cohen f),
        achieved power, Tukey HSD (all pairs) and Dunnett (vs C0); two-way ANOVA mixture x exposure
        with replicates for the durability file.

It refuses to run on mixture means: each mixture needs >= 2 valid specimens.
With the blank template it exits cleanly and states that no data are available.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st

CONTROL = "C0"


def one_way(groups: dict, alpha=0.05) -> dict:
    names = list(groups)
    arrs = [np.asarray(groups[k], float) for k in names]
    k, N = len(arrs), sum(len(a) for a in arrs)
    gm = np.concatenate(arrs).mean()
    ss_b = sum(len(a) * (a.mean() - gm) ** 2 for a in arrs)
    ss_w = sum(((a - a.mean()) ** 2).sum() for a in arrs)
    df_b, df_w = k - 1, N - k
    ms_b, ms_w = ss_b / df_b, ss_w / df_w
    F = ms_b / ms_w
    p = st.f.sf(F, df_b, df_w)
    eta2 = ss_b / (ss_b + ss_w)
    omega2 = max(0.0, (ss_b - df_b * ms_w) / (ss_b + ss_w + ms_w))
    f = np.sqrt(eta2 / (1 - eta2)) if eta2 < 1 else np.inf
    power = st.ncf.sf(st.f.ppf(1 - alpha, df_b, df_w), df_b, df_w, f ** 2 * N)  # observed (achieved) power
    # Welch ANOVA
    n_i = np.array([len(a) for a in arrs], float)
    v_i = np.array([a.var(ddof=1) for a in arrs])
    if np.all(v_i > 0):
        w = n_i / v_i
        mw = (w * np.array([a.mean() for a in arrs])).sum() / w.sum()
        lam = ((1 - w / w.sum()) ** 2 / (n_i - 1)).sum()
        Fw = ((w * (np.array([a.mean() for a in arrs]) - mw) ** 2).sum() / (k - 1)) / (1 + 2 * (k - 2) / (k ** 2 - 1) * lam)
        dfw2 = (k ** 2 - 1) / (3 * lam)
        pw = st.f.sf(Fw, k - 1, dfw2)
    else:
        Fw = dfw2 = pw = np.nan
    kw = st.kruskal(*arrs)
    resid = np.concatenate([a - a.mean() for a in arrs])
    sh = st.shapiro(resid) if len(resid) >= 3 else None
    lev = st.levene(*arrs, center="median")
    return {"k": k, "N": N, "F": F, "df_between": df_b, "df_within": df_w, "p": p, "eta2": eta2, "omega2": omega2,
            "cohen_f": f, "achieved_power": power, "Welch_F": Fw, "Welch_df2": dfw2, "Welch_p": pw,
            "KruskalWallis_H": kw.statistic, "KruskalWallis_p": kw.pvalue,
            "Shapiro_W": getattr(sh, "statistic", np.nan), "Shapiro_p": getattr(sh, "pvalue", np.nan),
            "Levene_W": lev.statistic, "Levene_p": lev.pvalue}


def posthoc(groups: dict) -> tuple:
    names = list(groups)
    arrs = [np.asarray(groups[k], float) for k in names]
    tk = st.tukey_hsd(*arrs)
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            rows.append({"group_1": names[i], "group_2": names[j], "mean_diff": arrs[i].mean() - arrs[j].mean(),
                         "p_Tukey": tk.pvalue[i, j]})
    dun = pd.DataFrame()
    if CONTROL in groups:
        others = [n for n in names if n != CONTROL]
        dr = st.dunnett(*[np.asarray(groups[n], float) for n in others], control=np.asarray(groups[CONTROL], float))
        dun = pd.DataFrame({"mix": others, "mean_diff_vs_C0": [np.mean(groups[n]) - np.mean(groups[CONTROL]) for n in others],
                            "p_Dunnett": dr.pvalue})
    return pd.DataFrame(rows), dun


def two_way_replicated(df: pd.DataFrame, value="value", a="Mix ID", b="Exposure type") -> pd.DataFrame:
    """Balanced or unbalanced two-way ANOVA with interaction (Type II via statsmodels)."""
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    d = df.rename(columns={a: "A", b: "B", value: "y"})
    m = smf.ols("y ~ C(A) * C(B)", data=d).fit()
    t = sm.stats.anova_lm(m, typ=2)
    t["eta2_partial"] = t["sum_sq"] / (t["sum_sq"] + t.loc["Residual", "sum_sq"])
    return t.reset_index().rename(columns={"index": "source"})


def run(specimen: Path, durability: Path | None, outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    wrote = 0
    sp = pd.read_csv(specimen) if specimen and specimen.exists() else pd.DataFrame()
    if not sp.empty:
        if "Valid? (Y/N)" in sp:
            sp = sp[sp["Valid? (Y/N)"].astype(str).str.upper().str.startswith("Y")]
        sp["Calculated result"] = pd.to_numeric(sp["Calculated result"], errors="coerce")
        sp = sp.dropna(subset=["Calculated result"])
        summ, ph_all, dn_all = [], [], []
        for test, g in sp.groupby("Test"):
            groups = {m: v["Calculated result"].values for m, v in g.groupby("Mix ID") if len(v) >= 2}
            if len(groups) < 2:
                print(f"[skip] {test}: need >= 2 mixtures with >= 2 valid specimens each")
                continue
            res = one_way(groups); res["Test"] = test
            summ.append(res)
            tk, dn = posthoc(groups)
            tk["Test"] = test; dn["Test"] = test
            ph_all.append(tk); dn_all.append(dn)
            desc = g.groupby("Mix ID")["Calculated result"].agg(["count", "mean", "std"])
            desc["cv_pct"] = 100 * desc["std"] / desc["mean"]
            desc.to_csv(outdir / f"descriptives_{str(test).replace(' ', '_')}.csv")
        if summ:
            pd.DataFrame(summ).to_csv(outdir / "oneway_anova_summary.csv", index=False)
            pd.concat(ph_all).to_csv(outdir / "tukey_hsd.csv", index=False)
            pd.concat(dn_all).to_csv(outdir / "dunnett_vs_C0.csv", index=False)
            wrote += 1
    if durability and durability.exists():
        du = pd.read_csv(durability)
        if not du.empty and "Retention (%) formula" in du:
            du["value"] = pd.to_numeric(du["Retention (%) formula"], errors="coerce")
            du = du.dropna(subset=["value"])
            if len(du):
                two_way_replicated(du).to_csv(outdir / "twoway_mixture_x_exposure.csv", index=False)
                wrote += 1
    if wrote == 0:
        msg = ("No specimen-level data found (the templates are blank). No inferential statistics were computed, "
               "by design. Enter the original individual measurements and re-run.")
        (outdir / "NO_DATA.txt").write_text(msg + "\n")
        print(msg)
    else:
        print(f"Replicate-based results written to {outdir}")
    return 0


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--specimen", type=Path, default=root / "data/templates/specimen_raw_data_template.csv")
    ap.add_argument("--durability", type=Path, default=root / "data/templates/durability_raw_data_template.csv")
    ap.add_argument("--outdir", type=Path, default=root / "results/replicate_anova")
    a = ap.parse_args()
    sys.exit(run(a.specimen, a.durability, a.outdir))
