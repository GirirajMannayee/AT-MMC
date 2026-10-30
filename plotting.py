"""All figures. Clean plots: no watermark, no logo, no decorative overlay.
Each figure is saved as PNG (300 dpi) and vector PDF in results/figures/."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from . import data as D, derived as X, stats as S
from .config import FIGURES, MIXES, PROPS, EXPOSURES

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.labelsize": 9, "axes.titlesize": 10,
                     "figure.dpi": 100, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "legend.frameon": False})
COL = {"C0": "#5b5b5b", "AT1": "#9ecae1", "AT2": "#6baed6", "AT3": "#3182bd", "AT4": "#2ca25f",
       "AT5": "#e6851e", "AT6": "#c0392b"}
BLUE, RED, GREY = "#1f4e79", "#c0392b", "#7f8c8d"
PAL = [COL[m] for m in MIXES]


def _save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}")
    plt.close(fig)
    return name


def _bars(ax, df, key, ylabel, fmt="{:g}", highlight="AT6"):
    v = df[key].astype(float).values
    b = ax.bar(MIXES, v, color=PAL, edgecolor="black", linewidth=0.5)
    for r, x in zip(b, v):
        ax.text(r.get_x() + r.get_width() / 2, x, fmt.format(x), ha="center", va="bottom", fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, v.max() * 1.15)


def _panel_figure(keys, df, name, titles=None, ncols=2, size=(7.2, 5.6)):
    n = len(keys)
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=size)
    axs = np.atleast_1d(axs).ravel()
    for i, k in enumerate(keys):
        lab, unit = PROPS[k]
        _bars(axs[i], df, k, f"{lab} ({unit})")
        axs[i].set_title(f"({'abcdefgh'[i]}) {lab}", loc="left")
    for j in range(n, len(axs)):
        axs[j].axis("off")
    fig.tight_layout()
    return _save(fig, name)


def fig03_xrf():
    ox = D.load_oxides().set_index("oxide")
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), sharey=True)
    for ax, (c, t, col) in zip(axs, [("ATM_mortar", "(a) Ancient temple mortar", "#1f77b4"),
                                     ("granite_waste_powder", "(b) Granite waste powder", "#2ca02c"),
                                     ("metakaolin", "(c) Metakaolin", "#7b3fa0")]):
        v = ox[c]
        bars = ax.bar(v.index, v.values, color=col)
        for r, x in zip(bars, v.values):
            ax.text(r.get_x() + r.get_width() / 2, x, f"{x:.2f}", ha="center", va="bottom", fontsize=6.5)
        ax.set_title(t, loc="left"); ax.tick_params(axis="x", rotation=45)
    axs[0].set_ylabel("Mass %")
    fig.tight_layout()
    return _save(fig, "fig03_xrf_oxide_composition")


def fig06_physical():
    return _panel_figure(["density_kg_m3", "water_absorption_pct", "apparent_porosity_pct", "sorptivity_mm_min05"],
                         D.load_all(), "fig06_physical_properties")


def fig07_mechanical():
    return _panel_figure(["compressive_MPa", "flexural_MPa", "splitting_tensile_MPa", "fracture_toughness_MPa_m05"],
                         D.load_all(), "fig07_mechanical_properties")


def fig10_durability():
    d = D.load_durability().set_index("mix_id")
    fig, ax = plt.subplots(figsize=(8.5, 4))
    w = 0.12
    x = np.arange(len(EXPOSURES))
    for i, m in enumerate(MIXES):
        ax.bar(x + (i - 3) * w, d.loc[m, list(EXPOSURES)].values.astype(float), w, label=m, color=COL[m],
               edgecolor="black", linewidth=0.4)
    ax.set_xticks(x); ax.set_xticklabels(list(EXPOSURES.values()), rotation=15)
    ax.set_ylabel("Residual compressive strength (%)"); ax.set_ylim(50, 102)
    ax.legend(ncol=7, loc="upper center", bbox_to_anchor=(0.5, 1.12))
    return _save(fig, "fig10_residual_strength_durability")


def fig12_thermal():
    df = D.load_all()
    df["alpha_reported_mm2_s"] = df["alpha_reported_mm2_s"]
    keys = ["thermal_conductivity_W_mK", "alpha_reported_mm2_s", "specific_heat_J_kgK", "cte_1e6_per_C"]
    labs = ["Thermal conductivity (W/m·K)", "Thermal diffusivity, reported (mm²/s)", "Specific heat (J/kg·K)", "CTE (10$^{-6}$/°C)"]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6))
    for i, (ax, k, l) in enumerate(zip(axs.ravel(), keys, labs)):
        _bars(ax, df, k, l)
        ax.set_title(f"({'abcd'[i]})", loc="left")
    fig.tight_layout()
    return _save(fig, "fig12_thermal_properties")


def fig14_descriptive():
    df = D.load_all().set_index("mix_id")
    fits = S.fibre_dose_fits().set_index("property")
    fig, axs = plt.subplots(2, 2, figsize=(9, 7))
    # (a)
    ax = axs[0, 0]; ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
    x = df.loc[["AT1", "AT2", "AT3", "AT4"], "basalt_fibre_pct"].astype(float)
    for a, k, c, mk in [(ax, "compressive_MPa", BLUE, "o"), (ax2, "fracture_toughness_MPa_m05", RED, "s")]:
        y = df.loc[x.index, k].astype(float)
        a.plot(x, y, mk, color=c)
        xx = np.array([0, 1.5]); a.plot(xx, fits.loc[k, "intercept_0pct_BF"] + fits.loc[k, "slope_per_1pct_BF"] * xx, "-", color=c, lw=1)
        a.tick_params(axis="y", colors=c)
    ax.set_ylabel("Compressive strength (MPa)", color=BLUE); ax2.set_ylabel("Fracture toughness (MPa·m$^{0.5}$)", color=RED)
    ax.set_xlabel("Basalt fibre (mass %), AT1–AT4")
    ax.set_title("(a) Descriptive linear trend, AT1–AT4 means", loc="left")
    ax.text(0.03, 0.95, f"R² (descriptive): {fits.loc['compressive_MPa','R2_descriptive']:.3f} / {fits.loc['fracture_toughness_MPa_m05','R2_descriptive']:.3f}",
            transform=ax.transAxes, fontsize=7, va="top")
    # (b)
    ax = axs[0, 1]
    d = D.load_durability().set_index("mix_id")
    for m in MIXES:
        ax.plot(range(5), d.loc[m, list(EXPOSURES)].values.astype(float), "-o", ms=3.5, color=COL[m], label=m)
    ax.set_xticks(range(5)); ax.set_xticklabels(["Water", "Wet–dry", "Sulfate", "Acid", "Freeze–thaw"], rotation=20)
    ax.set_ylabel("Residual strength (%)"); ax.set_title("(b) Mixture × exposure profiles", loc="left"); ax.legend(ncol=2, fontsize=7)
    # (c)
    ax = axs[1, 0]
    tw = S.two_way_unreplicated()
    shares = tw.loc[:2, "eta2"].values * 100
    bars = ax.bar(["Mixture\n(df 6)", "Exposure\n(df 4)", "Mix × exposure\n(df 24)"], shares, color=[BLUE, GREY, RED], width=0.6)
    for r, v in zip(bars, shares):
        ax.text(r.get_x() + r.get_width() / 2, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=8)
    kw = S.kendall_w_and_friedman()
    ax.text(0.03, 0.95, f"Kendall's W = {kw['kendall_W']:.3f} (descriptive)", transform=ax.transAxes, fontsize=7.5, va="top")
    ax.set_ylabel("Share of total SS (%)"); ax.set_ylim(0, 85)
    ax.set_title("(c) Unreplicated two-way decomposition", loc="left")
    # (d)
    ax = axs[1, 1]
    rp = S.replicate_precision()
    names = {"density_kg_m3": "Density", "water_absorption_pct": "Water abs.", "apparent_porosity_pct": "App. porosity",
             "sorptivity_mm_min05": "Sorptivity", "compressive_MPa": "Compressive", "flexural_MPa": "Flexural",
             "splitting_tensile_MPa": "Splitting tensile", "fracture_toughness_MPa_m05": "Fracture toughness",
             "acid_5pct_28d": "Acid residual"}
    rp = rp.sort_values("max_CV_pct_n5")
    y = np.arange(len(rp))
    ax.barh(y + 0.2, rp["max_CV_pct_n5"], 0.4, color=GREY, label="n = 5")
    ax.barh(y - 0.2, rp["max_CV_pct_n3"], 0.4, color=BLUE, label="n = 3")
    ax.set_yticks(y); ax.set_yticklabels([names[p] for p in rp["property"]])
    ax.set_xlabel("Hypothetical max within-mixture CV (%) for nominal p < 0.05")
    ax.set_title("(d) Replicate-precision planning scenarios", loc="left"); ax.legend(loc="lower right")
    fig.tight_layout()
    return _save(fig, "fig14_descriptive_summaries")


def fig15_compressive():
    df = D.load_all()
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.plot(MIXES, df["compressive_MPa"], "-o", color=BLUE)
    for m, v in zip(MIXES, df["compressive_MPa"]):
        ax.annotate(f"{v:g}", (m, v), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=8)
    ax.set_ylabel("28-day compressive strength (MPa)"); ax.set_ylim(35, 60)
    return _save(fig, "fig15_compressive_strength_by_mix")


def fig16_lca():
    l = D.load_lca()
    items = [("gwp_kgCO2eq_per_kg", "GWP (kg CO$_2$-eq/kg)"), ("embodied_energy_MJ_per_kg", "Embodied energy (MJ/kg)"),
             ("primary_energy_MJ_per_kg", "Primary energy (MJ/kg)"), ("waste_utilization_pct_unverified", "Waste utilization (%)\n(as reported; unverified)")]
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.3))
    for ax, (k, lab) in zip(axs, items):
        b = ax.bar(l["mix_id"], l[k], color=[COL[m] for m in l["mix_id"]], edgecolor="black", linewidth=0.5)
        for r, v in zip(b, l[k]):
            ax.text(r.get_x() + r.get_width() / 2, v, f"{v:g}", ha="center", va="bottom", fontsize=8)
        ax.set_ylabel(lab); ax.set_ylim(0, l[k].max() * 1.15)
        if k.startswith("waste"):
            for r in b: r.set_hatch("//")
    fig.tight_layout()
    return _save(fig, "fig16_life_cycle_impacts")


def fig17_improvement():
    ch = X.pct_change_vs_c0()["AT6"].drop("density_kg_m3")
    names = {k: PROPS[k][0] for k in ch.index}
    ch = ch.sort_values()
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.barh([names[k] for k in ch.index], ch.values, color=[BLUE if v > 0 else RED for v in ch.values])
    for i, v in enumerate(ch.values):
        ax.text(v + (2 if v > 0 else -2), i, f"{v:+.1f}", va="center", ha="left" if v > 0 else "right", fontsize=8)
    ax.axvline(0, color="black", lw=0.8); ax.set_xlabel("Change of AT6 relative to C0 (%), Eq. (14)")
    ax.set_xlim(-85, 130)
    return _save(fig, "fig17_AT6_change_vs_C0")


def fig18_porosity_acid():
    df = D.load_all(); f = S.porosity_acid_fit()
    fig, ax = plt.subplots(figsize=(5, 4))
    for m, p, a in zip(MIXES, df["apparent_porosity_pct"], df["acid_5pct_28d"]):
        ax.scatter(p, a, color=COL[m], edgecolor="black", zorder=3); ax.annotate(m, (p, a), xytext=(4, 4), textcoords="offset points", fontsize=8)
    xx = np.linspace(7.5, 14.6, 20); ax.plot(xx, f["intercept"] + f["slope"] * xx, "--", color=GREY)
    ax.set_xlabel("Apparent porosity (%)"); ax.set_ylabel("Residual strength after 5% H$_2$SO$_4$ (%)")
    ax.text(0.97, 0.95, f"R² = {f['R2']:.3f} (n = 7, descriptive)", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    return _save(fig, "fig18_porosity_vs_acid_residual")


def fig19_correlation():
    c = X.correlation_matrix()
    labels = [PROPS[k][0].replace("$_2$", "2").replace("$_4$", "4") for k in c.columns]
    fig, ax = plt.subplots(figsize=(8, 6.8))
    im = ax.imshow(c.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{c.values[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="white" if abs(c.values[i, j]) > 0.6 else "black")
    ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, label="Pearson r (n = 7 mixture summaries)", shrink=0.8)
    return _save(fig, "fig19_correlation_matrix")


def fig20_heatmap():
    ch = X.pct_change_vs_c0().drop(columns="C0")
    labels = [PROPS[k][0].replace("$_2$", "2").replace("$_4$", "4") for k in ch.index]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    lim = np.abs(ch.values).max()
    im = ax.imshow(ch.values, cmap="PuOr_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(6)); ax.set_xticklabels(ch.columns); ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    for i in range(ch.shape[0]):
        for j in range(ch.shape[1]):
            ax.text(j, i, f"{ch.values[i, j]:+.1f}", ha="center", va="center", fontsize=7)
    ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, label="Change vs C0 (%)")
    return _save(fig, "fig20_heatmap_change_vs_C0")


def fig21_fibre_dose():
    df = D.load_all().set_index("mix_id")
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.6))
    for ax, k, lab, t in [(axs[0], "fracture_toughness_MPa_m05", "Fracture toughness (MPa·m$^{0.5}$)", "(a)"),
                          (axs[1], "compressive_MPa", "Compressive strength (MPa)", "(b)")]:
        s1 = ["AT1", "AT2", "AT3", "AT4"]; s2 = ["AT5", "AT6"]
        ax.plot(df.loc[s1, "basalt_fibre_pct"], df.loc[s1, k], "-o", color=BLUE, label="Binder 20/10/20, w/b 0.42")
        ax.plot(df.loc[s2, "basalt_fibre_pct"], df.loc[s2, k], "-s", color=RED, label="Binder 25/15/25, w/b 0.40")
        ax.axhline(df.loc["C0", k], ls="--", color=GREY, label="C0")
        ax.set_xlabel("Basalt fibre (mass %)"); ax.set_ylabel(lab); ax.set_title(t, loc="left")
    axs[0].legend(fontsize=7)
    fig.tight_layout()
    return _save(fig, "fig21_fibre_dosage_series")


def fig22_ryshkewitch():
    df = D.load_all().set_index("mix_id"); fits = S.ryshkewitch_fits()
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 3.8))
    for ax, k, lab, t in [(axs[0], "compressive_MPa", "Compressive strength (MPa)", "(a)"),
                          (axs[1], "fracture_toughness_MPa_m05", "Fracture toughness (MPa·m$^{0.5}$)", "(b)")]:
        P = df["apparent_porosity_pct"].astype(float)
        for m in MIXES:
            ax.scatter(P[m], df.loc[m, k], color=COL[m], edgecolor="black", zorder=3)
        r = fits[(fits.response == k) & (fits.n == 7)].iloc[0]
        xx = np.linspace(7.5, 14.6, 50); ax.plot(xx, np.exp(r.intercept_a + r.slope_b_per_pct * xx), "--", color=GREY,
                                                  label=f"y = {r.y0_fitting_constant:.2f}·exp(−{r.b_Eq5_fraction:.2f}P), R² = {r.R2_of_ln_y:.3f}")
        ax.set_xlabel("Apparent porosity (%)"); ax.set_ylabel(lab); ax.set_title(t, loc="left"); ax.legend(fontsize=7)
    fig.tight_layout()
    return _save(fig, "fig22_ryshkewitch_fits")


def fig23_diffusivity():
    t = X.physical_thermal_consistency()
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    ax.scatter(t["alpha_calc_reported_rho_mm2_s"], t["alpha_reported_mm2_s"], marker="o", color=BLUE, label="Calc. with reported ρ")
    ax.scatter(t["alpha_calc_implied_rho_mm2_s"], t["alpha_reported_mm2_s"], marker="s", color=RED, label="Calc. with implied ρ$_d$")
    ax.plot([0.35, 0.82], [0.35, 0.82], "--", color="black", lw=0.8, label="1:1")
    ax.set_xlabel("Computed α = k/(ρ·c$_p$) (mm²/s)"); ax.set_ylabel("Reported α (mm²/s)"); ax.legend(fontsize=7)
    return _save(fig, "fig23_diffusivity_reported_vs_computed")


def fig24_schematic():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    ax = axs[0]; ax.axis("off"); ax.set_title("(a) Proposed matrix pathway (hypothesis only)", loc="left")
    boxes = ["Lime +\nmetakaolin", "Ca(OH)$_2$ + AS$_2$\n+ H$_2$O", "C–(A)–S–H +\nCa-aluminate\nhydrates", "Carbonation\n→ CaCO$_3$"]
    for i, b in enumerate(boxes):
        ax.add_patch(Rectangle((i * 0.25 + 0.01, 0.35), 0.2, 0.3, fill=False, lw=1)); ax.text(i * 0.25 + 0.11, 0.5, b, ha="center", va="center", fontsize=7)
        if i < 3: ax.add_patch(FancyArrowPatch((i * 0.25 + 0.21, 0.5), (i * 0.25 + 0.26, 0.5), arrowstyle="->", mutation_scale=10))
    ax.text(0.5, 0.15, "Reaction products not quantified in this study", ha="center", fontsize=7, style="italic"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax = axs[1]; ax.set_title("(b) Fibre bridging across a crack (hypothesis only)", loc="left"); ax.axis("off")
    ax.add_patch(Rectangle((0, 0.3), 0.45, 0.4, color="#d9d9d9")); ax.add_patch(Rectangle((0.55, 0.3), 0.45, 0.4, color="#d9d9d9"))
    for y in (0.4, 0.5, 0.6): ax.plot([0.3, 0.7], [y, y + 0.02], color=RED, lw=1.5)
    ax.text(0.5, 0.15, "crack", ha="center", fontsize=8); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.85, "Illustration only; no measured quantities", ha="center", fontsize=7, style="italic")
    return _save(fig, "fig24_conceptual_schematic")


def figS1_mix_composition():
    m = D.load_mix().set_index("mix_id").fillna(0)
    cols = ["lime_pct", "metakaolin_pct", "granite_waste_pct", "sand_pct", "basalt_fibre_pct"]
    labs = ["Lime", "Metakaolin", "Granite waste", "Sand", "Basalt fibre"]
    fig, ax = plt.subplots(figsize=(7, 3.8)); bottom = np.zeros(len(m))
    for c, l, col in zip(cols, labs, ["#c7e9b4", "#7fcdbb", "#bdbdbd", "#e6c88a", RED]):
        ax.bar(m.index.astype(str), m[c], bottom=bottom, label=l, color=col, edgecolor="black", linewidth=0.4); bottom += m[c].values
    ax.set_ylabel("Mass % of solids"); ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.13), fontsize=7)
    return _save(fig, "figS1_mix_composition")


def figS2_desirability():
    sc = S.desirability_scenarios()
    fig, ax = plt.subplots(figsize=(6.5, 3.8)); w = 0.38
    for i, (s, g) in enumerate(sc.groupby("scenario")):
        ax.bar(np.arange(7) + (i - 0.5) * w, g["D"].values, w, label=s, color=[BLUE, GREY][i])
    ax.set_xticks(range(7)); ax.set_xticklabels(MIXES); ax.set_ylabel("Overall desirability D (4 criteria)"); ax.legend(fontsize=7)
    return _save(fig, "figS2_desirability_limit_sensitivity")


ALL = [fig03_xrf, fig06_physical, fig07_mechanical, fig10_durability, fig12_thermal, fig14_descriptive,
       fig15_compressive, fig16_lca, fig17_improvement, fig18_porosity_acid, fig19_correlation, fig20_heatmap,
       fig21_fibre_dose, fig22_ryshkewitch, fig23_diffusivity, fig24_schematic, figS1_mix_composition, figS2_desirability]


def make_all():
    return [f() for f in ALL]
