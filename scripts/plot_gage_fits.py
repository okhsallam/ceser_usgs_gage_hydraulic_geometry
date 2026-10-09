#!/usr/bin/env python3
"""Plot the hydraulic-geometry power-law fits for one or more gages.

Each gage gets a four-panel card: the three fitted relations (width, depth,
velocity vs discharge) on log-log axes, plus the stage-discharge field
measurements with the derived riverbed stage where bathymetry exists.

    python3 scripts/plot_gage_fits.py                      # 4 representative gages
    python3 scripts/plot_gage_fits.py --sites 05411850 05586100
    python3 scripts/plot_gage_fits.py --best 6             # 6 highest-R2 gages
    python3 scripts/plot_gage_fits.py --sites 05411850 --out figures/one.png

Needs only numpy, pandas, matplotlib.
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FT_PER_M = 3.28084

# (panel title, observation column, coefficient col, exponent col, R2 col, y label)
RELATIONS = [
    ("Width",    "channel_width",    "coeff_a_width",    "exp_b_width",    "r2_width",    "Width  $w$  [ft]",    "w = a·Q^b"),
    ("Depth",    "_depth",           "coeff_c_depth",    "exp_f_depth",    "r2_depth",    "Mean depth  $d$  [ft]", "d = c·Q^f"),
    ("Velocity", "channel_velocity", "coeff_k_velocity", "exp_m_velocity", "r2_velocity", "Velocity  $v$  [ft/s]", "v = k·Q^m"),
]


def load():
    coef = pd.read_csv(os.path.join(DATA, "hydraulic_geometry_coefficients.csv"),
                       dtype={"site_id": str, "huc4": str}).set_index("site_id")
    bath_path = os.path.join(DATA, "channel_bathymetry.csv")
    bath = pd.read_csv(bath_path, dtype={"site_id": str, "huc4": str}).set_index("site_id")
    return coef, bath


def gage_observations(site):
    """Channel measurements used for the three power-law fits."""
    p = os.path.join(DATA, "channel_features", f"{site}.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p, dtype={"site_id": str})
    for c in ("channel_flow", "channel_width", "channel_area", "channel_velocity"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["_depth"] = df["channel_area"] / df["channel_width"]
    return df[(df.channel_flow > 0) & df.channel_width.notna() & df.channel_area.notna()]


def gage_stage_discharge(site):
    """Stage-discharge field measurements used for the bathymetry step."""
    p = os.path.join(DATA, "field_measurements", f"{site}.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p, dtype={"site_id": str})
    for c in ("MeanGageHeight", "Discharge"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["MeanGageHeight", "Discharge"])


def draw_gage(axes, site, coef, bath):
    row = coef.loc[site]
    obs = gage_observations(site)

    # --- three power-law panels ---------------------------------------------
    for ax, (_name, ycol, ck, ek, r2k, ylab, form) in zip(axes[:3], RELATIONS):
        ax.set_xlabel("Discharge  $Q$  [ft³/s]")
        ax.set_ylabel(ylab)
        ax.set_xscale("log")
        ax.set_yscale("log")

        c, e = row[ck], row[ek]
        y = obs[ycol] if obs is not None else pd.Series(dtype=float)
        good = y.notna() & (y > 0) if len(y) else y
        if len(y) and good.any():
            ax.scatter(obs.loc[good, "channel_flow"], y[good], s=16, c=_style.BLUE,
                       alpha=0.55, edgecolors=_style.SURFACE, linewidths=0.4, zorder=3)
        if pd.notna(c) and pd.notna(e) and len(y) and good.any():
            q = np.logspace(np.log10(obs.loc[good, "channel_flow"].min()),
                            np.log10(obs.loc[good, "channel_flow"].max()), 120)
            ax.plot(q, c * q ** e, color=_style.ORANGE, lw=2.0, zorder=4)
            lhs = form.split("=")[0].strip()
            r2 = row[r2k]
            ax.annotate(f"${lhs}$ = {c:.3g}·Q$^{{{e:.3f}}}$", xy=(0.04, 0.95),
                        xycoords="axes fraction", va="top", fontsize=9,
                        color=_style.INK_2)
            ax.annotate("R² = %s" % (f"{r2:.2f}" if pd.notna(r2) else "n/a"),
                        xy=(0.04, 0.845), xycoords="axes fraction", va="top",
                        fontsize=9.5, fontweight="bold",
                        color=_style.INK if pd.notna(r2) and r2 >= 0.7 else _style.MUTED)
            ax.annotate(f"n = {int(good.sum())}", xy=(0.04, 0.745),
                        xycoords="axes fraction", va="top", fontsize=8.5,
                        color=_style.MUTED)
        else:
            ax.annotate("fit unavailable", xy=(0.5, 0.5), xycoords="axes fraction",
                        ha="center", color=_style.MUTED, fontsize=9)
        ax.grid(True, which="both", alpha=0.45)
        ax.set_axisbelow(True)

    # --- stage-discharge panel ----------------------------------------------
    ax = axes[3]
    sd = gage_stage_discharge(site)
    ax.set_xlabel("Discharge  $Q$  [ft³/s]")
    ax.set_ylabel("Mean gage height  [ft]")
    if sd is not None and len(sd):
        ax.scatter(sd.Discharge, sd.MeanGageHeight, s=16, c=_style.BLUE, alpha=0.55,
                   edgecolors=_style.SURFACE, linewidths=0.4, zorder=3)
        ax.set_xscale("log")
        if site in bath.index:
            b = bath.loc[site]
            stage0 = (b.riverbed_elev_m - b.gage_datum_alt_va_m) * FT_PER_M
            ax.axhline(stage0, color=_style.ORANGE, lw=2.0, ls="--", zorder=4)
            lo = min(stage0, sd.MeanGageHeight.min())
            hi = max(stage0, sd.MeanGageHeight.max())
            pad = 0.08 * (hi - lo or 1.0)
            ax.set_ylim(lo - pad, hi + 3.2 * pad)
            ax.annotate(f"riverbed {stage0:.2f} ft\nburn {b.burn_value_m_30m:.2f} m\nn = {len(sd)}",
                        xy=(0.04, 0.95), xycoords="axes fraction", va="top",
                        fontsize=8.5, color=_style.INK_2, linespacing=1.5)
    else:
        ax.annotate("no stage–discharge\nmeasurements", xy=(0.5, 0.5),
                    xycoords="axes fraction", ha="center", va="center",
                    color=_style.MUTED, fontsize=9, linespacing=1.5)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        for sp in ax.spines.values():
            sp.set_color(_style.GRID)
        ax.grid(False)
        return _header(row, site)
    ax.grid(True, which="both", alpha=0.45)
    ax.set_axisbelow(True)

    return _header(row, site)


def _header(row, site):
    flag = "" if row.exponent_sum_in_range == 1 else "  ⚠ outside 0.85–1.15"
    return (f"USGS {site} — {row.station_name}\n"
            f"HUC4 {row.huc4} · stream order {int(row.stream_order)} · "
            f"{row.drain_area_sqkm:,.0f} km² · b+f+m = {row.exponent_sum:.3f}{flag}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", nargs="*", default=None, help="USGS site ids")
    ap.add_argument("--best", type=int, default=None,
                    help="instead, plot the N gages with the highest mean R²")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    _style.apply()
    coef, bath = load()

    if args.best:
        mean_r2 = coef[["r2_width", "r2_depth", "r2_velocity"]].mean(axis=1)
        sites = list(mean_r2.sort_values(ascending=False).head(args.best).index)
    elif args.sites:
        sites = [s.zfill(8) for s in args.sites]
    else:
        # default: a spread of stream orders among well-fit, bathymetry-bearing gages
        mean_r2 = coef[["r2_width", "r2_depth", "r2_velocity"]].mean(axis=1)
        pool = coef[(coef.index.isin(bath.index)) & (mean_r2 > 0.7)].copy()
        pool["_r2"] = mean_r2[pool.index]
        sites = [g.sort_values("_r2", ascending=False).index[0]
                 for _, g in pool.groupby("stream_order")]
        sites = sites[:4]

    missing = [s for s in sites if s not in coef.index]
    if missing:
        sys.exit(f"no coefficients for: {', '.join(missing)}")

    n = len(sites)
    titles = []
    fig = plt.figure(figsize=(15.0, 3.6 * n + 0.55), constrained_layout=True)
    gs = fig.add_gridspec(n + 1, 4, height_ratios=[1] * n + [0.1])
    axes = np.array([[fig.add_subplot(gs[r, c]) for c in range(4)] for r in range(n)])
    legend_ax = fig.add_subplot(gs[n, :])
    legend_ax.axis("off")
    for i, site in enumerate(sites):
        titles.append(draw_gage(axes[i], site, coef, bath))
    fig.suptitle("Hydraulic-geometry power-law fits  ·  CESER domain",
                 fontsize=13, fontweight="bold")
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, alpha=0.7,
               markerfacecolor=_style.BLUE, markeredgecolor=_style.SURFACE,
               label="Observed USGS measurement"),
        Line2D([], [], color=_style.ORANGE, lw=2.0, label="Fitted power law"),
        Line2D([], [], color=_style.ORANGE, lw=2.0, ls="--",
               label="Derived riverbed stage"),
    ]
    legend_ax.legend(handles=handles, loc="center", ncol=3, frameon=False,
                     handletextpad=0.7, columnspacing=3.0)
    fig.get_layout_engine().set(h_pad=0.30, w_pad=0.06, hspace=0.14, wspace=0.03)

    # Row headers go on the FIGURE, added after the layout settles: a wide text
    # artist parented to an axes would force constrained_layout to shrink every
    # panel to fit it.
    fig.canvas.draw()
    fig.set_layout_engine("none")
    for i, title in enumerate(titles):
        pos = axes[i][0].get_position()
        top = max(a.get_position().y1 for a in axes[i])
        fig.text(pos.x0, top + 0.012, title, fontsize=10.5, fontweight="bold",
                 color=_style.INK, ha="left", va="bottom", linespacing=1.45)

    out = args.out or os.path.join(ROOT, "figures", "gage_fits.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}  ({', '.join(sites)})")


if __name__ == "__main__":
    main()
