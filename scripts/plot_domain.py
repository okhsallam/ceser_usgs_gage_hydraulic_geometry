#!/usr/bin/env python3
"""Plot the CESER domain: HUC4 boundaries, order-5+ streams, and all 88 gages
coloured by what data each one carries.

    python3 scripts/plot_domain.py                 # -> figures/domain_overview.png
    python3 scripts/plot_domain.py --out my.png

Needs only numpy, pandas, matplotlib.
"""
import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _style

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def _rings(geom):
    """Yield coordinate rings from a (Multi)Polygon geometry dict."""
    if geom["type"] == "Polygon":
        yield geom["coordinates"]
    elif geom["type"] == "MultiPolygon":
        for poly in geom["coordinates"]:
            yield poly


def polygon_patch(geom, **kw):
    """Build one compound PathPatch per (Multi)Polygon, holes included."""
    verts, codes = [], []
    for poly in _rings(geom):
        for ring in poly:
            ring = np.asarray(ring)[:, :2]
            if len(ring) < 3:
                continue
            verts.extend(ring)
            codes.extend([Path.MOVETO] + [Path.LINETO] * (len(ring) - 1))
    if not verts:
        return None
    return PathPatch(Path(np.asarray(verts), codes), **kw)


def ring_centroid(geom):
    """Area-weighted centroid of the largest exterior ring (shoelace)."""
    best, best_area = None, -1.0
    for poly in _rings(geom):
        r = np.asarray(poly[0])[:, :2]
        if len(r) < 3:
            continue
        x, y = r[:, 0], r[:, 1]
        cross = x[:-1] * y[1:] - x[1:] * y[:-1]
        a = cross.sum() / 2.0
        if abs(a) < 1e-12:
            continue
        cx = ((x[:-1] + x[1:]) * cross).sum() / (6 * a)
        cy = ((y[:-1] + y[1:]) * cross).sum() / (6 * a)
        if abs(a) > best_area:
            best, best_area = (cx, cy), abs(a)
    return best


def line_segments(geom):
    """Yield coordinate arrays from a (Multi)LineString geometry dict."""
    if geom["type"] == "LineString":
        yield np.asarray(geom["coordinates"])[:, :2]
    elif geom["type"] == "MultiLineString":
        for part in geom["coordinates"]:
            yield np.asarray(part)[:, :2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "figures", "domain_overview.png"))
    args = ap.parse_args()

    _style.apply()
    gages = pd.read_csv(os.path.join(DATA, "gages.csv"), dtype={"site_id": str, "huc4": str})

    fig, ax = plt.subplots(figsize=(9.5, 9.0))

    # --- basemap: streams, then HUC4 boundaries (both recessive) -------------
    with open(os.path.join(DATA, "streams.geojson")) as fh:
        streams = json.load(fh)["features"]
    for f in streams:
        for seg in line_segments(f["geometry"]):
            ax.plot(seg[:, 0], seg[:, 1], color=_style.BLUE_200,
                    lw=0.5, zorder=1, solid_capstyle="round")

    with open(os.path.join(DATA, "domain_boundary.geojson")) as fh:
        bnd = json.load(fh)["features"]
    for f in bnd:
        p = polygon_patch(f["geometry"], facecolor="none",
                          edgecolor=_style.AXIS, lw=1.1, zorder=2)
        if p:
            ax.add_patch(p)
        c = ring_centroid(f["geometry"])
        if c:
            ax.annotate(f["properties"]["huc4"], c, color=_style.MUTED,
                        fontsize=13, fontweight="bold", ha="center", va="center",
                        alpha=0.55, zorder=4,
                        path_effects=[pe.withStroke(linewidth=3.5, foreground=_style.SURFACE)])

    # --- gages, split into the three data tiers ------------------------------
    tiers = [
        ("Coefficients + bathymetry", (gages.has_coefficients == 1) & (gages.has_bathymetry == 1), _style.BLUE),
        ("Coefficients only",         (gages.has_coefficients == 1) & (gages.has_bathymetry == 0), _style.ORANGE),
        ("No fitted coefficients",    gages.has_coefficients == 0,                                 _style.AQUA),
    ]
    # area-proportional marker size from drainage area
    da = gages.drain_area_sqkm.fillna(gages.drain_area_sqkm.median())
    size = 18 + 190 * np.sqrt(da / da.max())

    for label, mask, color in tiers:
        ax.scatter(gages.loc[mask, "lng"], gages.loc[mask, "lat"],
                   s=size[mask], c=color, edgecolors=_style.SURFACE,
                   linewidths=1.0, zorder=5, label=f"{label} (n={int(mask.sum())})")

    # --- frame ---------------------------------------------------------------
    lat0 = np.deg2rad(gages.lat.mean())
    ax.set_aspect(1.0 / np.cos(lat0))  # keep degrees visually square
    ax.set_xlabel("Longitude [°E]")
    ax.set_ylabel("Latitude [°N]")
    ax.set_title("CESER domain — USGS gages with hydraulic-geometry parameters", pad=12)
    ax.grid(True, alpha=0.5)
    ax.set_axisbelow(True)

    legend = ax.legend(loc="lower left", title="Data available per gage",
                       labelspacing=0.9, borderpad=0.8)
    legend.get_title().set_color(_style.INK_2)
    legend.get_title().set_fontsize(9)

    # size key, kept separate from the colour legend
    ref = [500, 5000, 20000]
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=np.sqrt(18 + 190 * np.sqrt(v / da.max())),
                      markerfacecolor=_style.MUTED, markeredgecolor=_style.SURFACE,
                      label=f"{v:,} km²") for v in ref]
    ax2 = ax.legend(handles=handles, loc="lower right", title="Drainage area",
                    labelspacing=1.3, borderpad=0.9, handletextpad=1.2)
    ax2.get_title().set_color(_style.INK_2)
    ax2.get_title().set_fontsize(9)
    ax.add_artist(legend)

    ax.annotate(f"{len(gages)} gages · HUC4 0706 · 0708 · 0709 · 0712 · 0713\n"
                f"Order-5+ NHDPlus flowlines shown in blue",
                xy=(0.5, -0.085), xycoords="axes fraction",
                ha="center", va="top", color=_style.MUTED, fontsize=8.5)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
