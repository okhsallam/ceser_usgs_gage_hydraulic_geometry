#!/usr/bin/env python3
"""
Build a domain-wide (CESER), HUC-agnostic data package of USGS gage
hydraulic-geometry parameters and the field measurements they were fit from.

Reads the per-HUC4 dataset tree and writes a flat, shareable package.
Re-runnable: output is overwritten deterministically.
"""
import os, re, csv, json, glob, shutil, collections

import numpy as np

SRC = "/lcrc/project/hydrosm/osallam/flood_ai/Illinois_HiPIMS/data/huc_ceser/dataset"
PKG = "/lcrc/project/hydrosm/osallam/flood_ai/ceser_usgs_gage_hydraulic_geometry"
OUT = os.path.join(PKG, "data")
HUCS = ["huc_0706", "huc_0708", "huc_0709", "huc_0712", "huc_0713"]
CRS84 = {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}


def sid(x):
    """Normalize any USGS id spelling to an 8-digit zero-padded string."""
    s = re.sub(r"[^0-9]", "", str(x))
    return s.zfill(8)


def num(x):
    if x is None or x == "":
        return None
    try:
        f = float(x)
        return int(f) if f.is_integer() else f
    except (TypeError, ValueError):
        return x


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  {os.path.relpath(path, PKG):<52} {len(rows):>5} rows")


def write_geojson(path, rows, lon="lng", lat="lat", name=None):
    feats = []
    for r in rows:
        if r.get(lon) in (None, "") or r.get(lat) in (None, ""):
            continue
        props = {k: v for k, v in r.items()}
        feats.append({
            "type": "Feature",
            "properties": props,
            "geometry": {"type": "Point", "coordinates": [float(r[lon]), float(r[lat])]},
        })
    fc = {"type": "FeatureCollection",
          "name": name or os.path.splitext(os.path.basename(path))[0],
          "crs": CRS84, "features": feats}
    with open(path, "w") as fh:
        json.dump(fc, fh, indent=1)
    print(f"  {os.path.relpath(path, PKG):<52} {len(feats):>5} features")


# ---------------------------------------------------------------------------
# 1. Read sources, deduplicating across (and within) HUC4 folders
# ---------------------------------------------------------------------------
gage_meta, gage_hucs = {}, collections.defaultdict(set)
for h in HUCS:
    h4 = h.replace("huc_", "")
    with open(f"{SRC}/{h}/geo/selected_gages_on_streams.csv") as fh:
        for r in csv.DictReader(fh):
            s = sid(r["staid"])
            gage_hucs[s].add(h4)
            gage_meta.setdefault(s, r)

coeffs = {}
for h in HUCS:
    h4 = h.replace("huc_", "")
    with open(f"{SRC}/{h}/geo/hydraulic_geometry_parameters.csv") as fh:
        for r in csv.DictReader(fh):
            s = sid(r["site_id"])
            gage_hucs[s].add(h4)
            coeffs[s] = r

bath = {}  # site -> {res: props}
for h in HUCS:
    for res in (10, 30):
        p = f"{SRC}/{h}/geo/gagues_field_measurements/bathymetry_parameter_summary_{res}.geojson"
        for f in json.load(open(p))["features"]:
            bath.setdefault(sid(f["properties"]["site_id"]), {})[res] = f["properties"]

fm_src, cf_src, peaks = {}, {}, set()
for h in HUCS:
    g = f"{SRC}/{h}/geo/gagues_field_measurements"
    for p in glob.glob(g + "/*_field_measurements.csv"):
        fm_src[sid(os.path.basename(p).split("_field")[0])] = p
    for p in glob.glob(g + "/*_channel_features_rich.csv"):
        cf_src[sid(os.path.basename(p).split("_channel")[0])] = p
    for p in glob.glob(f"{SRC}/{h}/annual_peak_flow/*.csv"):
        peaks.add(sid(os.path.basename(p)[:-4]))

all_sites = sorted(set(gage_meta) | set(coeffs) | set(bath) | set(fm_src) | set(cf_src))

# ---------------------------------------------------------------------------
# 2. Lay out the package
# ---------------------------------------------------------------------------
for sub in ("field_measurements", "channel_features"):
    d = os.path.join(OUT, sub)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
os.makedirs(OUT, exist_ok=True)
print("Writing package files:")

# --- gages.csv / .geojson : the master gage inventory -----------------------
GAGE_F = ["site_id", "monitoring_location_id", "station_name", "lat", "lng", "state",
          "drain_area_sqkm", "stream_order", "gage_datum_alt_ft", "gagesii_class",
          "huc02", "huc4", "hcdn_2009", "active_2009",
          "has_coefficients", "has_bathymetry", "n_field_measurements",
          "n_channel_features", "has_annual_peaks"]
gage_rows = []
for s in all_sites:
    m = gage_meta.get(s, {})
    c = coeffs.get(s, {})
    b = bath.get(s, {}).get(30) or bath.get(s, {}).get(10) or {}
    lat = num(m.get("lat_gage") or c.get("lat_gage") or b.get("lat_gage"))
    lng = num(m.get("lng_gage") or c.get("lng_gage") or b.get("lng_gage"))
    nfm = ncf = 0
    if s in fm_src:
        with open(fm_src[s]) as fh:
            nfm = max(0, sum(1 for _ in fh) - 1)
    if s in cf_src:
        with open(cf_src[s]) as fh:
            ncf = max(0, sum(1 for _ in fh) - 1)
    gage_rows.append({
        "site_id": s,
        "monitoring_location_id": c.get("monitoring_location_id") or f"USGS-{s}",
        "station_name": m.get("staname", ""),
        "lat": lat, "lng": lng, "state": m.get("state", ""),
        "drain_area_sqkm": num(m.get("drain_sqkm") or c.get("drain_sqkm")),
        "stream_order": num(m.get("streamorde") or c.get("stream_order")),
        "gage_datum_alt_ft": num(m.get("alt_va") or c.get("raw_alt_va")),
        "gagesii_class": m.get("class", ""),
        "huc02": m.get("huc02", ""),
        "huc4": ";".join(sorted(gage_hucs[s])),
        "hcdn_2009": m.get("hcdn_2009", ""),
        "active_2009": m.get("active09", ""),
        "has_coefficients": int(s in coeffs),
        "has_bathymetry": int(s in bath),
        "n_field_measurements": nfm,
        "n_channel_features": ncf,
        "has_annual_peaks": int(s in peaks),
    })
write_csv(os.path.join(OUT, "gages.csv"), GAGE_F, gage_rows)
write_geojson(os.path.join(OUT, "gages.geojson"), gage_rows)

def fit_quality(site, r):
    """Recompute goodness-of-fit of the stored power laws against the
    channel-feature observations they were fit from (w=aQ^b, d=cQ^f, v=kQ^m)."""
    def fl(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    path = cf_src.get(site)
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    Q = np.array([fl(x.get("channel_flow")) for x in rows])
    W = np.array([fl(x.get("channel_width")) for x in rows])
    A = np.array([fl(x.get("channel_area")) for x in rows])
    V = np.array([fl(x.get("channel_velocity")) for x in rows])
    keep = np.isfinite(Q) & np.isfinite(W) & np.isfinite(A) & (Q > 0)
    Q, W, A, V = Q[keep], W[keep], A[keep], V[keep]
    Dp = np.divide(A, W, out=np.full_like(A, np.nan), where=W != 0)

    out = {"n_fit_points": int(len(Q))}
    for tag, obs, ck, ek in [("width", W, "coeff_a_width", "exp_b_width"),
                             ("depth", Dp, "coeff_c_depth", "exp_f_depth"),
                             ("velocity", V, "coeff_k_velocity", "exp_m_velocity")]:
        c, e = fl(r.get(ck)), fl(r.get(ek))
        g = np.isfinite(obs) & np.isfinite(Q)
        out[f"r2_{tag}"] = None
        if np.isfinite(c) and np.isfinite(e) and g.sum() >= 3:
            o, pred = obs[g], c * np.power(Q[g], e)
            st = float(np.sum((o - o.mean()) ** 2))
            if st > 0:
                out[f"r2_{tag}"] = round(1.0 - float(np.sum((o - pred) ** 2)) / st, 4)
    return out


# --- hydraulic_geometry_coefficients --------------------------------------
COEF_F = ["site_id", "monitoring_location_id", "huc4", "station_name", "lat", "lng",
          "stream_order", "drain_area_sqkm",
          "coeff_a_width", "exp_b_width", "coeff_c_depth", "exp_f_depth",
          "coeff_k_velocity", "exp_m_velocity", "exponent_sum",
          "max_flow_cfs", "n_channel_features", "n_fit_points",
          "r2_width", "r2_depth", "r2_velocity",
          "n_relations_fitted", "exponent_sum_in_range",
          "width_unit", "depth_unit", "velocity_unit", "flow_unit"]
byid = {r["site_id"]: r for r in gage_rows}
coef_rows = []
for s in sorted(coeffs):
    r, g = coeffs[s], byid[s]
    coef_rows.append({
        "site_id": s, "monitoring_location_id": r.get("monitoring_location_id", f"USGS-{s}"),
        "huc4": g["huc4"], "station_name": g["station_name"],
        "lat": num(r.get("lat_gage")), "lng": num(r.get("lng_gage")),
        "stream_order": num(r.get("stream_order")), "drain_area_sqkm": num(r.get("drain_sqkm")),
        "coeff_a_width": num(r.get("coeff_a_width")), "exp_b_width": num(r.get("exp_b_width")),
        "coeff_c_depth": num(r.get("coeff_c_depth")), "exp_f_depth": num(r.get("exp_f_depth")),
        "coeff_k_velocity": num(r.get("coeff_k_velocity")), "exp_m_velocity": num(r.get("exp_m_velocity")),
        "exponent_sum": num(r.get("exponent_sum")), "max_flow_cfs": num(r.get("max_flow_raw")),
        "n_channel_features": g["n_channel_features"],
        "width_unit": "ft", "depth_unit": "ft", "velocity_unit": "ft/s", "flow_unit": "ft^3/s",
    })
    es = num(r.get("exponent_sum"))
    coef_rows[-1]["exponent_sum_in_range"] = int(
        isinstance(es, (int, float)) and 0.85 <= es <= 1.15)
    # A failed fit leaves a null coefficient but an exponent of exactly 0.0, so it
    # still contributes to exponent_sum. Count the relations that actually converged
    # so that can be screened on directly.
    coef_rows[-1]["n_relations_fitted"] = sum(
        1 for ck in ("coeff_a_width", "coeff_c_depth", "coeff_k_velocity")
        if isinstance(num(r.get(ck)), (int, float)))
    coef_rows[-1].update(fit_quality(s, r))
write_csv(os.path.join(OUT, "hydraulic_geometry_coefficients.csv"), COEF_F, coef_rows)
write_geojson(os.path.join(OUT, "hydraulic_geometry_coefficients.geojson"), coef_rows)

# --- channel_bathymetry : wide over DEM resolution ------------------------
BATH_F = ["site_id", "huc4", "station_name", "lat", "lng", "stream_order", "drain_area_sqkm",
          "max_width_m", "max_depth_m", "gage_datum_alt_va_m", "riverbed_elev_m",
          "dem_elev_m_10m", "burn_value_m_10m", "dem_elev_m_30m", "burn_value_m_30m",
          "vertical_datum_type"]
bath_rows = []
for s in sorted(bath):
    v = bath[s]
    ref = v.get(30) or v.get(10)
    g = byid[s]
    bath_rows.append({
        "site_id": s, "huc4": g["huc4"], "station_name": g["station_name"],
        "lat": num(ref.get("lat_gage")), "lng": num(ref.get("lng_gage")),
        "stream_order": num(ref.get("stream_order")), "drain_area_sqkm": num(ref.get("drain_sqkm")),
        "max_width_m": num(ref.get("max_width_m")), "max_depth_m": num(ref.get("max_depth_m")),
        "gage_datum_alt_va_m": num(ref.get("alt_va_m")),
        "riverbed_elev_m": num(ref.get("riverbed_elev_m")),
        "dem_elev_m_10m": num(v.get(10, {}).get("dem_elev_m")),
        "burn_value_m_10m": num(v.get(10, {}).get("burn_value_m")),
        "dem_elev_m_30m": num(v.get(30, {}).get("dem_elev_m")),
        "burn_value_m_30m": num(v.get(30, {}).get("burn_value_m")),
        "vertical_datum_type": ref.get("datum_type", ""),
    })
write_csv(os.path.join(OUT, "channel_bathymetry.csv"), BATH_F, bath_rows)
write_geojson(os.path.join(OUT, "channel_bathymetry.geojson"), bath_rows)

# --- per-gage measurement files + concatenated long tables ----------------
def copy_set(src_map, subdir, out_all, extra_first="site_id"):
    allrows, fields = [], None
    for s in sorted(src_map):
        with open(src_map[s]) as fh:
            rdr = csv.DictReader(fh)
            rows = list(rdr)
            if fields is None:
                fields = [extra_first] + rdr.fieldnames
        dst = os.path.join(OUT, subdir, f"{s}.csv")
        with open(dst, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                r[extra_first] = s
                w.writerow(r)
                allrows.append(r)
    write_csv(os.path.join(OUT, out_all), fields, allrows)
    print(f"  {os.path.join('data', subdir):<52} {len(src_map):>5} per-gage files")

copy_set(fm_src, "field_measurements", "field_measurements_all.csv")
copy_set(cf_src, "channel_features", "channel_features_all.csv")

print(f"\nDomain totals: {len(all_sites)} gages | {len(coeffs)} with coefficients | "
      f"{len(bath)} with bathymetry | {len(fm_src)} field-measurement files | "
      f"{len(cf_src)} channel-feature files")
