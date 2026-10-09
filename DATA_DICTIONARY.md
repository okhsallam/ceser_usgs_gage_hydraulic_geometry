# Data Dictionary

`site_id` (8-digit zero-padded USGS site number, string) is the join key in every file.

---

## `data/gages.csv` / `gages.geojson` — 88 rows

Master inventory. One row per gage in the domain, with flags for what else exists.

| Column | Type | Units | Description |
|---|---|---|---|
| `site_id` | string | — | USGS site number, 8-digit zero-padded |
| `monitoring_location_id` | string | — | USGS OGC-API id, `USGS-<site_id>` |
| `station_name` | string | — | USGS station name (GAGES-II `staname`) |
| `lat`, `lng` | float | deg | Gage location, WGS84 |
| `state` | string | — | US state postal code |
| `drain_area_sqkm` | float | km² | Upstream drainage area |
| `stream_order` | int | — | NHDPlus Strahler stream order |
| `gage_datum_alt_ft` | float | ft | Gage datum elevation (USGS `alt_va`) |
| `gagesii_class` | string | — | GAGES-II class: `Ref` (reference) or `Non-ref` |
| `huc02` | string | — | HUC2 region code |
| `huc4` | string | — | Source HUC4 folder(s), `;`-joined. Provenance only |
| `hcdn_2009` | string | — | HCDN-2009 membership, blank if not a member |
| `active_2009` | string | — | Active as of 2009 (`yes`/blank) |
| `has_coefficients` | 0/1 | — | Row present in `hydraulic_geometry_coefficients.csv` |
| `has_bathymetry` | 0/1 | — | Row present in `channel_bathymetry.csv` |
| `n_field_measurements` | int | — | Stage–discharge measurement count (0 if no file) |
| `n_channel_features` | int | — | Channel measurement count (0 if no file) |
| `has_annual_peaks` | 0/1 | — | Annual peak-flow series exists upstream (not shipped here) |

---

## `data/hydraulic_geometry_coefficients.csv` / `.geojson` — 84 rows

Power-law fits `w = a·Q^b`, `d = c·Q^f`, `v = k·Q^m`. **Imperial units.**

| Column | Type | Units | Description |
|---|---|---|---|
| `site_id` | string | — | Join key |
| `monitoring_location_id` | string | — | `USGS-<site_id>` |
| `huc4` | string | — | Source HUC4, provenance only |
| `station_name` | string | — | USGS station name |
| `lat`, `lng` | float | deg | Gage location, WGS84 |
| `stream_order` | int | — | NHDPlus Strahler stream order |
| `drain_area_sqkm` | float | km² | Upstream drainage area |
| `coeff_a_width` | float | ft / (ft³/s)^b | Width coefficient `a` |
| `exp_b_width` | float | — | Width exponent `b`, bounded [0, 1]. **0.0 exactly means the fit failed** |
| `coeff_c_depth` | float | ft / (ft³/s)^f | Depth coefficient `c`. **Null for 3 gages** |
| `exp_f_depth` | float | — | Depth exponent `f`, bounded [0, 1]. **0.0 exactly means the fit failed** |
| `coeff_k_velocity` | float | (ft/s) / (ft³/s)^m | Velocity coefficient `k`. **Null for 4 gages** |
| `exp_m_velocity` | float | — | Velocity exponent `m`, bounded [0, 1]. **0.0 exactly means the fit failed** |
| `exponent_sum` | float | — | `b + f + m`. Theory expects ≈ 1; observed 0.155–1.645. Includes the 0.0 of any failed fit |
| `max_flow_cfs` | float | ft³/s | Largest discharge in the fitting sample |
| `n_channel_features` | int | — | Rows in the gage's channel-features file |
| `n_fit_points` | int | — | Rows actually usable in the fit (finite, `Q > 0`) |
| `r2_width` | float | — | R² of `w = a·Q^b` vs observations. Median 0.62 |
| `r2_depth` | float | — | R² of `d = c·Q^f` vs observations. Median 0.83 |
| `r2_velocity` | float | — | R² of `v = k·Q^m` vs observations. Median 0.71 |
| `n_relations_fitted` | int | — | How many of the three fits converged (3 = all; 78 gages, 5 have 2, 1 has 1). **Screen on this** — see the warning below |
| `exponent_sum_in_range` | 0/1 | — | 1 if `exponent_sum` ∈ [0.85, 1.15]. True for 72 of 84. Not sufficient on its own |
| `width_unit`, `depth_unit`, `velocity_unit`, `flow_unit` | string | — | Unit declarations: `ft`, `ft`, `ft/s`, `ft^3/s` |

R² is recomputed by `build_package.py` against the same observations the fit used, in
linear space; it is a QC aid added here, not an upstream output. Null where the
coefficient is null or fewer than 3 usable points exist.

> **Failed fits are not null throughout.** When a relation fails to converge its
> coefficient is null but its exponent is stored as exactly 0.0, and `exponent_sum` still
> adds that zero in. Gage 05466500 has a failed velocity fit yet `exponent_sum` = 1.089,
> so `exponent_sum_in_range` = 1. Always pair it with `n_relations_fitted == 3`.

Full derivation of every column here: **[METHODS.md](METHODS.md)**.

---

## `data/channel_bathymetry.csv` / `.geojson` — 73 rows

Riverbed elevation and DEM burn depth, at two DEM resolutions. **Metric units.**

| Column | Type | Units | Description |
|---|---|---|---|
| `site_id` | string | — | Join key |
| `huc4` | string | — | Source HUC4, provenance only |
| `station_name` | string | — | USGS station name |
| `lat`, `lng` | float | deg | Gage location, WGS84 |
| `stream_order` | int | — | NHDPlus Strahler stream order |
| `drain_area_sqkm` | float | km² | Upstream drainage area |
| `max_width_m` | float | m | `a·Q_max^b` at the largest measured discharge — an extrapolation from the fit, not a measured maximum |
| `max_depth_m` | float | m | `c·Q_max^f` at the largest measured discharge — likewise extrapolated |
| `gage_datum_alt_va_m` | float | m | Gage datum elevation, `alt_va` × 0.3048 |
| `riverbed_elev_m` | float | m | `min` over measurements of (gage datum + stage − fitted depth). Inherits the depth fit's error |
| `dem_elev_m_10m` | float | m | Min DEM elevation within 150 m of the gage, 10 m DEM. Null for 1 gage |
| `burn_value_m_10m` | float | m | `dem_elev_m_10m − riverbed_elev_m`. Always > 0. Null for 1 gage |
| `dem_elev_m_30m` | float | m | Min DEM elevation within 150 m of the gage, 30 m DEM |
| `burn_value_m_30m` | float | m | `dem_elev_m_30m − riverbed_elev_m`. Always > 0 |
| `vertical_datum_type` | string | — | Reported datum. **`Unknown` for all rows** — see limitation 4 |

`dem_elev_m_*` is the **minimum** DEM elevation inside a 150 m-radius circle around the
gage, which finds the channel thalweg rather than the bank the gage house sits on.

Only `dem_elev_m` and `burn_value_m` differ between the 10 m and 30 m runs (69 of 73
gages), so all other columns are resolution-independent and stored once.

Full derivation: **[METHODS.md](METHODS.md)**.

---

## `data/field_measurements/<site_id>.csv` — 82 files, 9,183 rows total

USGS stage–discharge field measurements. Input to the bathymetry/burn step.
`field_measurements_all.csv` is these files concatenated.

| Column | Type | Units | Description |
|---|---|---|---|
| `site_id` | string | — | Join key (added by this package) |
| `MeanGageHeight` | float | ft | Mean gage height (stage) during the measurement |
| `Discharge` | float | ft³/s | Measured discharge |
| `alt_va` | float | ft | Gage datum elevation, repeated per row |
| `v_datum_type` | string | — | Vertical datum: `NAVD88` (8,725 rows) or `Unknown` (458) |
| `streamorde` | int | — | NHDPlus Strahler stream order, repeated per row |
| `drain_sqkm` | float | km² | Drainage area, repeated per row |
| `lat_gage`, `lng_gage` | float | deg | Gage location, repeated per row |

---

## `data/channel_features/<site_id>.csv` — 86 files, 44,251 rows total

USGS channel measurements. Input to the hydraulic-geometry fits.
`channel_features_all.csv` is these files concatenated. Columns are passed through from
the USGS OGC API unchanged apart from the added `site_id`.

| Column | Type | Units | Description |
|---|---|---|---|
| `site_id` | string | — | Join key (added by this package) |
| `geometry` | WKT | deg | `POINT (lng lat)`, WGS84 |
| `id`, `field_visit_id`, `measurement_number` | string | — | USGS record identifiers |
| `monitoring_location_id` | string | — | `USGS-<site_id>` |
| `time` | datetime | — | Measurement timestamp |
| `channel_name` | string | — | Channel/section label |
| `channel_flow` | float | ft³/s | Discharge — the `Q` in every fit |
| `channel_width` | float | ft | Measured top width — the `w` in the width fit |
| `channel_area` | float | ft² | Measured cross-section area; `d = area / width` |
| `channel_velocity` | float | ft/s | Mean velocity — the `v` in the velocity fit |
| `channel_*_unit` | string | — | Unit declaration for the matching value column |
| `channel_location_distance` | float | ft | Distance of the section from the gage |
| `channel_location_direction` | string | — | Direction of the section from the gage |
| `channel_stability` | string | — | Qualitative bed stability rating |
| `channel_material` | string | — | Qualitative bed material description |
| `channel_evenness` | string | — | Qualitative cross-section evenness rating |
| `horizontal_velocity_description` | string | — | Qualitative horizontal velocity distribution |
| `vertical_velocity_description` | string | — | Qualitative vertical velocity distribution |
| `longitudinal_velocity_description` | string | — | Qualitative longitudinal velocity distribution |
| `measurement_type` | string | — | USGS measurement type |
| `channel_measurement_type` | string | — | USGS channel measurement subtype |
| `last_modified` | datetime | — | USGS record last-modified timestamp |

Rows are used in a fit only where `channel_flow`, `channel_width` and `channel_area` are
all finite and `channel_flow > 0`; `n_fit_points` in the coefficients file reports the
surviving count per gage.
