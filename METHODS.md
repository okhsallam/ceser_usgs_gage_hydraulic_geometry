# Methods — the math behind the package

Everything here is transcribed from the code that produced the data
(`watersheds_ceser_hydro_prepare_new.ipynb`, cells 6 and 9) and from
`build_package.py` in this repo. Where the implementation departs from the
textbook form, that is stated rather than smoothed over.

**Notation.** $Q$ discharge, $w$ top width, $A$ cross-section area, $d$ mean depth,
$v$ mean velocity, $H$ mean gage height (stage), $z$ elevation. Subscript $i$ indexes
one field measurement. Imperial units unless a symbol carries $_m$ (metres).

---

## 1. At-a-station hydraulic geometry

### 1.1 The relations

Leopold & Maddock (1953) describe how a channel adjusts at a single cross-section as
discharge rises, through three power laws:

$$w = a\,Q^{b} \qquad d = c\,Q^{f} \qquad v = k\,Q^{m}$$

Mean depth is not measured directly; it is derived from the gaged area and width:

$$d_i = \frac{A_i}{w_i}$$

### 1.2 The continuity constraint

Discharge is width × depth × velocity, so for every $Q$:

$$Q = w\,d\,v = (a\,Q^{b})(c\,Q^{f})(k\,Q^{m}) = a\,c\,k\;Q^{\,b+f+m}$$

which holds for all $Q$ only if

$$\boxed{\;b+f+m = 1\quad\text{and}\quad a\,c\,k = 1\;}$$

The three relations are fit **independently**, so neither identity is imposed — they
are diagnostics, not constraints. `exponent_sum` stores $b+f+m$ and
`exponent_sum_in_range` flags $0.85 \le b+f+m \le 1.15$. A worked check on gage
05416900 at $Q = 5000$ ft³/s gives $w\,d\,v = 5070$ ft³/s, a closure error of 1.4%.

### 1.3 How the coefficients were fit

Bounded non-linear least squares (`scipy.optimize.curve_fit`, Trust Region Reflective)
on **untransformed** values. For each relation independently, with $y_i$ the observed
width, depth or velocity:

$$\min_{\theta=(\kappa,\,\epsilon)} \sum_i \bigl[\, y_i - \kappa\,Q_i^{\epsilon} \,\bigr]^2
\qquad \text{s.t.}\quad \kappa \in [10^{-3},\,\infty),\;\; \epsilon \in [0,\,1]$$

starting from $\theta_0 = (1.0,\,0.3)$, `maxfev=3000`.

Observations enter the fit only where $Q_i$, $w_i$ and $A_i$ are all finite and
$Q_i > 0$; a gage needs at least 3 such rows. The surviving count per gage is
`n_fit_points`.

> **This is least squares in linear space, not the usual log-space regression.**
> Fitting $\log y = \log\kappa + \epsilon \log Q$ by OLS — the conventional approach —
> weights all decades of discharge equally. Minimising squared error on untransformed
> values instead lets the largest flows dominate the objective, because their residuals
> are numerically largest. Expect these coefficients to track high flows more closely
> and low flows more loosely than a log-space fit of the same data would. For
> flood-modelling use that bias is benign, arguably helpful; for low-flow work it is not.

### 1.4 What a failed fit looks like

Each relation is fit in its own `try`. On failure the coefficient is recorded as
**null** but the exponent keeps its initialised value of **exactly 0.0** — not null.
Since `exponent_sum` is computed as $b+f+m$ regardless, a failed relation silently
contributes a zero.

This produces a trap worth knowing about. Gage **05466500** (Edwards River) has a
failed velocity fit — `coeff_k_velocity` is null and `exp_m_velocity` is 0.0 — yet
$b+f+m = 0.772 + 0.317 + 0 = 1.089$ falls inside the valid band, so
`exponent_sum_in_range` reads 1. The flag alone would wave it through.

Screen on **`n_relations_fitted`** (how many of the three converged: 78 gages have all
3, five have 2, one has 1) rather than on `exponent_sum_in_range` alone.

### 1.5 Goodness of fit

`r2_width`, `r2_depth` and `r2_velocity` are computed by `build_package.py` in this
repo — they are not upstream outputs. Each is the ordinary coefficient of
determination, in **linear space**, against the same observations the fit used:

$$R^2 = 1 - \frac{\sum_i (y_i - \hat{y}_i)^2}{\sum_i (y_i - \bar{y})^2},
\qquad \hat{y}_i = \kappa\,Q_i^{\epsilon}$$

Null where the coefficient is null, where fewer than 3 usable points exist, or where
$\sum_i (y_i-\bar y)^2 = 0$. Because the fit minimises exactly this numerator in this
space, $R^2$ here is a fair in-sample measure of that objective — it is not a
cross-validated or out-of-sample score, and it is not comparable to an $R^2$ from a
log-space fit.

---

## 2. Channel bathymetry and DEM burn depth

The goal is $z_{bed}$, the absolute riverbed elevation, and from it the **burn depth** —
how far a DEM must be lowered so its channel sits at the true bed.

### 2.1 Bed elevation from stage and fitted depth

This step uses the **stage–discharge** field measurements ($H_i$, $Q_i$) together with
the depth power law fit in §1.

For each measurement, the water-surface elevation is the gage datum plus the stage, and
the bed beneath it is that surface minus the water depth at that discharge. The depth is
not measured — it is predicted from the fitted relation:

$$\hat{d}_i = c\,Q_i^{\,f} \qquad\text{(ft, then} \times 0.3048 \text{ to m)}$$

$$z_{bed,i} = \bigl(z_{datum} + H_i\bigr) - \hat{d}_i \qquad \text{(all in m)}$$

where $z_{datum} = \texttt{alt\_va} \times 0.3048$ is the gage datum elevation. A single
bed elevation is then taken as the **minimum** over all measurements:

$$\boxed{\;z_{bed} = \min_i\Bigl[\bigl(z_{datum} + H_i\bigr) - c\,Q_i^{\,f}\cdot 0.3048\Bigr]\;}$$

Taking the minimum, rather than a mean, makes the estimate the deepest bed consistent
with any observation — a conservative choice for burning, and one sensitive to a single
outlying measurement.

> This is **not** a zero-flow stage extrapolation. The riverbed is reconstructed
> measurement by measurement from the fitted depth law. It therefore inherits the
> depth fit's error directly: where `r2_depth` is poor, `riverbed_elev_m` is poor.

### 2.2 DEM elevation at the gage

$z_{DEM}$ is the **minimum** DEM elevation within a 150 m radius of the gage. The gage
lat/lon is projected to EPSG:5070, a square window of ±150 m is sliced from the raster,
a circular mask $\sqrt{(X-x)^2+(Y-y)^2} \le 150$ m is applied, and the minimum of the
masked cells is taken. Using the minimum within a neighbourhood, rather than the value
at the gage pixel, is what finds the channel thalweg instead of the bank the gage
house sits on.

Computed at both 10 m and 30 m DEM resolution (`dem_elev_m_10m`, `dem_elev_m_30m`).

### 2.3 Burn depth

$$\boxed{\;\Delta z_{burn} = z_{DEM} - z_{bed}\;}$$

Positive means the DEM sits above the reconstructed bed and must be lowered by
$\Delta z_{burn}$ to seat the channel.

**A gage is retained only if $\Delta z_{burn} > 0$.** Gages where the DEM already sat at
or below the reconstructed bed were dropped upstream, which is most of the 84 → 73
attrition. So `burn_value_m` is positive everywhere in this file *by construction*, and
the 73 are a biased subset — not evidence that every channel needs burning.

### 2.4 Channel dimensions at the peak measured flow

Evaluated from the fitted laws at $Q_{max}$, the largest discharge in the gage's
stage–discharge record:

$$w_{max} = a\,Q_{max}^{\,b}\cdot 0.3048 \;\text{m}
\qquad d_{max} = c\,Q_{max}^{\,f}\cdot 0.3048 \;\text{m}$$

These are extrapolations from the fit, not measured maxima.

### 2.5 The riverbed stage plotted in the figures

`scripts/plot_gage_fits.py` draws the bed on the stage–discharge panel by inverting
§2.1 back to gage-height units:

$$H_{bed} = \frac{z_{bed} - z_{datum}}{0.3048} \quad\text{[ft]}$$

the gage height at which the fitted depth law gives zero water depth. It can be
negative where the bed lies below the gage datum (e.g. gage 05437050, $-0.77$ ft).

---

## 3. Units

| Quantity | Unit | Where |
|---|---|---|
| $Q$, `channel_flow`, `Discharge`, `max_flow_cfs` | ft³/s | coefficients, channel features, field measurements |
| $w$, $d$, $H$, `alt_va` | ft | coefficients, channel features, field measurements |
| $v$ | ft/s | coefficients, channel features |
| $a, c, k$ | ft/(ft³/s)$^{b,f}$, (ft/s)/(ft³/s)$^m$ | coefficients |
| $b, f, m$ | dimensionless | coefficients |
| all `*_m` columns | m | channel bathymetry |
| drainage area | km² | everywhere |

The single conversion used throughout is $1\ \text{ft} = 0.3048\ \text{m}$ exactly.
**The coefficients file is imperial and the bathymetry file is metric** — they are not
in the same unit system, which is the most likely way to misuse this package.

---

## 4. Known departures from textbook practice

1. Least squares on untransformed values, not log-space OLS (§1.3) — biases the fit
   toward high flows.
2. Continuity ($b+f+m=1$, $ack=1$) is diagnosed, never imposed (§1.2).
3. Failed fits store exponent 0.0, not null, and still enter `exponent_sum` (§1.4).
4. Bed elevation is a minimum over measurements, so one outlier can set it (§2.1).
5. Burn depths are positive-only by construction, making the 73-gage set biased (§2.3).
6. Vertical datum is `Unknown` for all 84 gages in the coefficient source, though the
   field measurements report NAVD88 for 8,725 of 9,183 rows. Bed and burn elevations are
   NAVD88-**assumed**, not datum-verified.

---

## References

Leopold, L.B. and Maddock, T. (1953). *The Hydraulic Geometry of Stream Channels and
Some Physiographic Implications.* USGS Professional Paper 252.
