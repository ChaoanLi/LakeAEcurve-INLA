# Lake Area-Elevation Curve Reconstruction via INLA-SPDE

**Authors**: Chaoan Li, Yinuo Zhu  
**Course**: STAT 647 Spatial Statistics  
**Institution**: Texas A&M University   
**Date**: December 2025

---

## 1. Introduction

### 1.1 Problem Statement

We aim to reconstruct lake bathymetry and estimate area-elevation (A-E) curves from remote sensing data. The A-E curve relates water surface elevation $h$ to inundated area:

$$A(h) = \int_{\mathcal{D}} \mathbb{1}\{Z(s) \leq h\} \, ds$$

where $Z(s)$ is the bottom elevation at location $s$.

### 1.2 Data Sources

| Data Type | Source | Resolution | Purpose |
|-----------|--------|------------|---------|
| Shoreline DEM | USGS | 30m | Elevation anchor points |
| Water Occurrence | Google Surface Water | 30m | Inundation frequency |
| Permanent Water Mask | Landsat | 30m | Deep water constraint |

---

## 2. Statistical Model

### 2.1 Latent Field Prior

The elevation field follows a Gaussian random field with Matérn covariance:

$$Z(s) \sim \text{GRF}(0, \mathcal{M}_\nu(\cdot; \rho, \sigma^2))$$

where $\rho$ is the practical correlation range and $\sigma^2$ is the marginal variance.

### 2.2 SPDE Representation

The Matérn field is approximated via the SPDE (Lindgren et al., 2011):

$$(\kappa^2 - \Delta)^{\alpha/2} Z(s) = \mathcal{W}(s)$$

with $\kappa = \sqrt{8\nu}/\rho$ and $\alpha = \nu + d/2$.

Discretization on a triangular mesh yields a sparse GMRF:

$$\mathbf{Z} \sim \mathcal{N}(\mathbf{0}, \mathbf{Q}^{-1}(\theta))$$

### 2.3 Observation Model

Water occurrence frequency is modeled as binomial:

$$Y_j \sim \text{Binomial}(N, p_j), \quad \text{logit}(p_j) = \alpha + f(s_j)$$

where $f(s)$ is the spatial random effect representing relative depth.

### 2.4 Priors

We use PC priors (Simpson et al., 2017):

- Range: $P(\rho < 500\text{m}) = 0.5$
- Sigma: $P(\sigma > 1\text{m}) = 0.01$
- Intercept: $\alpha \sim \mathcal{N}(0, 100)$

---

## 3. Calibration Framework

### 3.1 Identifiability Issue

The model estimates relative depth $f(s)$, not absolute elevation. Calibration maps $f(s) \mapsto Z(s)$ via:

$$Z(s) = a + b \cdot f(s)$$

### 3.2 Implemented Methods

**Method 1: Regression Calibration**

Fits $Z_j^{obs} = a + b \cdot f(s_j) + \epsilon_j$ using robust regression (rlm).

**Method 2: Endpoint-Constrained**

Forces:
- $\min(Z_{cal}) = z_{min}^{true}$
- $Z$ at max area $= z_{max}^{true}$

**Method 3: Hybrid**

Weighted average: $a = w \cdot a_{reg} + (1-w) \cdot a_{end}$ with $w = \min(R^2, 0.8)$.

**Method 4: Piecewise Affine**

Divides the elevation range into $K=4$ segments, fitting separate $(a_k, b_k)$ per segment to match the true A-E curve shape.

### 3.3 Model Selection

Best method selected by minimum MAE:

$$\text{best} = \arg\min_m \text{MAE}_m$$

### 3.4 Results (Belton Lake)

| Method | MAE (km²) | RMSE (km²) | R² |
|--------|-----------|------------|-----|
| Regression | 23.40 | — | — |
| Endpoint | 4.23 | — | — |
| Hybrid | 4.23 | — | — |
| **Piecewise** | **1.45** | — | **0.998** |

**Selected**: Piecewise (4 segments)

---

## 4. Implementation

### 4.1 Pipeline Overview

```
Data Loading → Mesh Construction → SPDE Definition → INLA Fitting
    ↓
Bathymetry Reconstruction → Calibration Framework → A-E Curve
```

### 4.2 Core Functions

| Module | Function | Purpose |
|--------|----------|---------|
| `data_generation.R` | `load_and_prep_data()` | Load and align rasters |
| `data_generation.R` | `build_observation_data()` | Construct binomial observations |
| `mesh_setup.R` | `build_mesh()` | SPDE triangulation |
| `spde_definition.R` | `define_spde()` | PC priors for spatial field |
| `fit_inlabru_model.R` | `fit_inla_model()` | INLA stack and fitting |
| `reconstruct_map.R` | `reconstruct_bathymetry()` | Posterior projection + calibration |
| `calibration_module.R` | `run_calibration_framework()` | 4-method comparison |
| `ae_curve_ppd.R` | `compute_ae_curve()` | A-E curve from bathymetry |

### 4.3 Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `max_edge` | (100m, 500m) | Balance resolution and computation |
| `prior_range` | (500m, 0.5) | Median range ~ lake scale |
| `prior_sigma` | (1m, 0.01) | Conservative variance bound |
| `N_trials` | 100 | Binomial discretization |
| `n_segments` | 4 | Piecewise calibration segments |

---

## 5. Results

### 5.1 Belton Lake

**Posterior Hyperparameters**:
- Spatial range: 1238 m
- Spatial sigma: 29.9 m

**A-E Curve Performance**:
- Elevation range: 147.09 – 181.19 m
- Area range: 0 – 50.09 km²
- MAE: 1.45 km²
- R²: 0.998

### 5.2 Validation

The predicted A-E curve closely matches TWDB survey data across the full elevation range. Piecewise calibration effectively captures the nonlinear relationship between model-predicted depth and true elevation.

---

## 6. Presentation Outline

### Slide 1: Title

**Lake Bathymetry Reconstruction via INLA-SPDE**

Chaoan Li, Yinuo Zhu  
Princeton University | STAT 647 | December 2025

---

### Slide 2: Motivation

**Problem**: Estimate lake depth from satellite observations

**Applications**:
- Reservoir capacity estimation
- Flood risk assessment
- Drought monitoring

**Challenge**: No direct underwater measurements available

---

### Slide 3: Data

Three raster inputs (30m resolution):

| Type | Description |
|------|-------------|
| Shoreline DEM | Elevation at water edge |
| Water Frequency | Fraction of time inundated |
| Permanent Water | Always-wet regions |

Study sites: Belton Lake, E.V. Spence Reservoir (Texas)

---

### Slide 4: Model

**Latent field**: $Z(s) \sim \text{GRF with Matérn covariance}$

**SPDE approximation**: 
$$(\kappa^2 - \Delta)^{\alpha/2} Z = \mathcal{W}$$

**Observation model**:
$$Y \sim \text{Binomial}(N, p), \quad \text{logit}(p) = \alpha + f(s)$$

**Priors**: PC priors on range and sigma

---

### Slide 5: Mesh and SPDE

- Constrained Delaunay triangulation
- 1523 vertices, adaptive resolution
- Inner edge: 100m, outer: 500m

[Figure: Mesh overlaid on lake boundary]

---

### Slide 6: Calibration

**Issue**: Model gives relative depth, need absolute elevation

**Solution**: Piecewise affine calibration

$$Z^{(k)} = a_k + b_k \cdot f, \quad k = 1, \ldots, 4$$

**Comparison**:

| Method | MAE (km²) |
|--------|-----------|
| Regression | 23.40 |
| Endpoint | 4.23 |
| Hybrid | 4.23 |
| **Piecewise** | **1.45** |

---

### Slide 7: Bathymetry Results

**Posterior mean elevation**:
- Range: 147 – 181 m
- Deeper regions (blue) near dam
- Shallow areas (yellow) near shore

**Posterior uncertainty**:
- Higher in deep water (far from observations)
- Lower near shoreline constraints

[Figure: Bathymetry mean and SD maps]

---

### Slide 8: A-E Curve Validation

**Predicted vs. True**:
- Black: TWDB survey data
- Blue: INLA-SPDE prediction

**Metrics**:
- MAE = 1.45 km²
- R² = 0.998
- Elevation: 147 – 181 m
- Area: 0 – 50 km²

[Figure: A-E curve comparison]

---

### Slide 9: Model Diagnostics

**Hyperparameters**:
- Range: 1238 m (lake-scale correlation)
- Sigma: 29.9 m (depth variation)

**Calibration quality**:
- 4-segment piecewise fit
- Captures nonlinear shape

**Uncertainty**:
- Increases with distance from shore
- 95% credible intervals available

---

### Slide 10: Conclusions

**Findings**:
1. INLA-SPDE effectively reconstructs bathymetry from water frequency
2. Piecewise calibration reduces MAE by 97% vs. regression
3. Final MAE = 1.45 km² (R² = 0.998)

**Contributions**:
- Integration of SPDE spatial modeling with remote sensing
- Automated calibration framework with model selection
- Reproducible workflow for any lake

**Limitations**:
- Requires validation data for calibration
- Assumes stationary spatial field

---

### Slide 11: References

1. Rue, H., Martino, S., & Chopin, N. (2009). INLA. *JRSSB*.
2. Lindgren, F., Rue, H., & Lindström, J. (2011). SPDE approach. *JRSSB*.
3. Simpson, D. et al. (2017). PC priors. *Statistical Science*.
4. Pekel, J.F. et al. (2016). Global surface water. *Nature*.

---

## 7. References

1. **Rue, H., Martino, S., & Chopin, N.** (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B*, 71(2), 319-392.

2. **Lindgren, F., Rue, H., & Lindström, J.** (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *Journal of the Royal Statistical Society: Series B*, 73(4), 423-498.

3. **Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.** (2017). Penalising model component complexity: A principled, practical approach to constructing priors. *Statistical Science*, 32(1), 1-28.

4. **Pekel, J. F., Cottam, A., Gorelick, N., & Belward, A. S.** (2016). High-resolution mapping of global surface water and its long-term changes. *Nature*, 540(7633), 418-422.

5. **Gao, H.** (2015). Satellite remote sensing of large lakes and reservoirs: from elevation and area to storage. *Wiley Interdisciplinary Reviews: Water*, 2(2), 147-157.

---

## Appendix A: Troubleshooting

**INLA convergence failure**:
- Check for extreme values in data
- Use `control.inla = list(int.strategy = "eb")`

**A-E curve offset**:
- Verify dam point coordinates
- Increase dam point weight in `build_observation_data()`

**High uncertainty**:
- Add more shoreline observations
- Refine mesh resolution in deep water

---

## Appendix B: Parameter Tuning

| Parameter | Range | Effect |
|-----------|-------|--------|
| `prior_range[1]` | 200-1000m | Spatial smoothness |
| `max_edge[1]` | 50-200m | Resolution |
| `n_segments` | 3-6 | Calibration flexibility |

---

*Document version: Final (December 2025)*
