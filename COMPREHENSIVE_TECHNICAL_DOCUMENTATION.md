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
| Dam Point | Manual digitization | Point | Minimum elevation anchor |

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

The model estimates relative depth $f(s)$, not absolute elevation. Calibration maps $f(s) \mapsto Z(s)$ via an affine transformation:

$$Z(s) = a + b \cdot (-f(s))$$

where the negative sign accounts for the inverse relationship between field values and elevation.

### 3.2 Design Principle

**Key insight**: DEM data covers the shoreline region (above water). All calibration methods use the same DEM-based calibration for this region. Methods differ only in how they extrapolate into the underwater region.

### 3.3 Implemented Methods

**Method 1: DEM-only (0 validation points)**

Uses quantile matching on shoreline DEM to derive $(a, b)$, then extrapolates to lake bottom using the same slope.

**Method 2: DEM+1pt (1 validation point)**

- Above water: Same DEM calibration
- Underwater: Interpolate from DEM boundary to dam point (minimum elevation)

**Method 3: DEM+2pt (2 validation points)**

- Above water: Same DEM calibration
- Underwater: Two segments using 0% and 50% area quantiles

**Method 4: DEM+4pt (4 validation points)**

- Above water: Same DEM calibration
- Underwater: Four segments using 0%, 25%, 50%, 75% area quantiles

### 3.4 Model Selection

Best method selected by minimum MAE:

$$\text{best} = \arg\min_m \text{MAE}_m$$

### 3.5 Results (Belton Lake)

| Method | Val. Pts | MAE (km²) | R² |
|--------|----------|-----------|-----|
| DEM-only | 0 | 23.40 | — |
| DEM+1pt | 1 | 4.23 | — |
| DEM+2pt | 2 | 2.89 | — |
| **DEM+4pt** | **4** | **1.45** | **0.998** |

**Selected**: DEM+4pt (piecewise with 4 segments)

---

## 4. Implementation

### 4.1 Pipeline Overview

```
Step 1:  load_and_prep_data()       → Load DEM, water frequency, permanent water
Step 2:  build_observation_data()   → Construct binomial observations + dam constraint
Step 3:  plot_original_data()       → Visualize input data
Step 4:  build_mesh()               → SPDE triangulation
Step 5:  plot_mesh()                → Visualize mesh
Step 6:  define_spde()              → PC priors for spatial field
Step 7:  build_projection_matrices()→ A matrices for observations
Step 8:  build_inla_stacks()        → Combine data for INLA
Step 9:  fit_inla_model()           → Binomial INLA fitting
Step 10: extract_spde_hyperpar()    → Get range/sigma posteriors
Step 11: reconstruct_bathymetry()   → Project + calibrate (4 methods)
Step 12: plot_bathymetry()          → Visualize mean/SD maps
Step 13: compute_ae_curve()         → Calculate A-E relationship
Step 14: plot_ae_curve()            → Compare with validation
Step 15: plot_calibration_comparison() → 4-method comparison
```

### 4.2 Core Functions

| Module | Function | Purpose |
|--------|----------|---------|
| `data_generation.R` | `load_and_prep_data()` | Load and align rasters, CRS projection |
| `data_generation.R` | `build_observation_data()` | Construct binomial obs + dam constraint |
| `mesh_setup.R` | `build_mesh()` | Constrained Delaunay triangulation |
| `mesh_setup.R` | `plot_mesh()` | Mesh visualization |
| `spde_definition.R` | `define_spde()` | PC priors on range/sigma |
| `spde_definition.R` | `build_projection_matrices()` | A matrices |
| `fit_inlabru_model.R` | `build_inla_stacks()` | Stack for binomial likelihood |
| `fit_inlabru_model.R` | `fit_inla_model()` | INLA fitting |
| `fit_inlabru_model.R` | `extract_spde_hyperpar()` | Transform θ to (ρ, σ) |
| `reconstruct_map.R` | `reconstruct_bathymetry()` | Posterior projection + calibration |
| `reconstruct_map.R` | `plot_bathymetry()` | Mean/SD raster visualization |
| `calibration_module.R` | `run_calibration_framework()` | 4-method comparison + auto-select |
| `calibration_module.R` | `compute_ae_metrics()` | MAE, RMSE, R², MAPE, NSE |
| `calibration_module.R` | `plot_calibration_comparison()` | 4-method A-E curve comparison |
| `ae_curve_ppd.R` | `compute_ae_curve()` | A-E curve from calibrated bathymetry |
| `ae_curve_ppd.R` | `plot_ae_curve()` | Compare predicted vs. true A-E |

### 4.3 Key Parameters

| Parameter | Value | Location | Rationale |
|-----------|-------|----------|-----------|
| `max_edge` | (100m, 500m) | `build_mesh()` | Inner/outer mesh resolution |
| `cutoff` | 50m | `build_mesh()` | Minimum node spacing |
| `offset` | (100m, 500m) | `build_mesh()` | Inner/outer boundary buffer |
| `prior_range` | (500m, 0.5) | `define_spde()` | P(ρ < 500) = 0.5 |
| `prior_sigma` | (1m, 0.01) | `define_spde()` | P(σ > 1) = 0.01 |
| `N_trials` | 100 | `build_observation_data()` | Binomial discretization |
| `dam_point_weight` | 50 | `build_observation_data()` | Dam point replication |
| `elevation_step` | 0.1m | `compute_ae_curve()` | A-E curve resolution |

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

### 5.2 Output Files

| File | Description |
|------|-------------|
| `Belton_original_data.png` | Input data visualization |
| `Belton_mesh.png` | SPDE mesh |
| `Belton_bathymetry_mean.png` | Posterior mean elevation |
| `Belton_bathymetry_sd.png` | Posterior standard deviation |
| `Belton_ae_curve.png` | Predicted vs. true A-E curve |
| `Belton_calibration_comparison.png` | 4-method calibration comparison |
| `Belton_ae_curve.csv` | Numerical A-E data |
| `Belton_results.RData` | Full R objects |

---

## 6. References

1. **Rue, H., Martino, S., & Chopin, N.** (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B*, 71(2), 319-392.

2. **Lindgren, F., Rue, H., & Lindström, J.** (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *Journal of the Royal Statistical Society: Series B*, 73(4), 423-498.

3. **Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.** (2017). Penalising model component complexity: A principled, practical approach to constructing priors. *Statistical Science*, 32(1), 1-28.

4. **Pekel, J. F., Cottam, A., Gorelick, N., & Belward, A. S.** (2016). High-resolution mapping of global surface water and its long-term changes. *Nature*, 540(7633), 418-422.

---

## Appendix A: Troubleshooting

**INLA convergence failure**:
- Check for extreme values in data
- Use `control.inla = list(int.strategy = "eb")`

**A-E curve offset**:
- Verify dam point coordinates
- Increase `dam_point_weight` in `build_observation_data()`

**High uncertainty in deep water**:
- Increase shoreline observations
- Decrease `max_edge[1]` for finer mesh

---

## Appendix B: Parameter Tuning Guide

| Parameter | Range | Effect |
|-----------|-------|--------|
| `prior_range[1]` | 200-1000m | Larger = smoother field |
| `max_edge[1]` | 50-200m | Smaller = higher resolution |
| `dam_point_weight` | 20-100 | Larger = stronger anchoring |

---

*Document version: Final (December 2025)*
