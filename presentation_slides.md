---
marp: true
theme: default
paginate: true
math: katex
---

# Lake Bathymetry Reconstruction via INLA-SPDE
## Bayesian Spatial Modeling of Area-Elevation Curves

**Chaoan Li, Yinuo Zhu**

STAT 647 Spatial Statistics | Texas A&M University

December 2025

<!-- notes -->
Good afternoon everyone. Today I'll present our project on reconstructing lake bathymetry using INLA-SPDE methods. This is joint work with Yinuo Zhu. We'll show how to estimate underwater topography from satellite imagery alone.

---

# Motivation: Why Area-Elevation Curves?

**Problem**: Lake depth is hard to measure directly
- Sonar surveys are expensive ($50k+ per lake)
- Many lakes have no bathymetric data
- Water managers need A-E curves for:
  - Reservoir capacity estimation
  - Flood risk assessment
  - Drought monitoring

**Our Goal**: Reconstruct A-E curve $A(h)$ from remote sensing

$$A(h) = \int_{\mathcal{D}} \mathbb{1}\{Z(s) \leq h\} \, ds$$

<!-- notes -->
Why do we care about this? Lake managers need to know how much water is stored at different levels. Traditional surveys cost tens of thousands of dollars. We're asking: can we estimate this from free satellite data?

---

# Study Site: Belton Lake, Texas

{{Lake Location Map - optional}}

- Reservoir on Leon River, Central Texas
- Surface area: ~50 km² at full pool
- Managed by U.S. Army Corps of Engineers
- **Validation data available** from Texas Water Development Board

**Why Belton?**
- Clear water occurrence patterns
- Survey data for validation
- Representative of Texas reservoirs

<!-- notes -->
We chose Belton Lake because we have ground truth data to validate our model. It's a typical Texas reservoir - not too large, not too small, and has well-documented elevation-area relationships from actual surveys.

---

# Data Sources (All Free & Public)

| Data | Source | Resolution | What it tells us |
|------|--------|------------|-----------------|
| **DEM** | USGS | 30m | Shoreline elevations |
| **Water Occurrence** | Google Surface Water | 30m | Inundation frequency $p \in [0,1]$ |
| **Permanent Water** | Landsat time series | 30m | Always-wet pixels |

{{Original Data Visualization - 3 panels}}

<!-- notes -->
All our input data is freely available. The DEM gives us elevation at the shoreline. Water occurrence from Google tells us what fraction of time each pixel was underwater over 30 years. Permanent water marks the deepest areas that never dry out.

---

# Key Statistical Challenges

1. **Non-Gaussian observations**
   - Water frequency $p \in [0,1]$ is not normal
   - Need binomial likelihood with logit link

2. **Spatial misalignment**
   - DEM covers land, water frequency covers lake
   - Different observation types at different locations

3. **Identifiability**
   - Model gives *relative* depth, not *absolute* elevation
   - Need calibration to anchor predictions

4. **Computational scale**
   - Thousands of spatial locations
   - Dense covariance matrix is $O(n^3)$

<!-- notes -->
These are the four main challenges. The data isn't Gaussian. Different data types are observed at different places. The model output needs calibration. And naively, spatial models are computationally expensive. Each of these drove our modeling decisions.

---

# Hierarchical Model Framework

**Layer 1: Latent Elevation Field**
$$Z(s) \sim \text{GRF}\big(0, \, \mathcal{M}_\nu(\cdot; \rho, \sigma^2)\big)$$

**Layer 2: Observation Model**
$$Y_j \sim \text{Binomial}(N, p_j), \quad \text{logit}(p_j) = \alpha + f(s_j)$$

**Layer 3: Priors**
- PC priors on $\rho$ (range) and $\sigma$ (marginal SD)
- Weak Gaussian prior on $\alpha$

<!-- notes -->
Here's our model. The true elevation is a Gaussian random field with Matérn covariance. We observe water frequency through a binomial model - higher elevations mean lower inundation probability. The logit link maps depth to probability.

---

# Why SPDE Instead of Covariance Matrix?

**Problem**: GRF requires $\Sigma \in \mathbb{R}^{n \times n}$
- Inversion is $O(n^3)$
- Storage is $O(n^2)$
- For $n = 50{,}000$ pixels: impossible

**Solution**: SPDE representation (Lindgren et al., 2011)

$$(\kappa^2 - \Delta)^{\alpha/2} Z(s) = \mathcal{W}(s)$$

This yields **sparse** precision matrix $\mathbf{Q}$:
- Only $O(n)$ non-zeros
- Inversion via Cholesky: $O(n^{3/2})$

<!-- notes -->
This is the key computational trick. Instead of working with a dense covariance matrix, we use a stochastic PDE whose solution has Matérn covariance. Discretizing this PDE on a mesh gives a sparse precision matrix. This is what makes the whole project feasible.

---

# Mesh Construction

**Finite Element Method** on triangular mesh:
- Basis functions $\phi_k(s)$: piecewise linear
- Field approximation: $Z(s) \approx \sum_k \phi_k(s) Z_k$

**Our mesh parameters**:
- Inner edge: 100m (high resolution in lake)
- Outer edge: 500m (coarse outside)
- 1,523 vertices

{{Mesh Plot}}

<!-- notes -->
We discretize the SPDE on a triangular mesh. Finer triangles inside the lake capture detail, coarser ones outside save computation. The mesh has about 1,500 nodes - much smaller than 50,000 pixels. This is the dimension reduction.

---

# Projection Matrix: Connecting Data to Mesh

**Observation locations** ≠ **Mesh nodes**

We need projection matrix $\mathbf{A}$:
$$Z_{\text{obs}} = \mathbf{A} \cdot Z_{\text{mesh}}$$

- $\mathbf{A}_{ij} = \phi_j(s_i)$ (barycentric coordinates)
- Each row has at most 3 non-zeros
- Sparse matrix: efficient computation

**Why this matters**: We can observe anywhere, predict anywhere, using the same mesh representation.

<!-- notes -->
The projection matrix A connects observations to mesh nodes. Each observation point lies in some triangle, and we interpolate using barycentric coordinates. This is beautifully sparse - each row has at most three non-zeros.

---

# INLA: Fast Bayesian Inference

**Integrated Nested Laplace Approximation** (Rue et al., 2009)

Instead of MCMC:
1. Approximate $p(\theta | y)$ on a grid
2. At each $\theta$: Gaussian approximation for $p(Z | \theta, y)$
3. Integrate: $p(Z_i | y) = \int p(Z_i | \theta, y) \, p(\theta | y) \, d\theta$

**Advantages**:
- Deterministic (no MCMC diagnostics)
- Fast: minutes instead of hours
- Accurate for latent Gaussian models

<!-- notes -->
INLA is our inference engine. Instead of sampling, it uses Laplace approximations. For our model, this takes about 5 minutes compared to potentially hours with MCMC. And we don't need to worry about convergence diagnostics.

---

# The Calibration Problem

**Issue**: Model estimates *relative* depth $f(s)$, not elevation $Z(s)$

The model is **not identifiable**:
$$\text{logit}(p) = \alpha + f(s) = (\alpha + c) + (f(s) - c)$$

**Solution**: Use shoreline DEM as anchor
- At shore: $Z_{\text{DEM}}(s_i) \approx a + b \cdot f(s_i)$
- Estimate $(a, b)$ from data
- Apply: $Z(s) = a + b \cdot f(s)$

<!-- notes -->
Here's a subtle issue. The model only tells us relative depth - it can't distinguish "intercept" from "field mean". We need external elevation data to anchor our predictions. The shoreline DEM provides these anchor points.

---

# Calibration Framework: 4 Methods

| Method | MAE (km²) | Approach |
|--------|-----------|----------|
| Regression | 23.40 | Robust linear fit on shore |
| Endpoint | 4.23 | Force min/max to match truth |
| Hybrid | 4.23 | Weighted blend |
| **Piecewise** | **1.45** | Segment-wise affine |

**Piecewise calibration** (4 segments):
$$Z^{(k)} = a_k + b_k \cdot f, \quad k = 1, \ldots, 4$$

Captures nonlinear relationship between model depth and true elevation.

<!-- notes -->
We implemented four calibration methods and compared them by MAE. Simple regression doesn't work well - it gives 23 km² error. Piecewise calibration with four segments achieves 1.45 km² error by allowing different slopes in different depth ranges.

---

# Result: Reconstructed Bathymetry

{{Bathymetry Posterior Mean}}

**Key observations**:
- Deepest regions (blue): near dam
- Shallowest (yellow): near shoreline
- Elevation range: 147 – 181 m
- Posterior uncertainty higher in deep water

{{Bathymetry Posterior SD}}

<!-- notes -->
Here's our reconstructed bathymetry. Blue is deep, yellow is shallow. The deepest part is near the dam as expected. Notice uncertainty is higher in deep water - we have less information there because those pixels are always underwater.

---

# Result: Area-Elevation Curve

{{AE Curve Comparison - True vs Predicted}}

**Validation against TWDB survey data**:
- MAE = **1.45 km²**
- R² = **0.998**
- Elevation range: 147.1 – 181.2 m
- Area range: 0 – 50.1 km²

The curves are nearly indistinguishable!

<!-- notes -->
This is the main result. Black line is ground truth from actual surveys, blue is our prediction. They're almost on top of each other. MAE of 1.45 km² out of 50 km² total - that's about 3% error. R-squared is 0.998.

---

# Model Diagnostics

**Posterior hyperparameters**:
- Spatial range: **1,238 m** (correlation decays over ~1 km)
- Marginal SD: **29.9 m** (substantial depth variation)

**Calibration quality**:
- 4-segment piecewise fit
- Each segment adjusts scale and shift
- Residuals show no spatial pattern

{{Calibration Comparison Plot - 4 methods}}

<!-- notes -->
The estimated spatial range is about 1.2 km, which makes physical sense for a lake of this size. The marginal standard deviation of 30 meters captures the depth variation from shore to dam.

---

# What We Contributed

1. **Data pipeline**: Automated preprocessing of DEM, water frequency, masks
   - CRS alignment, resolution matching, outlier removal

2. **Mesh construction**: Adaptive triangulation balancing accuracy vs. speed

3. **Calibration framework**: 4 methods with automatic selection
   - Piecewise affine captures nonlinear depth-elevation relationship

4. **Validation workflow**: Direct comparison with survey data
   - Metrics: MAE, RMSE, R²

5. **Reproducible code**: Modular R pipeline (~2,000 lines)

<!-- notes -->
Let me highlight our contributions. We built the entire pipeline from raw data to final curves. The calibration framework with four methods and automatic selection is novel. And everything is reproducible.

---

# Limitations

1. **Requires validation data for calibration**
   - Without ground truth, we can't calibrate properly
   - Future: use dam elevation as single anchor point?

2. **Stationarity assumption**
   - Assumes same spatial correlation everywhere
   - Reality: rivers, channels may have different structure

3. **Water frequency noise**
   - Clouds, seasonal vegetation affect observations
   - Median filter helps but doesn't eliminate

4. **Single time snapshot**
   - Assumes static bathymetry
   - Ignores sedimentation over decades

<!-- notes -->
We should be honest about limitations. The biggest one: we need some validation data to calibrate. The stationarity assumption is also restrictive - a river channel probably has different spatial structure than open water.

---

# Future Work

1. **Semi-supervised calibration**
   - Use only dam point + physical constraints
   - Remove need for full validation curves

2. **Non-stationary SPDE**
   - Barrier models for land boundaries
   - Spatially varying range parameter

3. **Multi-lake transfer learning**
   - Share hyperparameter priors across similar lakes
   - Hierarchical model over lake population

4. **Temporal extension**
   - Track bathymetry changes over decades
   - Detect sedimentation trends

<!-- notes -->
Several directions for future work. Most exciting is removing the need for validation data entirely - using only physical constraints like the dam elevation. Also, non-stationary models could capture how spatial correlation differs between deep and shallow regions.

---

# Summary

**We showed**:
- INLA-SPDE enables efficient bathymetry reconstruction
- Water occurrence + DEM → underwater topography
- Piecewise calibration reduces MAE by **94%** vs. regression
- Final error: **1.45 km²** (R² = 0.998)

**Key takeaways**:
1. SPDE makes large-scale spatial modeling tractable
2. Calibration is essential for identifiability
3. Validation against ground truth confirms accuracy

---

# Thank You!

**Questions?**

---

**Code**: Available upon request

**Data**: Google Earth Engine, USGS NED, TWDB

**References**:
- Rue et al. (2009). INLA. *JRSSB*.
- Lindgren et al. (2011). SPDE-GMRF link. *JRSSB*.
- Pekel et al. (2016). Global Surface Water. *Nature*.

<!-- notes -->
That's all I have. Happy to take questions about the model, the calibration, or anything else.

