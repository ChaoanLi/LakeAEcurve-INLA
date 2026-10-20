# Lake Bathymetry Reconstruction via INLA-SPDE

**Course**: STAT 647 Spatial Statistics  
**Authors**: Chaoan Li, Yinuo Zhu  
**Institution**: Texas A&M University  
**Date**: December 2025

---

## Overview

This project implements a Bayesian spatial model for reconstructing lake bathymetry and estimating area-elevation (A-E) curves from multi-source remote sensing data. The methodology uses the INLA-SPDE framework to efficiently estimate continuous depth fields from discrete water occurrence observations.

### Key Results (Belton Lake)

| Metric | Value |
|--------|-------|
| MAE | 1.1 km² |
| R² | 0.98 |
| Calibration Method | DEM+4pt (piecewise) |

## Repository Structure

```
LakeAEcurve-INLA/
├── 01_Data_Prep/               # Data preprocessing
│   ├── data/                   # Input data (DEM, water frequency, masks)
│   ├── data_generation.R       # Raster loading and observation building
│   └── mesh_setup.R            # SPDE mesh construction
├── 02_Model_Implementation/    # Model fitting
│   ├── spde_definition.R       # SPDE prior specification
│   └── fit_inlabru_model.R     # INLA stack and model fitting
├── 03_Result_Analysis/         # Post-processing
│   ├── reconstruct_map.R       # Bathymetry reconstruction
│   ├── ae_curve_ppd.R          # A-E curve computation
│   └── calibration_module.R    # Calibration framework (4 methods)
├── 04_Validation/              # Ground truth A-E curves
├── outputs/                    # Generated outputs
├── cache/                      # Cached intermediate results
├── main_runner.R               # Main execution script
└── COMPREHENSIVE_TECHNICAL_DOCUMENTATION.md
```

## Quick Start

### Requirements

```r
install.packages(c("terra", "sf", "ggplot2", "viridis", "patchwork"))

# INLA (from special repository)
install.packages("INLA", 
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
```

### Data and working directory

Run from the repository root: `main_runner.R` uses `getwd()` to locate the modules and input data. Its default Belton Lake configuration expects:

- `01_Data_Prep/data/DEM/Belton_DEM_100m_Buffer.tif`
- `01_Data_Prep/data/occurrence/Belton_GSW_Occurrence_Buffer100m.tif`
- `01_Data_Prep/data/permanent_water/Belton_Min_Permanent_Water.shp`
- `01_Data_Prep/data/DAM_point/Belton_dam.shp`
- `04_Validation/Belton_AVE.csv`

Keep the companion files required by each shapefile together. Check these paths before running; update the `lake_info` configuration if your inputs are stored elsewhere.

### Running the Analysis

To check the observation-building code without running the full spatial fit:

```sh
Rscript tests/observation-data.R
```

These checks cover a dam point supplied in a different coordinate system and
permanent-water cells both present in and absent from the occurrence data.
After updating observation-building code, remove the generated
`cache/Belton_obs_data.RData` before recomputing observations; an existing cache
bypasses that code.

The full fit also requires a working INLA native executable, not just a successful
`library(INLA)` call. The validation environment (conda R 4.5.3 on macOS ARM)
reached native fitting but the official macOS binary expected framework R's
`libRblas.dylib`. Full-model reproduction was therefore not verified in that
environment. Use an R/INLA installation with compatible native libraries.

```r
source("main_runner.R")
```

The script will:
1. Load and preprocess raster data
2. Build SPDE mesh
3. Fit INLA model with binomial likelihood
4. Run calibration framework (4 methods, auto-select best)
5. Generate bathymetry maps and A-E curves

### Output Files

- `outputs/Belton_bathymetry_mean.png` - Posterior mean elevation
- `outputs/Belton_bathymetry_sd.png` - Posterior standard deviation
- `outputs/Belton_ae_curve.png` - A-E curve with validation
- `outputs/Belton_calibration_comparison.png` - Calibration method comparison

## Methodology

### Model

Water occurrence frequency $p(s)$ is modeled via:

$$\text{logit}(p(s)) = \alpha + f(s)$$

where $f(s)$ is a Gaussian random field with Matérn covariance, approximated via SPDE.

### Calibration

Four DEM-based calibration methods are implemented:

| Method | Val. Points | Description |
|--------|-------------|-------------|
| DEM-only | 0 | Extrapolate using DEM slope |
| DEM+1pt | 1 | DEM + dam point anchor |
| DEM+2pt | 2 | DEM + 2-segment interpolation |
| DEM+4pt | 4 | DEM + 4-segment interpolation (best) |

All methods share the same DEM-based calibration above the waterline; differences occur only in underwater extrapolation.

## References

- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models. *JRSSB*.
- Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields. *JRSSB*.
- Simpson, D., et al. (2017). Penalising model component complexity. *Statistical Science*.

## License

MIT License
