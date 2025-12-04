################################################################################
# Lake Area-Elevation Curve Estimation via INLA-SPDE
#
# Authors: Chaoan Li, Yinuo Zhu
# Course: STAT 647 Spatial Statistics, Texas A&M University
# Date: December 2025
################################################################################

# ============================
# 0. Environment Setup
# ============================

cat("\n")
cat("================================================================================\n")
cat("  Lake Bathymetry Reconstruction using INLA-SPDE\n")
cat("================================================================================\n")
cat("\n")

rm(list = ls())
gc()

required_packages <- c("terra", "sf", "ggplot2", "viridis", "patchwork", "INLA")

cat("Loading required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("   Installing: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("All packages loaded.\n\n")

if (!require("INLA", quietly = TRUE)) {
  cat("Installing INLA...\n")
  install.packages("INLA", 
                   repos = c(getOption("repos"), 
                             INLA = "https://inla.r-inla-download.org/R/stable"), 
                   dep = TRUE)
  library(INLA)
}

# ============================
# 1. Configuration
# ============================

cat("================================================================================\n")
cat("1. Configuration\n")
cat("================================================================================\n\n")

project_dir <- getwd()
cat(sprintf("Project directory: %s\n\n", project_dir))

source(file.path(project_dir, "01_Data_Prep", "data_generation.R"))
source(file.path(project_dir, "01_Data_Prep", "mesh_setup.R"))
source(file.path(project_dir, "02_Model_Implementation", "spde_definition.R"))
source(file.path(project_dir, "02_Model_Implementation", "fit_inlabru_model.R"))
source(file.path(project_dir, "03_Result_Analysis", "reconstruct_map.R"))
source(file.path(project_dir, "03_Result_Analysis", "ae_curve_ppd.R"))

cat("All modules loaded.\n\n")

data_dir <- file.path(project_dir, "01_Data_Prep", "data")

output_dir <- file.path(project_dir, "outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================
# 2. Lake Configuration
# ============================

cat("================================================================================\n")
cat("2. Belton Lake Configuration\n")
cat("================================================================================\n\n")

lake_info <- list(
  name = "Belton",
  dem_path = file.path(data_dir, "DEM", "Belton_DEM_100m_Buffer.tif"),
  water_freq_path = file.path(data_dir, "occurrence", 
                              "Belton_GSW_Occurrence_Buffer100m.tif"),
  perm_water_path = file.path(data_dir, "permanent_water", 
                              "Belton_Min_Permanent_Water.shp"),
  dam_point_path = file.path(data_dir, "DAM_point", "Belton_dam.shp"),
  true_ae_path = "04_Validation/Belton_AVE.csv"
)

lakes_to_process <- list(lake_info)

cat("Processing: Belton Lake\n\n")

# ============================
# 3. Main Processing Loop
# ============================

cat("================================================================================\n")
cat("3. Processing\n")
cat("================================================================================\n\n")

all_results <- list()

for (lake_idx in seq_along(lakes_to_process)) {
  
  lake_info <- lakes_to_process[[lake_idx]]
  lake_name <- lake_info$name
  
  cat("\n")
  cat("################################################################################\n")
  cat(sprintf("# Processing: %s (%d/%d)\n", lake_name, lake_idx, length(lakes_to_process)))
  cat("################################################################################\n")
  cat("\n")
  
  # 3.1 Load and preprocess data
  data_list <- load_and_prep_data(
    dem_path = lake_info$dem_path,
    water_freq_path = lake_info$water_freq_path,
    perm_water_path = lake_info$perm_water_path,
    lake_name = lake_name
  )
  
  # 3.2 Build observation data
  obs_data <- build_observation_data(
    data_list = data_list,
    N_trials = 100,
    use_shore_elev = TRUE,
    dam_point_path = lake_info$dam_point_path,
    dam_point_weight = 50
  )
  
  data_list$obs_data <- obs_data
  
  # 3.3 Visualize original data
  plot_original_data(
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_original_data.png", lake_name))
  )
  
  # 3.4 Build SPDE mesh
  mesh <- build_mesh(
    data_list = data_list,
    max_edge = c(100, 500),
    cutoff = 50,
    offset = c(100, 500)
  )
  
  plot_mesh(
    mesh = mesh,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_mesh.png", lake_name))
  )
  
  # 3.5 Define SPDE model
  spde <- define_spde(
    mesh = mesh,
    prior_range = c(500, 0.5),
    prior_sigma = c(1, 0.01)
  )
  
  # 3.6 Build projection matrices
  proj_matrices <- build_projection_matrices(
    mesh = mesh,
    obs_data = obs_data
  )
  
  # 3.7 Build INLA stacks
  stack_list <- build_inla_stacks(
    spde = spde,
    obs_data = obs_data,
    proj_matrices = proj_matrices
  )
  
  # 3.8 Fit INLA model
  result <- fit_inla_model(
    stack_list = stack_list,
    spde = spde,
    obs_data = obs_data,
    use_elev = TRUE,
    beta_prior = c(0, 0.01)
  )
  
  # 3.9 Extract hyperparameters
  hyperpar <- extract_spde_hyperpar(
    result = result,
    spde = spde
  )
  
  # 3.10 Reconstruct bathymetry
  bathy_result <- reconstruct_bathymetry(
    result = result,
    mesh = mesh,
    data_list = data_list,
    template_raster = data_list$dem
  )
  
  plot_bathymetry(
    bathy_result = bathy_result,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_bathymetry.png", lake_name))
  )
  
  # 3.11 Compute A-E curve
  ae_df <- compute_ae_curve(
    bathy_result = bathy_result,
    data_list = data_list,
    elevation_step = 0.1
  )
  
  plot_ae_curve(
    ae_df = ae_df,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_ae_curve.png", lake_name)),
    true_ae_path = lake_info$true_ae_path
  )
  
  # 3.12 Generate calibration comparison plot
  if (!is.null(bathy_result$calibration_results)) {
    cat("\nGenerating calibration comparison plot...\n")
    source("03_Result_Analysis/calibration_module.R")
    
    plot_calibration_comparison(
      calib_results = bathy_result$calibration_results,
      data_list = data_list,
      save_path = file.path(output_dir, sprintf("%s_calibration_comparison.png", lake_name))
    )
    
    calib_report <- generate_calibration_report(bathy_result$calibration_results)
    cat("\n", calib_report, "\n")
  }
  
  # 3.13 Save results
  cat("\nSaving results...\n")
  
  lake_results <- list(
    lake_name = lake_name,
    data_list = data_list,
    obs_data = obs_data,
    mesh = mesh,
    spde = spde,
    result = result,
    hyperpar = hyperpar,
    bathy_result = bathy_result,
    ae_df = ae_df
  )
  
  save(lake_results, 
       file = file.path(output_dir, sprintf("%s_results.RData", lake_name)))
  
  cat(sprintf("Results saved: %s_results.RData\n", lake_name))
  
  write.csv(ae_df, 
            file = file.path(output_dir, sprintf("%s_ae_curve.csv", lake_name)),
            row.names = FALSE)
  
  cat(sprintf("A-E curve saved: %s_ae_curve.csv\n", lake_name))
  
  all_results[[lake_name]] <- lake_results
  
  cat("\n")
  cat(sprintf("Lake %s processing complete.\n", lake_name))
  cat("\n")
}

# ============================
# 4. Summary
# ============================

cat("\n")
cat("================================================================================\n")
cat("4. Summary\n")
cat("================================================================================\n\n")

cat(sprintf("Processed %d lake(s):\n", length(all_results)))
for (lake_name in names(all_results)) {
  cat(sprintf("  - %s\n", lake_name))
}

cat(sprintf("\nOutput directory: %s\n", output_dir))
cat("\nOutput files:\n")
cat("  - *_original_data.png: Input data visualization\n")
cat("  - *_mesh.png: SPDE mesh visualization\n")
cat("  - *_bathymetry_mean.png: Posterior mean elevation\n")
cat("  - *_bathymetry_sd.png: Posterior standard deviation\n")
cat("  - *_ae_curve.png: Area-Elevation curve\n")
cat("  - *_calibration_comparison.png: Calibration method comparison\n")
cat("  - *_ae_curve.csv: A-E curve data (CSV)\n")
cat("  - *_results.RData: Complete results (R objects)\n")

cat("\n")
cat("================================================================================\n")
cat("  Analysis Complete\n")
cat("================================================================================\n")
cat("\n")

# ============================
# 5. Key Results
# ============================

cat("\n")
cat("================================================================================\n")
cat("5. Key Results\n")
cat("================================================================================\n\n")

for (lake_name in names(all_results)) {
  lake_res <- all_results[[lake_name]]
  
  cat(sprintf("\n--- %s ---\n", lake_name))
  
  fixed_summary <- lake_res$result$summary.fixed
  
  if ("intercept_water" %in% rownames(fixed_summary)) {
    alpha <- fixed_summary["intercept_water", "mean"]
    cat(sprintf("  Alpha (intercept): %.4f\n", alpha))
  }
  
  if ("neg_field_water" %in% rownames(fixed_summary)) {
    beta <- fixed_summary["neg_field_water", "mean"]
    cat(sprintf("  Beta (elevation coef): %.4f\n", beta))
  }
  
  cat(sprintf("  Spatial range: %.2f m\n", lake_res$hyperpar$range$mean))
  cat(sprintf("  Spatial sigma: %.2f m\n", lake_res$hyperpar$sigma$mean))
  
  ae_df <- lake_res$ae_df
  cat(sprintf("  A-E curve: %.2f to %.2f km2 (at elev %.2f to %.2f m)\n",
              min(ae_df$area_km2), max(ae_df$area_km2),
              min(ae_df$elevation), max(ae_df$elevation)))
}

cat("\n")
cat("================================================================================\n")
cat("Program finished successfully.\n")
cat("================================================================================\n")
cat("\n")
