################################################################################
# INLA Model Fitting Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' Build INLA stacks for model fitting
#' @param spde SPDE object
#' @param obs_data List from build_observation_data()
#' @param proj_matrices List from build_projection_matrices()
#' @return List with stack objects
build_inla_stacks <- function(spde, obs_data, proj_matrices) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Building INLA stacks...\n")
  cat("====================================\n")
  
  elev_df <- obs_data$elev_df
  water_df <- obs_data$water_df
  A_elev <- proj_matrices$A_elev
  A_water <- proj_matrices$A_water
  
  spde_index <- inla.spde.make.index("field", n.spde = spde$n.spde)
  
  # Using water-only model for stability
  # Shore elevation is used for calibration in post-processing
  
  cat("\n1. Building water frequency stack (Binomial likelihood)...\n")
  cat("   Note: Using water-only model for stability\n")
  cat("   Shore elevation will be used for calibration in post-processing\n")
  
  stack_water <- inla.stack(
    data = list(
      y = water_df$y_water,
      Ntrials = water_df$N
    ),
    A = list(1, A_water),
    effects = list(
      list(intercept = rep(1, nrow(water_df))),
      spde_index
    ),
    tag = "water"
  )
  
  cat(sprintf("   Water observations: %d\n", nrow(water_df)))
  cat(sprintf("   Shore observations available for calibration: %d\n", nrow(elev_df)))
  
  stack_full <- stack_water
  stack_elev <- NULL
  
  cat("\nINLA stacks built.\n\n")
  
  return(list(
    stack_full = stack_full,
    stack_elev = stack_elev,
    stack_water = stack_water
  ))
}


#' Fit INLA model
#' @param stack_list List from build_inla_stacks()
#' @param spde SPDE object
#' @param obs_data List from build_observation_data()
#' @param use_elev Whether to use elevation observations
#' @param beta_prior Prior mean and precision for beta: c(mean, precision)
#' @return INLA result object
fit_inla_model <- function(stack_list, 
                           spde,
                           obs_data,
                           use_elev = TRUE,
                           beta_prior = c(0, 0.01)) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Fitting INLA model...\n")
  cat("====================================\n")
  
  stack_full <- stack_list$stack_full
  
  # 1. Build model formula
  cat("\n1. Building model formula...\n")
  
  # Spatial field interpretation: positive values = deeper water
  
  formula <- y ~ -1 + intercept + f(field, model = spde)
  family_vec <- "binomial"
  
  cat("   Formula: y ~ intercept + f(field)\n")
  cat("   Family: Binomial\n")
  cat("   Field interpretation: relative bathymetry (positive = deeper)\n")
  
  # 2. Set priors
  cat("\n2. Setting priors...\n")
  
  control_fixed <- list(
    mean = list(intercept = 0),
    prec = list(intercept = 0.01)
  )
  
  cat("   Intercept prior: N(0, precision = 0.01)\n")
  
  # 3. Fit model
  cat("\n3. Fitting model with INLA...\n")
  cat("   (This may take several minutes...)\n\n")
  
  stack_data <- inla.stack.data(stack_full)
  stack_A <- inla.stack.A(stack_full)
  
  result <- inla(
    formula = formula,
    data = stack_data,
    family = family_vec,
    Ntrials = stack_data$Ntrials,
    control.predictor = list(
      A = stack_A,
      compute = TRUE
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE
    ),
    control.fixed = control_fixed,
    verbose = TRUE
  )
  
  cat("\nModel fitting complete.\n\n")
  
  # 4. Model summary
  cat("====================================\n")
  cat("Model Summary\n")
  cat("====================================\n\n")
  
  print(summary(result))
  
  cat("\n====================================\n")
  cat("Fixed Effects Summary\n")
  cat("====================================\n")
  
  fixed_summary <- result$summary.fixed
  print(fixed_summary)
  
  if ("intercept" %in% rownames(fixed_summary)) {
    intercept_mean <- fixed_summary["intercept", "mean"]
    intercept_sd <- fixed_summary["intercept", "sd"]
    cat(sprintf("\nIntercept (baseline logit): %.4f +/- %.4f\n", 
                intercept_mean, intercept_sd))
  }
  
  cat("\n====================================\n")
  cat("Hyperparameters (SPDE)\n")
  cat("====================================\n")
  
  hyper_summary <- result$summary.hyperpar
  print(hyper_summary)
  
  if ("Theta1 for field" %in% rownames(hyper_summary)) {
    cat("\n   Note: Theta parameters correspond to log(tau) and log(kappa).\n")
    cat("   Use inla.spde2.result() to extract range and sigma.\n")
  }
  
  cat("\n====================================\n\n")
  
  return(result)
}


#' Extract spatial field hyperparameters (range and sigma)
#' @param result INLA fit result
#' @param spde SPDE object
#' @return List with range and sigma posterior statistics
extract_spde_hyperpar <- function(result, spde) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Extracting SPDE hyperparameters...\n")
  cat("====================================\n\n")
  
  spde_result <- inla.spde2.result(result, "field", spde)
  
  # Range posterior
  tryCatch({
    range_summary <- inla.zmarginal(spde_result$marginals.range.nominal[[1]])
    cat("Range (spatial correlation distance):\n")
    print(range_summary)
    cat("\n")
  }, error = function(e) {
    cat("Warning: Could not extract range summary\n")
    range_summary <<- list(mean = NA, sd = NA)
  })
  
  # Variance posterior
  tryCatch({
    variance_marginal <- spde_result$marginals.variance.nominal[[1]]
    variance_marginal <- variance_marginal[variance_marginal[,1] >= 0, ]
    
    sigma_summary <- inla.zmarginal(variance_marginal)
    cat("Variance (marginal variance of the field):\n")
    print(sigma_summary)
    cat("\n")
  }, error = function(e) {
    cat("Warning: Could not extract variance summary\n")
    sigma_summary <<- list(mean = NA, sd = NA)
  })
  
  # Compute sigma (standard deviation)
  tryCatch({
    variance_marginal <- spde_result$marginals.variance.nominal[[1]]
    variance_marginal <- variance_marginal[variance_marginal[,1] >= 0, ]
    
    sigma_marginal <- inla.tmarginal(
      function(x) sqrt(pmax(x, 0)),
      variance_marginal
    )
    
    sigma_marginal <- sigma_marginal[is.finite(sigma_marginal[,1]) & is.finite(sigma_marginal[,2]), ]
    
    if (nrow(sigma_marginal) > 0) {
      sigma_summary2 <- inla.zmarginal(sigma_marginal)
      cat("Sigma (standard deviation):\n")
      print(sigma_summary2)
      cat("\n")
    } else {
      cat("Warning: Could not compute valid sigma values\n")
      sigma_summary2 <- list(mean = NA, sd = NA)
    }
  }, error = function(e) {
    cat(sprintf("Warning: Error computing sigma: %s\n", e$message))
    sigma_summary2 <- list(mean = NA, sd = NA)
  })
  
  cat("====================================\n\n")
  
  return(list(
    range = range_summary,
    variance = sigma_summary,
    sigma = sigma_summary2,
    spde_result = spde_result
  ))
}
