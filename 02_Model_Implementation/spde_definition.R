################################################################################
# SPDE Model Definition and Prior Specification
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' Define SPDE model with PC priors
#' @param mesh INLA mesh object
#' @param prior_range Prior range parameters: c(range0, P(range < range0))
#' @param prior_sigma Prior sigma parameters: c(sigma0, P(sigma > sigma0))
#' @return SPDE object
define_spde <- function(mesh,
                        prior_range = c(500, 0.5),
                        prior_sigma = c(1, 0.01)) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Defining SPDE model with PC priors...\n")
  cat("====================================\n")
  
  cat(sprintf("   Prior for range: P(range < %.1f) = %.2f\n", 
              prior_range[1], prior_range[2]))
  cat(sprintf("   Prior for sigma: P(sigma > %.2f) = %.2f\n", 
              prior_sigma[1], prior_sigma[2]))
  
  spde <- inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = prior_range,
    prior.sigma = prior_sigma
  )
  
  cat("\nSPDE model defined.\n\n")
  
  return(spde)
}


#' Build projection matrices from mesh to observation locations
#' @param mesh INLA mesh object
#' @param obs_data List from build_observation_data()
#' @return List with A_elev and A_water matrices
build_projection_matrices <- function(mesh, obs_data) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Building projection matrices...\n")
  cat("====================================\n")
  
  elev_df <- obs_data$elev_df
  water_df <- obs_data$water_df
  
  # 1. Projection matrix for elevation observations
  if (nrow(elev_df) > 0) {
    cat(sprintf("\n1. Building A matrix for elevation observations (%d points)...\n", 
                nrow(elev_df)))
    
    A_elev <- inla.spde.make.A(
      mesh = mesh,
      loc = as.matrix(elev_df[, c("x", "y")])
    )
    
    cat(sprintf("   A_elev dimensions: %d x %d\n", nrow(A_elev), ncol(A_elev)))
  } else {
    cat("\n1. No elevation observations, A_elev set to NULL.\n")
    A_elev <- NULL
  }
  
  # 2. Projection matrix for water frequency observations
  cat(sprintf("\n2. Building A matrix for water observations (%d points)...\n", 
              nrow(water_df)))
  
  A_water <- inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(water_df[, c("x", "y")])
  )
  
  cat(sprintf("   A_water dimensions: %d x %d\n", nrow(A_water), ncol(A_water)))
  
  cat("\nProjection matrices built.\n\n")
  
  return(list(
    A_elev = A_elev,
    A_water = A_water
  ))
}
