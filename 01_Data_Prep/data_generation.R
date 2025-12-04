################################################################################
# Data Loading and Preprocessing Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' Load and preprocess raster data for one lake
#' @param dem_path Path to shoreline DEM raster
#' @param water_freq_path Path to water occurrence frequency raster
#' @param perm_water_path Path to permanent water shapefile
#' @param lake_name Lake name (for output messages)
#' @param use_cache Whether to use cache (default TRUE)
#' @param cache_dir Cache directory (default "cache")
#' @return List containing all preprocessed data
load_and_prep_data <- function(dem_path, 
                               water_freq_path, 
                               perm_water_path,
                               lake_name = "Lake",
                               use_cache = TRUE,
                               cache_dir = "cache") {
  
  library(terra)
  library(sf)
  
  # Check cache
  if (use_cache) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    
    cache_file <- file.path(cache_dir, sprintf("%s_prep_data.RData", lake_name))
    
    if (file.exists(cache_file)) {
      cat("====================================\n")
      cat(sprintf("Loading cached data for %s...\n", lake_name))
      cat("====================================\n")
      cat(sprintf("   Cache file: %s\n", cache_file))
      
      load(cache_file)
      
      if ("dem_file" %in% names(cached_data)) {
        cat("   Loading raster files from cache...\n")
        result_data <- list(
          dem = rast(cached_data$dem_file),
          water_freq = rast(cached_data$water_freq_file),
          perm_water = rast(cached_data$perm_water_file),
          lake_boundary = cached_data$lake_boundary,
          perm_water_polygon = cached_data$perm_water_polygon,
          cell_area = cached_data$cell_area,
          lake_name = cached_data$lake_name
        )
        cat("Cached data loaded successfully.\n\n")
        return(result_data)
      } else {
        cat("Old cache format detected, reprocessing...\n\n")
      }
    }
  }
  
  cat("====================================\n")
  cat(sprintf("Loading data for %s...\n", lake_name))
  cat("====================================\n")
  
  # 1. Read raster data
  cat("1. Reading raster data...\n")
  
  dem_rast <- rast(dem_path)
  cat(sprintf("   DEM loaded: %s\n", dem_path))
  
  water_freq_rast <- rast(water_freq_path)
  cat(sprintf("   Water frequency loaded: %s\n", water_freq_path))
  
  perm_water_sf <- st_read(perm_water_path, quiet = TRUE)
  cat(sprintf("   Permanent water mask loaded: %s\n", perm_water_path))
  
  # 2. Check CRS and resolution alignment
  cat("\n2. Checking CRS and resolution alignment...\n")
  
  dem_crs <- crs(dem_rast)
  water_freq_crs <- crs(water_freq_rast)
  perm_water_crs <- st_crs(perm_water_sf)
  
  cat(sprintf("   Original DEM CRS: EPSG:%s\n", 
              st_crs(dem_crs)$epsg %||% "unknown"))
  
  is_geographic <- st_is_longlat(st_crs(dem_crs))
  
  if (is_geographic) {
    cat("   Geographic CRS detected (degrees)\n")
    cat("   Projecting to UTM (meters)...\n")
    
    bbox <- st_bbox(st_as_sfc(st_bbox(dem_rast)))
    center_lon <- mean(c(bbox["xmin"], bbox["xmax"]))
    center_lat <- mean(c(bbox["ymin"], bbox["ymax"]))
    
    utm_zone <- floor((center_lon + 180) / 6) + 1
    utm_epsg <- ifelse(center_lat >= 0, 
                       32600 + utm_zone,
                       32700 + utm_zone)
    
    cat(sprintf("   Target CRS: UTM Zone %d (EPSG:%d)\n", utm_zone, utm_epsg))
    
    dem_rast <- project(dem_rast, paste0("EPSG:", utm_epsg), method = "bilinear")
    water_freq_rast <- project(water_freq_rast, paste0("EPSG:", utm_epsg), 
                               method = "bilinear")
    perm_water_sf <- st_transform(perm_water_sf, crs = utm_epsg)
    
    dem_crs <- crs(dem_rast)
    cat("   Projection complete.\n")
  } else {
    cat("   CRS is already projected.\n")
  }
  
  if (dem_crs != crs(water_freq_rast)) {
    cat("   Reprojecting water frequency to DEM CRS...\n")
    water_freq_rast <- project(water_freq_rast, dem_rast, method = "bilinear")
  }
  
  if (!st_crs(perm_water_sf) == st_crs(dem_crs)) {
    cat("   Reprojecting permanent water to DEM CRS...\n")
    perm_water_sf <- st_transform(perm_water_sf, crs = dem_crs)
  }
  
  dem_res <- res(dem_rast)
  water_freq_res <- res(water_freq_rast)
  
  cat(sprintf("   DEM resolution: %.2f x %.2f\n", dem_res[1], dem_res[2]))
  cat(sprintf("   Water frequency resolution: %.2f x %.2f\n", 
              water_freq_res[1], water_freq_res[2]))
  
  if (!all(abs(dem_res - water_freq_res) < 1e-6)) {
    cat("   Resampling water frequency to DEM resolution...\n")
    water_freq_rast <- resample(water_freq_rast, dem_rast, method = "bilinear")
  }
  
  # 3. Create permanent water raster
  cat("\n3. Creating permanent water raster mask...\n")
  
  perm_water_rast <- rasterize(perm_water_sf, dem_rast, field = 1, 
                               background = 0, touches = FALSE)
  cat("   Permanent water rasterized.\n")
  
  # 4. Normalize water frequency to [0, 1]
  cat("\n4. Normalizing water frequency to [0, 1]...\n")
  
  water_freq_max <- global(water_freq_rast, "max", na.rm = TRUE)[[1]]
  cat(sprintf("   Water frequency max value: %.2f\n", water_freq_max))
  
  if (water_freq_max > 1.5) {
    cat("   Converting from [0, 100] to [0, 1]...\n")
    water_freq_rast <- water_freq_rast / 100
  }
  
  water_freq_rast <- clamp(water_freq_rast, lower = 0, upper = 1, values = TRUE)
  
  # 5. Create lake boundary polygon
  cat("\n5. Creating lake boundary polygon...\n")
  
  lake_boundary <- st_union(perm_water_sf)
  lake_boundary <- st_buffer(lake_boundary, dist = 500)
  lake_boundary <- st_as_sf(data.frame(id = 1, geom = lake_boundary))
  st_geometry(lake_boundary) <- "geom"
  
  cat("   Lake boundary created with 500m buffer.\n")
  
  # 6. Mask rasters to lake boundary
  cat("\n6. Masking rasters to lake boundary...\n")
  
  dem_masked <- crop(dem_rast, lake_boundary)
  dem_masked <- mask(dem_masked, lake_boundary)
  
  water_freq_masked <- crop(water_freq_rast, lake_boundary)
  water_freq_masked <- mask(water_freq_masked, lake_boundary)
  
  perm_water_masked <- crop(perm_water_rast, lake_boundary)
  perm_water_masked <- mask(perm_water_masked, lake_boundary)
  
  cat("   Rasters masked to lake boundary.\n")
  
  # 7. Summary
  cat("\n7. Data summary:\n")
  cat(sprintf("   DEM range: %.2f to %.2f m\n", 
              global(dem_masked, "min", na.rm = TRUE)[[1]],
              global(dem_masked, "max", na.rm = TRUE)[[1]]))
  cat(sprintf("   Water frequency range: %.3f to %.3f\n", 
              global(water_freq_masked, "min", na.rm = TRUE)[[1]],
              global(water_freq_masked, "max", na.rm = TRUE)[[1]]))
  cat(sprintf("   Permanent water pixels: %d\n", 
              sum(values(perm_water_masked) == 1, na.rm = TRUE)))
  
  cell_area <- prod(res(dem_masked))
  cat(sprintf("   Cell area: %.2f m2\n", cell_area))
  
  cat("\nData loading and preprocessing complete.\n\n")
  
  # 8. Return result
  result_data <- list(
    dem = dem_masked,
    water_freq = water_freq_masked,
    perm_water = perm_water_masked,
    lake_boundary = lake_boundary,
    perm_water_polygon = perm_water_sf,
    cell_area = cell_area,
    lake_name = lake_name
  )
  
  # 9. Save cache
  if (use_cache) {
    cat("Saving preprocessed data to cache...\n")
    
    dem_cache_file <- file.path(cache_dir, sprintf("%s_dem_cache.tif", lake_name))
    water_freq_cache_file <- file.path(cache_dir, sprintf("%s_water_freq_cache.tif", lake_name))
    perm_water_cache_file <- file.path(cache_dir, sprintf("%s_perm_water_cache.tif", lake_name))
    
    writeRaster(dem_masked, dem_cache_file, overwrite = TRUE)
    writeRaster(water_freq_masked, water_freq_cache_file, overwrite = TRUE)
    writeRaster(perm_water_masked, perm_water_cache_file, overwrite = TRUE)
    
    cached_data <- list(
      dem_file = dem_cache_file,
      water_freq_file = water_freq_cache_file,
      perm_water_file = perm_water_cache_file,
      lake_boundary = lake_boundary,
      perm_water_polygon = perm_water_sf,
      cell_area = cell_area,
      lake_name = lake_name
    )
    
    save(cached_data, file = cache_file)
    cat(sprintf("Cache saved: %s\n\n", cache_file))
  }
  
  return(result_data)
}


#' Build observation data frames for INLA
#' @param data_list List from load_and_prep_data()
#' @param N_trials Binomial trial count (default 100)
#' @param use_shore_elev Whether to use shore elevation as observations (default TRUE)
#' @param shore_buffer_m Shore buffer width in meters
#' @param dam_point_path Dam point shapefile path (optional)
#' @param dam_point_weight Dam point replication count (default 50)
#' @param use_cache Whether to use cache (default TRUE)
#' @param cache_dir Cache directory (default "cache")
#' @return List with elev_df and water_df
build_observation_data <- function(data_list, 
                                   N_trials = 100,
                                   use_shore_elev = TRUE,
                                   shore_buffer_m = 200,
                                   dam_point_path = NULL,
                                   dam_point_weight = 50,
                                   use_cache = TRUE,
                                   cache_dir = "cache") {
  
  library(terra)
  
  # Check cache
  if (use_cache) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    
    lake_name <- data_list$lake_name
    cache_file <- file.path(cache_dir, sprintf("%s_obs_data.RData", lake_name))
    
    if (file.exists(cache_file)) {
      cat("====================================\n")
      cat("Loading cached observation data...\n")
      cat("====================================\n")
      cat(sprintf("   Cache file: %s\n", cache_file))
      
      load(cache_file)
      
      cat("Cached observation data loaded successfully.\n\n")
      
      return(cached_obs_data)
    }
  }
  
  cat("====================================\n")
  cat("Building observation data frames...\n")
  cat("====================================\n")
  
  dem <- data_list$dem
  water_freq <- data_list$water_freq
  perm_water <- data_list$perm_water
  
  # 1. Build elevation observation data frame
  cat("\n1. Building elevation observation data...\n")
  
  if (use_shore_elev) {
    # Shore mask: DEM valid and not permanent water
    shore_mask <- !is.na(dem) & (is.na(perm_water) | perm_water == 0)
    
    elev_coords_all <- xyFromCell(dem, which(values(shore_mask)))
    elev_values_all <- values(dem)[values(shore_mask)]
    
    valid_idx <- !is.na(elev_values_all)
    elev_coords_all <- elev_coords_all[valid_idx, ]
    elev_values_all <- elev_values_all[valid_idx]
    
    # 3-sigma outlier removal
    if (length(elev_values_all) > 10) {
      z_scores <- abs(scale(elev_values_all))
      keep_idx <- z_scores < 3
      
      n_removed <- sum(!keep_idx)
      if (n_removed > 0) {
        cat(sprintf("   Removing %d outliers (%.1f%%) using 3-sigma criterion\n",
                    n_removed, 100 * n_removed / length(elev_values_all)))
      }
      
      elev_coords <- elev_coords_all[keep_idx, ]
      elev_values <- elev_values_all[keep_idx]
    } else {
      elev_coords <- elev_coords_all
      elev_values <- elev_values_all
    }
    
    elev_df <- data.frame(
      x = elev_coords[, 1],
      y = elev_coords[, 2],
      elev = elev_values
    )
    
    cat(sprintf("   Total shore points: %d\n", length(elev_values_all)))
    cat(sprintf("   Shore points after 3-sigma filtering: %d\n", nrow(elev_df)))
    cat(sprintf("   Elevation range: %.2f to %.2f m\n", 
                min(elev_df$elev), max(elev_df$elev)))
  } else {
    elev_df <- data.frame(x = numeric(0), y = numeric(0), elev = numeric(0))
    cat("   No elevation observations used.\n")
  }
  
  # 1.5. Add dam point as lowest elevation constraint
  if (!is.null(dam_point_path) && file.exists(dam_point_path)) {
    cat("\n   Adding dam point as lowest elevation constraint...\n")
    
    tryCatch({
      library(sf)
      
      dam_point <- st_read(dam_point_path, quiet = TRUE)
      
      dam_crs <- st_crs(dam_point)
      dem_crs <- crs(dem)
      
      if (!st_is_longlat(dam_crs) && st_is_longlat(st_crs(dem_crs))) {
        cat("   Transforming dam point CRS to match DEM...\n")
        dam_point <- st_transform(dam_point, crs = dem_crs)
      }
      
      dam_coords <- st_coordinates(dam_point)
      
      cat(sprintf("   Dam point location: (%.6f, %.6f)\n", 
                  dam_coords[1], dam_coords[2]))
      
      dam_elev_raw <- extract(dem, dam_coords)
      
      if (is.data.frame(dam_elev_raw)) {
        dam_elevation <- dam_elev_raw[1, ncol(dam_elev_raw)]
      } else if (is.matrix(dam_elev_raw)) {
        dam_elevation <- dam_elev_raw[1, ncol(dam_elev_raw)]
      } else {
        dam_elevation <- as.numeric(dam_elev_raw)[1]
      }
      
      if (length(dam_elevation) > 0 && !is.na(dam_elevation)) {
        cat(sprintf("   Dam elevation from DEM: %.2f m\n", dam_elevation))
        cat(sprintf("   Replicating %d times to increase constraint weight\n", 
                    dam_point_weight))
        
        dam_obs <- data.frame(
          x = rep(dam_coords[1], dam_point_weight),
          y = rep(dam_coords[2], dam_point_weight),
          elev = rep(dam_elevation, dam_point_weight)
        )
        
        elev_df <- rbind(elev_df, dam_obs)
        
        cat(sprintf("   Total elevation observations: %d (including %d dam replicates)\n",
                    nrow(elev_df), dam_point_weight))
      } else {
        cat("   Warning: Could not extract valid elevation at dam point location\n")
        cat("   Continuing without dam constraint...\n")
      }
    }, error = function(e) {
      cat(sprintf("   Warning: Error processing dam point: %s\n", e$message))
      cat("   Continuing without dam constraint...\n")
    })
  }
  
  # 2. Build water frequency observation data frame
  cat("\n2. Building water frequency observation data...\n")
  
  water_valid <- !is.na(water_freq)
  water_coords <- xyFromCell(water_freq, which(values(water_valid)))
  water_p <- values(water_freq)[values(water_valid)]
  
  valid_idx <- !is.na(water_p)
  water_coords <- water_coords[valid_idx, ]
  water_p <- water_p[valid_idx]
  
  # For permanent water pixels, set p = 1
  perm_water_cells <- which(values(perm_water) == 1)
  perm_water_coords <- xyFromCell(perm_water, perm_water_cells)
  
  water_df <- data.frame(
    x = water_coords[, 1],
    y = water_coords[, 2],
    p = water_p
  )
  
  # Add permanent water points (p = 1)
  if (length(perm_water_cells) > 0) {
    perm_df <- data.frame(
      x = perm_water_coords[, 1],
      y = perm_water_coords[, 2],
      p = 1.0
    )
    
    for (i in 1:nrow(perm_df)) {
      match_idx <- which(water_df$x == perm_df$x[i] & water_df$y == perm_df$y[i])
      if (length(match_idx) > 0) {
        water_df$p[match_idx] <- 1.0
      } else {
        water_df <- rbind(water_df, perm_df[i, ])
      }
    }
  }
  
  # Convert to binomial counts
  water_df$N <- N_trials
  water_df$y_water <- round(water_df$p * water_df$N)
  
  water_df$y_water <- pmax(0, pmin(water_df$N, water_df$y_water))
  
  cat(sprintf("   Water frequency observations: %d points\n", nrow(water_df)))
  cat(sprintf("   Water occurrence range: %.3f to %.3f\n", 
              min(water_df$p), max(water_df$p)))
  cat(sprintf("   Binomial counts (y_water): %d to %d (N = %d)\n", 
              min(water_df$y_water), max(water_df$y_water), N_trials))
  
  n_permanent <- sum(water_df$y_water == water_df$N)
  cat(sprintf("   Permanent water points (y = N): %d\n", n_permanent))
  
  cat("\nObservation data frames built.\n\n")
  
  result_obs <- list(
    elev_df = elev_df,
    water_df = water_df,
    N_trials = N_trials
  )
  
  # Save cache
  if (use_cache) {
    cat("Saving observation data to cache...\n")
    cached_obs_data <- result_obs
    save(cached_obs_data, file = cache_file)
    cat(sprintf("Cache saved: %s\n\n", cache_file))
  }
  
  return(result_obs)
}
