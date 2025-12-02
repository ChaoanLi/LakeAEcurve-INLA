################################################################################
# Data Loading and Preprocessing Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' 加载和预处理栅格数据
#' Load and preprocess raster data for one lake
#' 
#' @param dem_path 岸边 DEM 栅格路径
#' @param water_freq_path 水出现频率栅格路径
#' @param perm_water_path 永久水域 shapefile 路径
#' @param lake_name 湖泊名称（用于输出信息）
#' @param use_cache 是否使用缓存（默认 TRUE）
#' @param cache_dir 缓存目录（默认 "cache"）
#' @return 包含所有预处理数据的列表
load_and_prep_data <- function(dem_path, 
                                water_freq_path, 
                                perm_water_path,
                                lake_name = "Lake",
                                use_cache = TRUE,
                                cache_dir = "cache") {
  
  library(terra)
  library(sf)
  
  # ---- 检查缓存 ----
  if (use_cache) {
    # 创建缓存目录
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    
    # 缓存文件名
    cache_file <- file.path(cache_dir, sprintf("%s_prep_data.RData", lake_name))
    
    # 如果缓存存在，直接加载
    if (file.exists(cache_file)) {
      cat("====================================\n")
      cat(sprintf("Loading cached data for %s...\n", lake_name))
      cat("====================================\n")
      cat(sprintf("   Cache file: %s\n", cache_file))
      
      load(cache_file)  # 加载 cached_data
      
      # 检查是否是新格式缓存(包含栅格文件路径)
      if ("dem_file" %in% names(cached_data)) {
        # 新格式：从文件读取栅格
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
        cat("✓ Cached data loaded successfully!\n")
        cat(sprintf("   Skipping preprocessing (saved time!)\n\n"))
        return(result_data)
      } else {
        # 旧格式缓存：外部指针可能失效，需要重新处理
        cat("⚠ Old cache format detected (may have invalid pointers)\n")
        cat("   Reprocessing data...\n\n")
        # 继续执行下面的正常处理流程
      }
    }
  }
  
  cat("====================================\n")
  cat(sprintf("Loading data for %s...\n", lake_name))
  cat("====================================\n")
  
  # ---- 1. 读取栅格数据 ----
  cat("1. Reading raster data...\n")
  
  # 读取 DEM（岸边高程）
  dem_rast <- rast(dem_path)
  cat(sprintf("   DEM loaded: %s\n", dem_path))
  
  # 读取水出现频率栅格
  water_freq_rast <- rast(water_freq_path)
  cat(sprintf("   Water frequency loaded: %s\n", water_freq_path))
  
  # 读取永久水域多边形
  perm_water_sf <- st_read(perm_water_path, quiet = TRUE)
  cat(sprintf("   Permanent water mask loaded: %s\n", perm_water_path))
  
  # ---- 2. 检查 CRS 和分辨率一致性 ----
  cat("\n2. Checking CRS and resolution alignment...\n")
  
  # 获取 CRS
  dem_crs <- crs(dem_rast)
  water_freq_crs <- crs(water_freq_rast)
  perm_water_crs <- st_crs(perm_water_sf)
  
  cat(sprintf("   Original DEM CRS: EPSG:%s\n", 
              st_crs(dem_crs)$epsg %||% "未知"))
  
  # 检查是否为地理坐标系（度单位）
  is_geographic <- st_is_longlat(st_crs(dem_crs))
  
  if (is_geographic) {
    cat("   ⚠ 检测到地理坐标系（单位：度）\n")
    cat("   正在投影到 UTM 坐标系（单位：米）...\n")
    
    # 自动确定合适的 UTM zone
    # 获取数据中心点
    bbox <- st_bbox(st_as_sfc(st_bbox(dem_rast)))
    center_lon <- mean(c(bbox["xmin"], bbox["xmax"]))
    center_lat <- mean(c(bbox["ymin"], bbox["ymax"]))
    
    # 计算 UTM zone
    utm_zone <- floor((center_lon + 180) / 6) + 1
    utm_epsg <- ifelse(center_lat >= 0, 
                       32600 + utm_zone,  # 北半球
                       32700 + utm_zone)  # 南半球
    
    cat(sprintf("   目标坐标系: UTM Zone %d (EPSG:%d)\n", utm_zone, utm_epsg))
    
    # 投影所有数据
    dem_rast <- project(dem_rast, paste0("EPSG:", utm_epsg), method = "bilinear")
    water_freq_rast <- project(water_freq_rast, paste0("EPSG:", utm_epsg), 
                                method = "bilinear")
    perm_water_sf <- st_transform(perm_water_sf, crs = utm_epsg)
    
    dem_crs <- crs(dem_rast)
    cat("   ✓ 投影完成\n")
  } else {
    cat("   CRS 已经是投影坐标系\n")
  }
  
  # 检查 CRS 是否一致
  if (dem_crs != crs(water_freq_rast)) {
    cat("   Reprojecting water frequency to DEM CRS...\n")
    water_freq_rast <- project(water_freq_rast, dem_rast, method = "bilinear")
  }
  
  # 如果永久水域 CRS 不同，重投影
  if (!st_crs(perm_water_sf) == st_crs(dem_crs)) {
    cat("   Reprojecting permanent water to DEM CRS...\n")
    perm_water_sf <- st_transform(perm_water_sf, crs = dem_crs)
  }
  
  # 检查分辨率
  dem_res <- res(dem_rast)
  water_freq_res <- res(water_freq_rast)
  
  cat(sprintf("   DEM resolution: %.2f x %.2f\n", dem_res[1], dem_res[2]))
  cat(sprintf("   Water frequency resolution: %.2f x %.2f\n", 
              water_freq_res[1], water_freq_res[2]))
  
  # 如果分辨率不同，重采样水频率栅格到 DEM 分辨率
  if (!all(abs(dem_res - water_freq_res) < 1e-6)) {
    cat("   Resampling water frequency to DEM resolution...\n")
    water_freq_rast <- resample(water_freq_rast, dem_rast, method = "bilinear")
  }
  
  # ---- 3. 创建永久水域栅格 ----
  cat("\n3. Creating permanent water raster mask...\n")
  
  # 将永久水域矢量栅格化到 DEM 网格
  perm_water_rast <- rasterize(perm_water_sf, dem_rast, field = 1, 
                                 background = 0, touches = FALSE)
  cat("   Permanent water rasterized.\n")
  
  # ---- 4. 确保水频率值在 [0, 1] 范围内 ----
  cat("\n4. Normalizing water frequency to [0, 1]...\n")
  
  # 检查最大值
  water_freq_max <- global(water_freq_rast, "max", na.rm = TRUE)[[1]]
  cat(sprintf("   Water frequency max value: %.2f\n", water_freq_max))
  
  # 如果是百分比形式 [0, 100]，转换为 [0, 1]
  if (water_freq_max > 1.5) {
    cat("   Converting from [0, 100] to [0, 1]...\n")
    water_freq_rast <- water_freq_rast / 100
  }
  
  # 确保在 [0, 1] 范围内
  water_freq_rast <- clamp(water_freq_rast, lower = 0, upper = 1, values = TRUE)
  
  # ---- 5. 创建湖区边界（基于永久水域 + 缓冲区）----
  cat("\n5. Creating lake boundary polygon...\n")
  
  # 使用永久水域作为核心，加上适当缓冲
  lake_boundary <- st_union(perm_water_sf)
  lake_boundary <- st_buffer(lake_boundary, dist = 500)  # 500m 缓冲
  lake_boundary <- st_as_sf(data.frame(id = 1, geom = lake_boundary))
  st_geometry(lake_boundary) <- "geom"
  
  cat("   Lake boundary created with 500m buffer.\n")
  
  # ---- 6. 裁剪栅格到湖区范围 ----
  cat("\n6. Masking rasters to lake boundary...\n")
  
  # 裁剪和 mask
  dem_masked <- crop(dem_rast, lake_boundary)
  dem_masked <- mask(dem_masked, lake_boundary)
  
  water_freq_masked <- crop(water_freq_rast, lake_boundary)
  water_freq_masked <- mask(water_freq_masked, lake_boundary)
  
  perm_water_masked <- crop(perm_water_rast, lake_boundary)
  perm_water_masked <- mask(perm_water_masked, lake_boundary)
  
  cat("   Rasters masked to lake boundary.\n")
  
  # ---- 7. 统计摘要 ----
  cat("\n7. Data summary:\n")
  cat(sprintf("   DEM range: %.2f to %.2f m\n", 
              global(dem_masked, "min", na.rm = TRUE)[[1]],
              global(dem_masked, "max", na.rm = TRUE)[[1]]))
  cat(sprintf("   Water frequency range: %.3f to %.3f\n", 
              global(water_freq_masked, "min", na.rm = TRUE)[[1]],
              global(water_freq_masked, "max", na.rm = TRUE)[[1]]))
  cat(sprintf("   Permanent water pixels: %d\n", 
              sum(values(perm_water_masked) == 1, na.rm = TRUE)))
  
  # 计算像元面积（用于后续 A-E 曲线）
  cell_area <- prod(res(dem_masked))  # m^2
  cat(sprintf("   Cell area: %.2f m²\n", cell_area))
  
  cat("\n✓ Data loading and preprocessing complete!\n\n")
  
  # ---- 8. 返回所有数据 ----
  result_data <- list(
    dem = dem_masked,
    water_freq = water_freq_masked,
    perm_water = perm_water_masked,
    lake_boundary = lake_boundary,
    perm_water_polygon = perm_water_sf,
    cell_area = cell_area,
    lake_name = lake_name
  )
  
  # ---- 9. 保存缓存 ----
  # 注意：由于 terra 对象包含外部指针,不能直接保存到 RData
  # 我们将栅格保存为临时文件,然后保存文件路径
  if (use_cache) {
    cat("Saving preprocessed data to cache...\n")
    
    # 为每个栅格创建临时文件
    dem_cache_file <- file.path(cache_dir, sprintf("%s_dem_cache.tif", lake_name))
    water_freq_cache_file <- file.path(cache_dir, sprintf("%s_water_freq_cache.tif", lake_name))
    perm_water_cache_file <- file.path(cache_dir, sprintf("%s_perm_water_cache.tif", lake_name))
    
    # 保存栅格到文件
    writeRaster(dem_masked, dem_cache_file, overwrite = TRUE)
    writeRaster(water_freq_masked, water_freq_cache_file, overwrite = TRUE)
    writeRaster(perm_water_masked, perm_water_cache_file, overwrite = TRUE)
    
    # 保存非栅格对象和栅格路径
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
    cat(sprintf("✓ Cache saved: %s\n\n", cache_file))
  }
  
  return(result_data)
}


#' 构建观测数据框
#' Build observation data frames for INLA
#' 
#' @param data_list load_and_prep_data() 返回的列表
#' @param N_trials Binomial 试验次数（默认 100）
#' @param use_shore_elev 是否使用岸边高程作为观测（默认 TRUE）
#' @param shore_buffer_m 岸边缓冲区宽度（米），用于定义"岸边"区域
#' @param dam_point_path 大坝点shapefile路径（可选，用作最低点约束）
#' @param dam_point_weight 大坝点重复次数（增加约束权重，默认50）
#' @param use_cache 是否使用缓存（默认 TRUE）
#' @param cache_dir 缓存目录（默认 "cache"）
#' @return 包含 elev_df 和 water_df 的列表
build_observation_data <- function(data_list, 
                                    N_trials = 100,
                                    use_shore_elev = TRUE,
                                    shore_buffer_m = 200,
                                    dam_point_path = NULL,
                                    dam_point_weight = 50,
                                    use_cache = TRUE,
                                    cache_dir = "cache") {
  
  library(terra)
  
  # ---- 检查缓存 ----
  if (use_cache) {
    # 创建缓存目录
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }
    
    # 缓存文件名
    lake_name <- data_list$lake_name
    cache_file <- file.path(cache_dir, sprintf("%s_obs_data.RData", lake_name))
    
    # 如果缓存存在，直接加载
    if (file.exists(cache_file)) {
      cat("====================================\n")
      cat("Loading cached observation data...\n")
      cat("====================================\n")
      cat(sprintf("   Cache file: %s\n", cache_file))
      
      load(cache_file)  # 加载 cached_obs_data
      
      cat("✓ Cached observation data loaded successfully!\n")
      cat(sprintf("   Skipping observation data construction (saved time!)\n\n"))
      
      return(cached_obs_data)
    }
  }
  
  cat("====================================\n")
  cat("Building observation data frames...\n")
  cat("====================================\n")
  
  dem <- data_list$dem
  water_freq <- data_list$water_freq
  perm_water <- data_list$perm_water
  
  # ---- 1. 构建高程观测数据框（Gaussian likelihood）----
  cat("\n1. Building elevation observation data...\n")
  
  if (use_shore_elev) {
    # 使用岸边区域的 DEM 作为观测
    # 改进策略：只使用接近水域的低海拔岸边点，避免高地污染校准
    
    # 基础岸边 mask：DEM 有效 且 非永久水域
    shore_mask <- !is.na(dem) & (is.na(perm_water) | perm_water == 0)
    
    # 提取初步岸边像元
    elev_coords_all <- xyFromCell(dem, which(values(shore_mask)))
    elev_values_all <- values(dem)[values(shore_mask)]
    
    # 去除 NA
    valid_idx <- !is.na(elev_values_all)
    elev_coords_all <- elev_coords_all[valid_idx, ]
    elev_values_all <- elev_values_all[valid_idx]
    
    # ===== 修复：使用统计准则（3σ）去除异常值 =====
    # 理论依据：3σ准则是统计学标准异常值检测方法
    # 保留所有在均值±3标准差范围内的观测
    # 这比之前的"最低10%"更有统计依据
    if (length(elev_values_all) > 10) {
      z_scores <- abs(scale(elev_values_all))
      keep_idx <- z_scores < 3  # 3σ准则
      
      n_removed <- sum(!keep_idx)
      if (n_removed > 0) {
        cat(sprintf("   Removing %d outliers (%.1f%%) using 3σ criterion\n",
                    n_removed, 100 * n_removed / length(elev_values_all)))
      }
      
      elev_coords <- elev_coords_all[keep_idx, ]
      elev_values <- elev_values_all[keep_idx]
    } else {
      # 样本太少，不进行过滤
      elev_coords <- elev_coords_all
      elev_values <- elev_values_all
    }
    
    elev_df <- data.frame(
      x = elev_coords[, 1],
      y = elev_coords[, 2],
      elev = elev_values
    )
    
    cat(sprintf("   Total shore points: %d\n", length(elev_values_all)))
    cat(sprintf("   Shore points after 3σ filtering: %d\n", nrow(elev_df)))
    cat(sprintf("   Elevation range: %.2f to %.2f m\n", 
                min(elev_df$elev), max(elev_df$elev)))
  } else {
    # 如果不使用岸边高程，返回空数据框
    elev_df <- data.frame(x = numeric(0), y = numeric(0), elev = numeric(0))
    cat("   No elevation observations used.\n")
  }
  
  # ---- 1.5. 添加大坝点作为最低点约束（如果提供）----
  if (!is.null(dam_point_path) && file.exists(dam_point_path)) {
    cat("\n   Adding dam point as lowest elevation constraint...\n")
    
    tryCatch({
      library(sf)
      
      # 读取大坝点
      dam_point <- st_read(dam_point_path, quiet = TRUE)
      
      # 转换到与DEM相同的CRS
      dam_crs <- st_crs(dam_point)
      dem_crs <- crs(dem)
      
      if (!st_is_longlat(dam_crs) && st_is_longlat(st_crs(dem_crs))) {
        cat("   Transforming dam point CRS to match DEM...\n")
        dam_point <- st_transform(dam_point, crs = dem_crs)
      }
      
      dam_coords <- st_coordinates(dam_point)
      
      cat(sprintf("   Dam point location: (%.6f, %.6f)\n", 
                  dam_coords[1], dam_coords[2]))
      
      # 从DEM提取大坝位置的高程（使用更健壮的方法）
      dam_elev_raw <- extract(dem, dam_coords)
      
      # 处理不同的返回格式
      if (is.data.frame(dam_elev_raw)) {
        dam_elevation <- dam_elev_raw[1, ncol(dam_elev_raw)]
      } else if (is.matrix(dam_elev_raw)) {
        dam_elevation <- dam_elev_raw[1, ncol(dam_elev_raw)]
      } else {
        dam_elevation <- as.numeric(dam_elev_raw)[1]
      }
      
      # 检查是否成功提取
      if (length(dam_elevation) > 0 && !is.na(dam_elevation)) {
        cat(sprintf("   Dam elevation from DEM: %.2f m\n", dam_elevation))
        cat(sprintf("   Replicating %d times to increase constraint weight\n", 
                    dam_point_weight))
        
        # 创建大坝点的重复观测（增加权重）
        dam_obs <- data.frame(
          x = rep(dam_coords[1], dam_point_weight),
          y = rep(dam_coords[2], dam_point_weight),
          elev = rep(dam_elevation, dam_point_weight)
        )
        
        # 合并到高程观测中
        elev_df <- rbind(elev_df, dam_obs)
        
        cat(sprintf("   ✓ Total elevation observations: %d (including %d dam replicates)\n",
                    nrow(elev_df), dam_point_weight))
      } else {
        cat("   ⚠ Warning: Could not extract valid elevation at dam point location\n")
        cat("   Continuing without dam constraint...\n")
      }
    }, error = function(e) {
      cat(sprintf("   ⚠ Warning: Error processing dam point: %s\n", e$message))
      cat("   Continuing without dam constraint...\n")
    })
  }
  
  # ---- 2. 构建水频率观测数据框（Binomial likelihood）----
  cat("\n2. Building water frequency observation data...\n")
  
  # 获取所有有水频率数据的像元
  water_valid <- !is.na(water_freq)
  water_coords <- xyFromCell(water_freq, which(values(water_valid)))
  water_p <- values(water_freq)[values(water_valid)]
  
  # 去除 NA
  valid_idx <- !is.na(water_p)
  water_coords <- water_coords[valid_idx, ]
  water_p <- water_p[valid_idx]
  
  # 对永久水域像元，强制设置 p = 1
  perm_water_cells <- which(values(perm_water) == 1)
  perm_water_coords <- xyFromCell(perm_water, perm_water_cells)
  
  # 合并：先用水频率数据，然后覆盖永久水域
  water_df <- data.frame(
    x = water_coords[, 1],
    y = water_coords[, 2],
    p = water_p
  )
  
  # 添加永久水域点（p = 1）
  if (length(perm_water_cells) > 0) {
    perm_df <- data.frame(
      x = perm_water_coords[, 1],
      y = perm_water_coords[, 2],
      p = 1.0
    )
    
    # 去除重复（如果水频率已经包含这些点）
    # 简单起见，我们直接覆盖
    # 找到永久水域在 water_df 中的位置
    for (i in 1:nrow(perm_df)) {
      match_idx <- which(water_df$x == perm_df$x[i] & water_df$y == perm_df$y[i])
      if (length(match_idx) > 0) {
        water_df$p[match_idx] <- 1.0
      } else {
        water_df <- rbind(water_df, perm_df[i, ])
      }
    }
  }
  
  # 转换为 Binomial 计数
  water_df$N <- N_trials
  water_df$y_water <- round(water_df$p * water_df$N)
  
  # 确保 y_water 在 [0, N] 范围内
  water_df$y_water <- pmax(0, pmin(water_df$N, water_df$y_water))
  
  cat(sprintf("   Water frequency observations: %d points\n", nrow(water_df)))
  cat(sprintf("   Water occurrence range: %.3f to %.3f\n", 
              min(water_df$p), max(water_df$p)))
  cat(sprintf("   Binomial counts (y_water): %d to %d (N = %d)\n", 
              min(water_df$y_water), max(water_df$y_water), N_trials))
  
  # 统计永久水域点数
  n_permanent <- sum(water_df$y_water == water_df$N)
  cat(sprintf("   Permanent water points (y = N): %d\n", n_permanent))
  
  cat("\n✓ Observation data frames built!\n\n")
  
  result_obs <- list(
    elev_df = elev_df,
    water_df = water_df,
    N_trials = N_trials
  )
  
  # ---- 保存缓存 ----
  if (use_cache) {
    cat("Saving observation data to cache...\n")
    cached_obs_data <- result_obs
    save(cached_obs_data, file = cache_file)
    cat(sprintf("✓ Cache saved: %s\n\n", cache_file))
  }
  
  return(result_obs)
}

