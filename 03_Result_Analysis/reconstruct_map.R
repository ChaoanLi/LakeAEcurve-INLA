################################################################################
# Bathymetry Reconstruction Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' 从 INLA 后验场重建栅格 DEM
#' Reconstruct raster DEM from INLA posterior field
#' 
#' @param result INLA 拟合结果
#' @param mesh INLA mesh 对象
#' @param data_list load_and_prep_data() 返回的列表
#' @param template_raster 模板栅格（用于定义输出栅格的范围和分辨率）
#' @return 包含后验均值和标准差栅格的列表
reconstruct_bathymetry <- function(result, 
                                     mesh, 
                                     data_list,
                                     template_raster = NULL) {
  
  library(INLA)
  library(terra)
  
  cat("====================================\n")
  cat("Reconstructing bathymetry from posterior field...\n")
  cat("====================================\n")
  
  # ---- 1. 使用模板栅格或 DEM 作为参考 ----
  if (is.null(template_raster)) {
    template_raster <- data_list$dem
  }
  
  cat("\n1. Setting up prediction grid...\n")
  cat(sprintf("   Template raster resolution: %.2f x %.2f\n", 
              res(template_raster)[1], res(template_raster)[2]))
  cat(sprintf("   Template raster extent: (%.2f, %.2f) to (%.2f, %.2f)\n",
              ext(template_raster)[1], ext(template_raster)[2],
              ext(template_raster)[3], ext(template_raster)[4]))
  
  # ---- 2. 获取所有像元的中心坐标 ----
  cat("\n2. Extracting prediction locations...\n")
  
  # 获取所有非 NA 像元的坐标（基于湖区 mask）
  template_mask <- !is.na(template_raster)
  pred_cells <- which(values(template_mask))
  pred_coords <- xyFromCell(template_raster, pred_cells)
  
  cat(sprintf("   Number of prediction cells: %d\n", length(pred_cells)))
  
  # ---- 3. 构建预测的 A 矩阵 ----
  cat("\n3. Building projection matrix for prediction...\n")
  
  A_pred <- inla.spde.make.A(
    mesh = mesh,
    loc = pred_coords
  )
  
  cat(sprintf("   A_pred dimensions: %d x %d\n", nrow(A_pred), ncol(A_pred)))
  
  # ---- 4. 提取后验场的均值和标准差 ----
  cat("\n4. Extracting posterior field statistics...\n")
  
  # 获取空间场的后验
  # field 在 result$summary.random$field 中
  field_summary <- result$summary.random$field
  
  # 提取均值和标准差
  field_mean <- field_summary$mean
  field_sd <- field_summary$sd
  
  cat(sprintf("   Field posterior: %d nodes\n", length(field_mean)))
  cat(sprintf("   Field mean range: %.2f to %.2f\n", 
              min(field_mean), max(field_mean)))
  cat(sprintf("   Field SD range: %.2f to %.2f\n", 
              min(field_sd), max(field_sd)))
  
  # ---- 5. 投影到预测网格 ----
  cat("\n5. Projecting to prediction grid...\n")
  
  # 后验均值
  Z_pred_mean <- as.vector(A_pred %*% field_mean)
  
  # 后验标准差（注意：这是近似，因为协方差矩阵很大）
  Z_pred_sd <- as.vector(A_pred %*% field_sd)
  
  cat(sprintf("   Predicted field mean range: %.2f to %.2f\n", 
              min(Z_pred_mean), max(Z_pred_mean)))
  cat(sprintf("   Predicted field SD range: %.2f to %.2f\n", 
              min(Z_pred_sd), max(Z_pred_sd)))
  
  # ---- 6. 校准：使用完整校准框架 ----
  # 运行 4 种校准方法，自动选择 MAE 最低的方案
  cat("\n6. Running Calibration Framework...\n")
  
  # 加载校准模块
  source("03_Result_Analysis/calibration_module.R")
  
  # 检查真实 AE 曲线是否存在
  true_ae_path <- file.path("04_Validation", sprintf("%s_AVE.csv", data_list$lake_name))
  
  if (file.exists(true_ae_path)) {
    # 运行完整校准框架
    calib_results <- run_calibration_framework(
      result = result,
      mesh = mesh,
      data_list = data_list,
      Z_pred_raw = Z_pred_mean
    )
    
    # 使用最佳方法的校准结果
    Z_pred_mean_calibrated <- calib_results$best$Z_calibrated
    calibration_offset <- ifelse(is.null(calib_results$best$params$a), 0, 
                                  calib_results$best$params$a)
    calibration_scale <- ifelse(is.null(calib_results$best$params$b), 1, 
                                 calib_results$best$params$b)
    
    # 保存校准结果（供后续使用）
    calibration_results <- calib_results
    
    cat(sprintf("\n   ✓ Best calibration method: %s\n", calib_results$best$method))
    cat(sprintf("   ✓ Final MAE: %.4f km²\n", calib_results$best$mae))
    cat(sprintf("   ✓ Calibrated range: %.2f to %.2f m\n",
                min(Z_pred_mean_calibrated, na.rm = TRUE),
                max(Z_pred_mean_calibrated, na.rm = TRUE)))
    
  } else {
    cat("   Warning: True A-E curve not found at:", true_ae_path, "\n")
    cat("   Cannot perform calibration without validation data\n")
    cat("   Returning uncalibrated relative bathymetry\n")
    Z_pred_mean_calibrated <- Z_pred_mean
    calibration_offset <- 0
    calibration_scale <- 1
    calibration_results <- NULL
  }
  
  # ---- 7. 转换为栅格 ----
  cat("\n7. Converting to raster format...\n")
  
  # 创建空栅格
  bathy_mean_rast <- rast(template_raster)
  bathy_sd_rast <- rast(template_raster)
  
  # 填充值
  values(bathy_mean_rast) <- NA
  values(bathy_sd_rast) <- NA
  
  bathy_mean_rast[pred_cells] <- Z_pred_mean_calibrated
  bathy_sd_rast[pred_cells] <- Z_pred_sd
  
  # 命名
  names(bathy_mean_rast) <- "elevation_mean"
  names(bathy_sd_rast) <- "elevation_sd"
  
  cat("\n✓ Bathymetry reconstruction complete!\n")
  cat(sprintf("   Calibration applied: %.3f m offset, %.3f scale\n\n", 
              calibration_offset, calibration_scale))
  
  # 返回结果（包含校准信息）
  result_list <- list(
    mean = bathy_mean_rast,
    sd = bathy_sd_rast,
    pred_coords = pred_coords,
    Z_pred_mean = Z_pred_mean,
    Z_pred_sd = Z_pred_sd,
    calibration = list(
      offset = calibration_offset,
      scale = calibration_scale
    )
  )
  
  # 如果有完整校准结果，也包含
  if (exists("calibration_results") && !is.null(calibration_results)) {
    result_list$calibration_results <- calibration_results
  }
  
  return(result_list)
}


#' 可视化重建的湖底 DEM
#' Visualize reconstructed bathymetry
#' 
#' @param bathy_result reconstruct_bathymetry() 返回的列表
#' @param data_list load_and_prep_data() 返回的列表
#' @param save_path 保存图片的路径（可选）
plot_bathymetry <- function(bathy_result, data_list, save_path = NULL) {
  
  library(ggplot2)
  library(terra)
  library(sf)
  library(viridis)
  
  cat("Plotting bathymetry...\n")
  
  # 检查输入
  if (is.null(bathy_result) || is.null(bathy_result$mean) || is.null(bathy_result$sd)) {
    stop("Invalid bathy_result: mean or sd is NULL")
  }
  
  bathy_mean <- bathy_result$mean
  bathy_sd <- bathy_result$sd
  
  # 转换栅格为数据框（用于 ggplot）
  bathy_df <- as.data.frame(bathy_mean, xy = TRUE, na.rm = TRUE)
  if (ncol(bathy_df) < 3) {
    stop("Failed to convert bathy_mean to data.frame")
  }
  colnames(bathy_df) <- c("x", "y", "elevation")
  
  sd_df <- as.data.frame(bathy_sd, xy = TRUE, na.rm = TRUE)
  if (ncol(sd_df) < 3) {
    stop("Failed to convert bathy_sd to data.frame")
  }
  colnames(sd_df) <- c("x", "y", "sd")
  
  cat(sprintf("   Mean elevation data: %d pixels\n", nrow(bathy_df)))
  cat(sprintf("   SD data: %d pixels\n", nrow(sd_df)))
  
  # ---- 图 1：后验均值（高程）----
  cat("   Creating mean elevation plot...\n")
  
  p1 <- ggplot() +
    geom_raster(data = bathy_df, aes(x = x, y = y, fill = elevation)) +
    scale_fill_viridis_c(option = "viridis", name = "Elevation (m)") +
    coord_sf() +
    labs(
      title = sprintf("Reconstructed Bathymetry: %s", data_list$lake_name),
      subtitle = "Posterior Mean Elevation",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  # 添加边界（如果可用）
  if (!is.null(data_list$perm_water_polygon)) {
    p1 <- p1 + geom_sf(data = data_list$perm_water_polygon, 
                       fill = NA, color = "red", linewidth = 0.8, alpha = 0.7)
  }
  if (!is.null(data_list$lake_boundary)) {
    p1 <- p1 + geom_sf(data = data_list$lake_boundary, 
                       fill = NA, color = "black", linewidth = 1, linetype = "dashed")
  }
  
  print(p1)
  cat("   Mean plot created successfully\n")
  
  # ---- 图 2：后验标准差（不确定性）----
  cat("   Creating uncertainty plot...\n")
  
  p2 <- ggplot() +
    geom_raster(data = sd_df, aes(x = x, y = y, fill = sd)) +
    scale_fill_viridis_c(option = "magma", name = "Std Dev (m)") +
    coord_sf() +
    labs(
      title = sprintf("Bathymetry Uncertainty: %s", data_list$lake_name),
      subtitle = "Posterior Standard Deviation",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  # 添加边界（如果可用）
  if (!is.null(data_list$perm_water_polygon)) {
    p2 <- p2 + geom_sf(data = data_list$perm_water_polygon, 
                       fill = NA, color = "cyan", linewidth = 0.8, alpha = 0.7)
  }
  if (!is.null(data_list$lake_boundary)) {
    p2 <- p2 + geom_sf(data = data_list$lake_boundary, 
                       fill = NA, color = "black", linewidth = 1, linetype = "dashed")
  }
  
  print(p2)
  cat("   Uncertainty plot created successfully\n")
  
  # 保存
  if (!is.null(save_path)) {
    # 分别保存两张图
    save_path_mean <- gsub("\\.png$", "_mean.png", save_path)
    save_path_sd <- gsub("\\.png$", "_sd.png", save_path)
    
    ggsave(save_path_mean, p1, width = 10, height = 8, dpi = 300)
    ggsave(save_path_sd, p2, width = 10, height = 8, dpi = 300)
    
    cat(sprintf("   Bathymetry plots saved to:\n"))
    cat(sprintf("     - %s\n", save_path_mean))
    cat(sprintf("     - %s\n", save_path_sd))
  }
  
  return(list(mean_plot = p1, sd_plot = p2))
}


#' 可视化原始数据
#' Visualize original data layers
#' 
#' @param data_list load_and_prep_data() 返回的列表
#' @param save_path 保存图片的路径（可选）
plot_original_data <- function(data_list, save_path = NULL) {
  
  library(ggplot2)
  library(terra)
  library(sf)
  library(viridis)
  library(patchwork)
  
  cat("Plotting original data layers...\n")
  
  dem <- data_list$dem
  water_freq <- data_list$water_freq
  perm_water <- data_list$perm_water
  
  # 转换为数据框
  dem_df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
  colnames(dem_df) <- c("x", "y", "elevation")
  
  water_df <- as.data.frame(water_freq, xy = TRUE, na.rm = TRUE)
  colnames(water_df) <- c("x", "y", "frequency")
  
  perm_df <- as.data.frame(perm_water, xy = TRUE, na.rm = TRUE)
  colnames(perm_df) <- c("x", "y", "permanent")
  perm_df$permanent <- factor(perm_df$permanent, levels = c(0, 1), 
                               labels = c("No", "Yes"))
  
  # ---- 图 1：DEM ----
  p1 <- ggplot() +
    geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
    scale_fill_viridis_c(option = "viridis", name = "Elevation (m)") +
    geom_sf(data = data_list$lake_boundary, 
            fill = NA, color = "black", linewidth = 0.8) +
    coord_sf() +
    labs(title = "Shoreline DEM") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # ---- 图 2：水出现频率 ----
  p2 <- ggplot() +
    geom_raster(data = water_df, aes(x = x, y = y, fill = frequency)) +
    scale_fill_viridis_c(option = "plasma", name = "Frequency") +
    geom_sf(data = data_list$lake_boundary, 
            fill = NA, color = "black", linewidth = 0.8) +
    coord_sf() +
    labs(title = "Water Occurrence Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # ---- 图 3：永久水域 ----
  p3 <- ggplot() +
    geom_raster(data = perm_df, aes(x = x, y = y, fill = permanent)) +
    scale_fill_manual(values = c("No" = "gray90", "Yes" = "blue"), 
                      name = "Permanent\nWater") +
    geom_sf(data = data_list$perm_water_polygon, 
            fill = NA, color = "darkblue", linewidth = 1) +
    geom_sf(data = data_list$lake_boundary, 
            fill = NA, color = "black", linewidth = 0.8) +
    coord_sf() +
    labs(title = "Permanent Water Mask") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 组合图
  p_combined <- (p1 | p2 | p3) +
    plot_annotation(
      title = sprintf("Original Data Layers: %s", data_list$lake_name),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  
  print(p_combined)
  
  # 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p_combined, width = 18, height = 6, dpi = 300)
    cat(sprintf("   Original data plot saved to: %s\n", save_path))
  }
  
  return(p_combined)
}

