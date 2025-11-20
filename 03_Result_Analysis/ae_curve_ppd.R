################################################################################
# Area-Elevation 曲线计算模块
# Area-Elevation Curve Computation Module
################################################################################

#' 计算 Area-Elevation 曲线
#' Compute Area-Elevation curve from bathymetry
#' 
#' @param bathy_result reconstruct_bathymetry() 返回的列表
#' @param data_list load_and_prep_data() 返回的列表
#' @param elevation_step 高程步长（米）
#' @return 包含 A-E 曲线数据的数据框
compute_ae_curve <- function(bathy_result, 
                               data_list,
                               elevation_step = 0.1) {
  
  library(terra)
  
  cat("====================================\n")
  cat("Computing Area-Elevation curve...\n")
  cat("====================================\n")
  
  bathy_mean <- bathy_result$mean
  cell_area <- data_list$cell_area
  
  # ---- 1. 获取高程范围 ----
  cat("\n1. Determining elevation range...\n")
  
  elev_values <- values(bathy_mean)
  elev_values <- elev_values[!is.na(elev_values)]
  
  elev_min <- min(elev_values)
  elev_max <- max(elev_values)
  
  cat(sprintf("   Elevation range: %.2f to %.2f m\n", elev_min, elev_max))
  cat(sprintf("   Elevation step: %.2f m\n", elevation_step))
  
  # ---- 2. 定义水位高度序列 ----
  elevation_levels <- seq(from = elev_min, to = elev_max, by = elevation_step)
  n_levels <- length(elevation_levels)
  
  cat(sprintf("   Number of elevation levels: %d\n", n_levels))
  
  # ---- 3. 计算每个水位下的淹没面积 ----
  cat("\n2. Computing inundated area for each elevation level...\n")
  cat("   (This may take a moment...)\n")
  
  areas <- numeric(n_levels)
  
  # 对每个水位高度 h，计算 Z <= h 的像元数量
  for (i in 1:n_levels) {
    h <- elevation_levels[i]
    
    # 淹没区域：高程 <= h
    inundated <- elev_values <= h
    n_inundated <- sum(inundated)
    
    # 计算面积（像元数 × 像元面积）
    areas[i] <- n_inundated * cell_area
    
    # 进度提示（每 10% 输出一次）
    if (i %% max(1, floor(n_levels / 10)) == 0) {
      cat(sprintf("   Progress: %d%%\n", round(100 * i / n_levels)))
    }
  }
  
  # ---- 4. 构建 A-E 曲线数据框 ----
  ae_df <- data.frame(
    elevation = elevation_levels,
    area = areas,
    area_km2 = areas / 1e6  # 转换为平方公里
  )
  
  cat("\n3. Area-Elevation curve summary:\n")
  cat(sprintf("   Minimum area (at lowest elevation): %.2f m² (%.4f km²)\n", 
              min(ae_df$area), min(ae_df$area_km2)))
  cat(sprintf("   Maximum area (at highest elevation): %.2f m² (%.4f km²)\n", 
              max(ae_df$area), max(ae_df$area_km2)))
  
  cat("\n✓ Area-Elevation curve computed!\n\n")
  
  return(ae_df)
}


#' 绘制 Area-Elevation 曲线
#' Plot Area-Elevation curve
#' 
#' @param ae_df compute_ae_curve() 返回的数据框
#' @param data_list load_and_prep_data() 返回的列表
#' @param save_path 保存图片的路径（可选）
#' @param true_ae_path 真实A-E曲线CSV路径（用于验证对比，可选）
#' @return ggplot 对象
plot_ae_curve <- function(ae_df, data_list, save_path = NULL, true_ae_path = NULL) {
  
  library(ggplot2)
  
  cat("Plotting Area-Elevation curve...\n")
  
  # 准备预测数据
  ae_df$type <- "Predicted (INLA-SPDE)"
  plot_data <- ae_df
  
  # 如果提供了真实A-E曲线CSV，读取并对比
  if (!is.null(true_ae_path) && file.exists(true_ae_path)) {
    cat("   Loading true A-E curve from validation data...\n")
    
    tryCatch({
      # 读取CSV文件（跳过第2行单位说明）
      ae_true_raw <- read.csv(true_ae_path, skip = 1)
      
      # 提取Elevation (m)和Area (km2)列
      # 根据观察，第4列是Elevation (m)，第6列是Area (km2)
      ae_true <- data.frame(
        elevation = ae_true_raw[, 4],
        area_km2 = ae_true_raw[, 6],
        type = "True (from validation data)"
      )
      
      # 移除无效数据
      ae_true <- ae_true[!is.na(ae_true$elevation) & !is.na(ae_true$area_km2), ]
      ae_true <- ae_true[ae_true$area_km2 > 0, ]  # 只保留面积>0的点
      
      cat(sprintf("   Loaded %d true A-E data points\n", nrow(ae_true)))
      cat(sprintf("   True elevation range: %.2f to %.2f m\n", 
                  min(ae_true$elevation), max(ae_true$elevation)))
      
      # 合并数据
      plot_data <- rbind(
        ae_df[, c("elevation", "area_km2", "type")],
        ae_true[, c("elevation", "area_km2", "type")]
      )
      
      # 计算误差统计（在重叠的高程范围内）
      elev_min_common <- max(min(ae_df$elevation), min(ae_true$elevation))
      elev_max_common <- min(max(ae_df$elevation), max(ae_true$elevation))
      
      if (elev_max_common > elev_min_common) {
        # 插值到相同的高程点
        common_elevations <- seq(elev_min_common, elev_max_common, length.out = 100)
        
        pred_interp <- approx(ae_df$elevation, ae_df$area_km2, 
                              xout = common_elevations, rule = 2)$y
        true_interp <- approx(ae_true$elevation, ae_true$area_km2, 
                              xout = common_elevations, rule = 2)$y
        
        # 计算误差
        mae <- mean(abs(pred_interp - true_interp), na.rm = TRUE)
        rmse <- sqrt(mean((pred_interp - true_interp)^2, na.rm = TRUE))
        # 避免除以0
        valid_idx <- true_interp > 0
        mape <- mean(abs((pred_interp[valid_idx] - true_interp[valid_idx]) / 
                         true_interp[valid_idx]) * 100, na.rm = TRUE)
        
        cat(sprintf("   Validation errors (in common elevation range):\n"))
        cat(sprintf("     MAE  = %.4f km²\n", mae))
        cat(sprintf("     RMSE = %.4f km²\n", rmse))
        cat(sprintf("     MAPE = %.2f%%\n", mape))
      }
      
    }, error = function(e) {
      cat(sprintf("   ⚠ Warning: Could not load true A-E curve: %s\n", e$message))
      cat("   Continuing with predicted curve only...\n")
    })
  }
  
  # 绘制曲线（横轴：高程，纵轴：面积）
  if ("type" %in% names(plot_data) && length(unique(plot_data$type)) > 1) {
    # 有对比数据
    p <- ggplot(plot_data, aes(x = elevation, y = area_km2, color = type, linetype = type)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("True (from validation data)" = "black", 
                                      "Predicted (INLA-SPDE)" = "blue")) +
      scale_linetype_manual(values = c("True (from validation data)" = "solid", 
                                        "Predicted (INLA-SPDE)" = "dashed")) +
      labs(
        title = sprintf("Area-Elevation Curve: %s", data_list$lake_name),
        subtitle = "Predicted vs. True comparison",
        x = "Elevation (m)",
        y = expression("Area (km"^2*")"),
        color = "Source",
        linetype = "Source"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95")
      )
  } else {
    # 只有预测数据
    p <- ggplot(plot_data, aes(x = elevation, y = area_km2)) +
      geom_line(color = "blue", linewidth = 1.2) +
      geom_point(color = "darkblue", size = 0.5, alpha = 0.5) +
      labs(
        title = sprintf("Area-Elevation Curve: %s", data_list$lake_name),
        subtitle = "Relationship between water level and inundated area",
        x = "Elevation (m)",
        y = expression("Area (km"^2*")")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray95")
      )
  }
  
  print(p)
  
  # 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat(sprintf("   A-E curve plot saved to: %s\n", save_path))
  }
  
  return(p)
}


#' 计算并绘制 A-E 曲线的不确定性（基于后验样本，可选）
#' Compute and plot A-E curve uncertainty using posterior samples
#' 
#' @param result INLA 拟合结果
#' @param mesh INLA mesh 对象
#' @param data_list load_and_prep_data() 返回的列表
#' @param bathy_result reconstruct_bathymetry() 返回的列表
#' @param n_samples 后验样本数量
#' @param elevation_step 高程步长
#' @param save_path 保存图片的路径（可选）
#' @return 包含均值和置信区间的 A-E 曲线数据框
compute_ae_curve_with_uncertainty <- function(result,
                                               mesh,
                                               data_list,
                                               bathy_result,
                                               n_samples = 100,
                                               elevation_step = 0.5,
                                               save_path = NULL) {
  
  library(INLA)
  library(terra)
  library(ggplot2)
  
  cat("====================================\n")
  cat("Computing A-E curve with uncertainty...\n")
  cat("====================================\n")
  
  cat(sprintf("\n   Number of posterior samples: %d\n", n_samples))
  cat(sprintf("   Elevation step: %.2f m\n", elevation_step))
  cat("\n   (This will take some time...)\n\n")
  
  # ---- 1. 从后验中抽样 ----
  cat("1. Sampling from posterior...\n")
  
  # 使用 inla.posterior.sample 抽样
  samples <- inla.posterior.sample(n = n_samples, result = result)
  
  cat(sprintf("   Drew %d posterior samples.\n", n_samples))
  
  # ---- 2. 对每个样本计算 A-E 曲线 ----
  cat("\n2. Computing A-E curve for each sample...\n")
  
  # 获取预测位置
  pred_coords <- bathy_result$pred_coords
  A_pred <- inla.spde.make.A(mesh = mesh, loc = pred_coords)
  
  # 高程范围
  elev_min <- min(bathy_result$Z_pred_mean)
  elev_max <- max(bathy_result$Z_pred_mean)
  elevation_levels <- seq(from = elev_min, to = elev_max, by = elevation_step)
  n_levels <- length(elevation_levels)
  
  # 存储每个样本的 A-E 曲线
  ae_samples <- matrix(NA, nrow = n_levels, ncol = n_samples)
  
  cell_area <- data_list$cell_area
  
  for (s in 1:n_samples) {
    # 提取当前样本的空间场
    sample_s <- samples[[s]]
    
    # 找到 field 对应的索引
    field_idx <- grep("field:", rownames(sample_s$latent))
    field_sample <- sample_s$latent[field_idx, 1]
    
    # 投影到预测网格
    Z_sample <- as.vector(A_pred %*% field_sample)
    
    # 计算 A-E 曲线
    for (i in 1:n_levels) {
      h <- elevation_levels[i]
      n_inundated <- sum(Z_sample <= h)
      ae_samples[i, s] <- n_inundated * cell_area
    }
    
    # 进度
    if (s %% max(1, floor(n_samples / 10)) == 0) {
      cat(sprintf("   Sample %d/%d (%.0f%%)\n", s, n_samples, 100 * s / n_samples))
    }
  }
  
  # ---- 3. 计算均值和置信区间 ----
  cat("\n3. Computing summary statistics...\n")
  
  ae_mean <- rowMeans(ae_samples)
  ae_lower <- apply(ae_samples, 1, quantile, probs = 0.025)
  ae_upper <- apply(ae_samples, 1, quantile, probs = 0.975)
  
  ae_df_unc <- data.frame(
    elevation = elevation_levels,
    area = ae_mean,
    area_km2 = ae_mean / 1e6,
    lower = ae_lower / 1e6,
    upper = ae_upper / 1e6
  )
  
  # ---- 4. 绘图（横轴：高程，纵轴：面积）----
  cat("\n4. Plotting A-E curve with uncertainty bands...\n")
  
  p <- ggplot(ae_df_unc, aes(x = elevation, y = area_km2)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "lightblue", alpha = 0.5) +
    geom_line(color = "blue", linewidth = 1.2) +
    labs(
      title = sprintf("Area-Elevation Curve with Uncertainty: %s", 
                      data_list$lake_name),
      subtitle = "Mean curve with 95% credible interval",
      x = "Elevation (m)",
      y = expression("Area (km"^2*")")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(face = "bold")
    )
  
  print(p)
  
  # 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat(sprintf("   A-E curve with uncertainty saved to: %s\n", save_path))
  }
  
  cat("\n✓ A-E curve with uncertainty computed!\n\n")
  
  return(list(
    ae_df = ae_df_unc,
    ae_samples = ae_samples,
    plot = p
  ))
}

