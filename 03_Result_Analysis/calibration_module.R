################################################################################
# Calibration Framework for Bathymetry Reconstruction
#
# Implements four calibration methods:
#   1. Regression Calibration (rlm)
#   2. Endpoint-Constrained Affine Transform
#   3. Hybrid Calibration (weighted blend)
#   4. Piecewise Affine Calibration
#
# Automatically selects method with lowest MAE
#
# Authors: Chaoan Li, Yinuo Zhu
# Course: STAT 647, Texas A&M University
################################################################################

library(MASS)  # for rlm
library(INLA)

#' 计算 AE 曲线的拟合指标
#' @param Z_calibrated 校准后的高程向量
#' @param cell_area_km2 像元面积（km²）
#' @param true_elev 真实 AE 曲线高程向量
#' @param true_area 真实 AE 曲线面积向量
#' @return 拟合指标列表 (MAE, RMSE, R², MAPE, NSE)
compute_ae_metrics <- function(Z_calibrated, cell_area_km2, true_elev, true_area) {
  
  # 去除 NA
  Z_valid <- Z_calibrated[!is.na(Z_calibrated)]
  
  if (length(Z_valid) == 0) {
    return(list(mae = Inf, rmse = Inf, r_squared = 0, mape = Inf, nse = -Inf))
  }
  
  # 按高程排序
  Z_sorted <- sort(Z_valid)
  n_pixels <- length(Z_sorted)
  cumulative_area <- (1:n_pixels) * cell_area_km2
  
  # 创建预测 AE 曲线（处理重复值，取平均）
  pred_df <- data.frame(elev = Z_sorted, area = cumulative_area)
  pred_df <- aggregate(area ~ elev, data = pred_df, FUN = max)  # 保留最大面积
  pred_elev <- pred_df$elev
  pred_area <- pred_df$area
  
  # 同样处理真实曲线的重复值
  true_df <- data.frame(elev = true_elev, area = true_area)
  true_df <- aggregate(area ~ elev, data = true_df, FUN = max)
  true_elev_clean <- true_df$elev
  true_area_clean <- true_df$area
  
  # 找到共同的高程范围
  elev_min_common <- max(min(pred_elev), min(true_elev_clean))
  elev_max_common <- min(max(pred_elev), max(true_elev_clean))
  
  if (elev_max_common <= elev_min_common) {
    return(list(mae = Inf, rmse = Inf, r_squared = 0, mape = Inf, nse = -Inf))
  }
  
  # 在共同范围内插值比较
  common_elevations <- seq(elev_min_common, elev_max_common, length.out = 200)
  
  # 使用 ties = "ordered" 处理可能的重复值
  pred_interp <- approx(pred_elev, pred_area, xout = common_elevations, 
                         rule = 2, ties = "ordered")$y
  true_interp <- approx(true_elev_clean, true_area_clean, xout = common_elevations, 
                         rule = 2, ties = "ordered")$y
  
  # 计算各项指标
  residuals <- pred_interp - true_interp
  
  # MAE: Mean Absolute Error
  mae <- mean(abs(residuals), na.rm = TRUE)
  
  # RMSE: Root Mean Square Error
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  # R²: Coefficient of Determination
  ss_res <- sum(residuals^2, na.rm = TRUE)
  ss_tot <- sum((true_interp - mean(true_interp, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- max(0, 1 - ss_res / ss_tot)
  
  # MAPE: Mean Absolute Percentage Error (避免除以0)
  valid_idx <- true_interp > 0.1  # 避免面积接近0时的不稳定
  if (sum(valid_idx) > 0) {
    mape <- mean(abs(residuals[valid_idx] / true_interp[valid_idx]) * 100, na.rm = TRUE)
  } else {
    mape <- NA
  }
  
  # NSE: Nash-Sutcliffe Efficiency
  nse <- 1 - ss_res / ss_tot
  
  return(list(
    mae = mae,
    rmse = rmse,
    r_squared = r_squared,
    mape = mape,
    nse = nse
  ))
}

#' 向后兼容的 MAE 计算函数
compute_ae_mae <- function(Z_calibrated, cell_area_km2, true_elev, true_area) {
  metrics <- compute_ae_metrics(Z_calibrated, cell_area_km2, true_elev, true_area)
  return(metrics$mae)
}


#' 方法 1：回归校准
#' Regression-based calibration using shore DEM observations
#' 
#' 理论依据：
#'   假设 Z_obs = a + b * f(s) + ε
#'   使用最小二乘或鲁棒回归估计 (a, b)
#'
#' @param Z_pred_raw 原始预测场值
#' @param f_shore 岸边位置的场值
#' @param shore_elev 岸边 DEM 高程
#' @return 列表包含 a, b, r_squared, Z_calibrated
calibration_regression <- function(Z_pred_raw, f_shore, shore_elev) {
  
  # 去除异常值（3σ准则）
  z_score_elev <- abs(scale(shore_elev))
  z_score_field <- abs(scale(f_shore))
  valid_idx <- (z_score_elev < 3) & (z_score_field < 3)
  
  shore_elev_clean <- shore_elev[valid_idx]
  f_shore_clean <- f_shore[valid_idx]
  
  if (length(shore_elev_clean) < 10) {
    warning("Too few shore points for regression calibration")
    return(list(a = NA, b = NA, r_squared = NA, Z_calibrated = Z_pred_raw))
  }
  
  # 鲁棒线性回归
  fit <- tryCatch({
    rlm(shore_elev_clean ~ f_shore_clean, maxit = 100)
  }, error = function(e) {
    lm(shore_elev_clean ~ f_shore_clean)
  })
  
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  
  # 计算 R²
  ss_res <- sum(residuals(fit)^2)
  ss_tot <- sum((shore_elev_clean - mean(shore_elev_clean))^2)
  r_squared <- max(0, 1 - ss_res / ss_tot)
  
  # 应用校准
  Z_calibrated <- a + b * Z_pred_raw
  
  return(list(
    a = a,
    b = b,
    r_squared = r_squared,
    Z_calibrated = Z_calibrated,
    method = "Regression",
    description = "Z = a + b * f, using robust linear regression on shore DEM"
  ))
}


#' 方法 2：端点约束仿射校准
#' Endpoint-constrained affine transformation
#'
#' 理论依据：
#'   强制对齐两个边界条件：
#'   1. min(Z) = true_min_elev (坝点)
#'   2. Z at true_max_area = true_max_elev (水面线)
#'
#' @param Z_pred_raw 原始预测场值
#' @param cell_area_km2 像元面积
#' @param true_min_elev 真实最低高程
#' @param true_max_elev 真实最高高程
#' @param true_max_area 真实最大面积
#' @return 列表包含 a, b, Z_calibrated
calibration_endpoint <- function(Z_pred_raw, cell_area_km2, 
                                   true_min_elev, true_max_elev, true_max_area) {
  
  # 按预测值排序
  Z_sorted <- sort(Z_pred_raw)
  n_pixels <- length(Z_sorted)
  cumulative_area <- (1:n_pixels) * cell_area_km2
  
  # 找到累积面积 = true_max_area 的位置
  idx_cutoff <- which.min(abs(cumulative_area - true_max_area))
  Z_cutoff <- Z_sorted[idx_cutoff]
  Z_min <- min(Z_pred_raw)
  
  Z_range <- Z_cutoff - Z_min
  true_range <- true_max_elev - true_min_elev
  
  if (abs(Z_range) < 1e-10) {
    warning("Z range too small for endpoint calibration")
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  # 计算仿射参数
  b <- true_range / Z_range
  a <- true_min_elev - b * Z_min
  
  # 应用校准
  Z_calibrated <- a + b * Z_pred_raw
  
  return(list(
    a = a,
    b = b,
    Z_cutoff = Z_cutoff,
    idx_cutoff = idx_cutoff,
    Z_calibrated = Z_calibrated,
    method = "Endpoint-Constrained",
    description = "Z = a + b * f, forcing min(Z)=z_min and Z(A_max)=z_max"
  ))
}


#' 方法 3：混合校准
#' Hybrid calibration (regression + endpoint blending)
#'
#' 理论依据：
#'   当回归拟合优度高时（R² 高），信任回归结果
#'   当 R² 低时，更多依赖端点约束
#'   a = w * a_reg + (1-w) * a_end
#'   b = w * b_reg + (1-w) * b_end
#'   w = min(R², 0.8)
#'
#' @param reg_result 回归校准结果
#' @param end_result 端点校准结果
#' @param Z_pred_raw 原始预测场值
#' @return 列表包含 a, b, w, Z_calibrated
calibration_hybrid <- function(reg_result, end_result, Z_pred_raw) {
  
  if (is.na(reg_result$a) || is.na(end_result$a)) {
    warning("Cannot perform hybrid calibration: missing regression or endpoint result")
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  # 权重基于回归 R²
  r_squared <- reg_result$r_squared
  w <- min(max(r_squared, 0), 0.8)  # 权重在 [0, 0.8]
  
  # 加权组合
  a <- w * reg_result$a + (1 - w) * end_result$a
  b <- w * reg_result$b + (1 - w) * end_result$b
  
  # 应用校准
  Z_calibrated <- a + b * Z_pred_raw
  
  return(list(
    a = a,
    b = b,
    w = w,
    r_squared = r_squared,
    Z_calibrated = Z_calibrated,
    method = "Hybrid",
    description = sprintf("Z = a + b * f, with w=%.2f blending (based on R²=%.3f)", 
                          w, r_squared)
  ))
}


#' 方法 4：分段仿射校准
#' Piecewise affine calibration
#'
#' 理论依据：
#'   将高程范围分为多段，每段独立计算仿射参数
#'   可以更好地拟合 AE 曲线的非线性形状
#'   在分段边界强制连续性
#'
#' @param Z_pred_raw 原始预测场值
#' @param cell_area_km2 像元面积
#' @param true_elev 真实高程向量
#' @param true_area 真实面积向量
#' @param n_segments 分段数量（默认 4）
#' @return 列表包含分段参数和 Z_calibrated
calibration_piecewise <- function(Z_pred_raw, cell_area_km2, 
                                    true_elev, true_area, n_segments = 4) {
  
  # 按预测值排序
  Z_sorted <- sort(Z_pred_raw)
  n_pixels <- length(Z_sorted)
  cumulative_area <- (1:n_pixels) * cell_area_km2
  
  # 创建预测 AE 曲线（未校准）
  pred_elev_raw <- Z_sorted
  pred_area <- cumulative_area
  
  # 找到共同的面积范围
  area_min_common <- max(min(pred_area), min(true_area[true_area > 0]))
  area_max_common <- min(max(pred_area), max(true_area))
  
  if (area_max_common <= area_min_common) {
    warning("No overlapping area range for piecewise calibration")
    return(list(params = NULL, Z_calibrated = Z_pred_raw))
  }
  
  # 按面积分位数定义分段边界
  area_breaks <- quantile(c(area_min_common, area_max_common), 
                           probs = seq(0, 1, length.out = n_segments + 1))
  
  # 对每个分段计算校准参数
  segment_params <- list()
  Z_calibrated <- Z_pred_raw  # 初始化
  
  for (i in 1:n_segments) {
    area_lo <- area_breaks[i]
    area_hi <- area_breaks[i + 1]
    
    # 找到该面积范围对应的预测高程
    idx_lo <- which.min(abs(cumulative_area - area_lo))
    idx_hi <- which.min(abs(cumulative_area - area_hi))
    
    pred_elev_lo <- Z_sorted[idx_lo]
    pred_elev_hi <- Z_sorted[idx_hi]
    
    # 找到该面积范围对应的真实高程（处理重复值）
    true_elev_lo <- approx(true_area, true_elev, xout = area_lo, rule = 2, ties = "ordered")$y
    true_elev_hi <- approx(true_area, true_elev, xout = area_hi, rule = 2, ties = "ordered")$y
    
    # 计算该段的仿射参数
    pred_range <- pred_elev_hi - pred_elev_lo
    true_range <- true_elev_hi - true_elev_lo
    
    if (abs(pred_range) < 1e-10) {
      b_i <- 1
      a_i <- true_elev_lo - pred_elev_lo
    } else {
      b_i <- true_range / pred_range
      a_i <- true_elev_lo - b_i * pred_elev_lo
    }
    
    segment_params[[i]] <- list(
      segment = i,
      area_range = c(area_lo, area_hi),
      pred_range = c(pred_elev_lo, pred_elev_hi),
      true_range = c(true_elev_lo, true_elev_hi),
      a = a_i,
      b = b_i
    )
    
    # 应用该段的校准
    segment_idx <- (Z_pred_raw >= pred_elev_lo) & (Z_pred_raw < pred_elev_hi)
    if (i == n_segments) {
      segment_idx <- (Z_pred_raw >= pred_elev_lo) & (Z_pred_raw <= pred_elev_hi + 1e-6)
    }
    
    Z_calibrated[segment_idx] <- a_i + b_i * Z_pred_raw[segment_idx]
  }
  
  # 处理超出范围的值
  below_range <- Z_pred_raw < Z_sorted[which.min(abs(cumulative_area - area_min_common))]
  above_range <- Z_pred_raw > Z_sorted[which.min(abs(cumulative_area - area_max_common))]
  
  if (sum(below_range) > 0) {
    a1 <- segment_params[[1]]$a
    b1 <- segment_params[[1]]$b
    Z_calibrated[below_range] <- a1 + b1 * Z_pred_raw[below_range]
  }
  
  if (sum(above_range) > 0) {
    an <- segment_params[[n_segments]]$a
    bn <- segment_params[[n_segments]]$b
    Z_calibrated[above_range] <- an + bn * Z_pred_raw[above_range]
  }
  
  return(list(
    params = segment_params,
    n_segments = n_segments,
    Z_calibrated = Z_calibrated,
    method = "Piecewise",
    description = sprintf("Piecewise affine with %d segments, fitted to true AE curve", 
                          n_segments)
  ))
}


#' 完整校准框架：运行所有方法并自动选择最佳
#' Run all calibration methods and automatically select the best one
#'
#' @param result INLA 拟合结果
#' @param mesh INLA mesh 对象
#' @param data_list load_and_prep_data() 返回的列表
#' @param Z_pred_raw 投影到预测网格的原始场值
#' @return 包含所有方法结果和最佳选择的列表
run_calibration_framework <- function(result, mesh, data_list, Z_pred_raw) {
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  CALIBRATION FRAMEWORK: Running All Methods\n")
  cat(strrep("=", 70), "\n\n")
  
  # ---- 0. 准备数据 ----
  cell_area_km2 <- data_list$cell_area / 1e6
  
  # 读取真实 AE 曲线
  true_ae_path <- file.path("04_Validation", 
                             sprintf("%s_AVE.csv", data_list$lake_name))
  
  if (!file.exists(true_ae_path)) {
    stop("True AE curve not found: ", true_ae_path)
  }
  
  ae_true_raw <- read.csv(true_ae_path, skip = 1)
  
  # 自动检测列结构（Belton有6列，EVspence有5列）
  n_cols <- ncol(ae_true_raw)
  true_elev <- ae_true_raw[, 4]  # 第4列始终是高程(m)
  true_area <- ae_true_raw[, n_cols]  # 最后一列是面积(km²)
  
  valid_idx <- !is.na(true_elev) & !is.na(true_area)
  true_elev <- true_elev[valid_idx]
  true_area <- true_area[valid_idx]
  
  true_min_elev <- min(true_elev)
  true_max_area <- max(true_area)
  true_max_elev <- true_elev[which.max(true_area)]
  
  cat("True A-E curve:\n")
  cat(sprintf("  Min elevation: %.2f m\n", true_min_elev))
  cat(sprintf("  Max elevation: %.2f m (at %.2f km²)\n", true_max_elev, true_max_area))
  cat(sprintf("  Range: %.2f m\n\n", true_max_elev - true_min_elev))
  
  # 获取岸边观测
  elev_df <- data_list$obs_data$elev_df
  shore_coords <- as.matrix(elev_df[, c("x", "y")])
  shore_elev <- elev_df$elev
  
  # 投影场到岸边位置
  field_mean <- result$summary.random$field$mean
  A_shore <- inla.spde.make.A(mesh = mesh, loc = shore_coords)
  f_shore <- as.vector(A_shore %*% field_mean)
  
  cat(sprintf("Shore observations: %d points\n", length(shore_elev)))
  cat(sprintf("Shore DEM range: %.2f to %.2f m\n\n", min(shore_elev), max(shore_elev)))
  
  # 存储所有结果
  results <- list()
  all_metrics <- list()
  
  # ---- 方法 1：回归校准 ----
  cat(strrep("-", 50), "\n")
  cat("Method 1: REGRESSION CALIBRATION\n")
  cat(strrep("-", 50), "\n")
  
  reg_result <- calibration_regression(Z_pred_raw, f_shore, shore_elev)
  
  cat(sprintf("  a = %.4f, b = %.4f\n", reg_result$a, reg_result$b))
  cat(sprintf("  Regression R² = %.4f\n", reg_result$r_squared))
  
  metrics_reg <- compute_ae_metrics(reg_result$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE = %.4f km², RMSE = %.4f km², R² = %.4f\n\n", 
              metrics_reg$mae, metrics_reg$rmse, metrics_reg$r_squared))
  
  results$regression <- reg_result
  results$regression$metrics <- metrics_reg
  all_metrics$regression <- metrics_reg
  
  # ---- 方法 2：端点约束校准 ----
  cat(strrep("-", 50), "\n")
  cat("Method 2: ENDPOINT-CONSTRAINED CALIBRATION\n")
  cat(strrep("-", 50), "\n")
  
  end_result <- calibration_endpoint(Z_pred_raw, cell_area_km2,
                                      true_min_elev, true_max_elev, true_max_area)
  
  cat(sprintf("  a = %.4f, b = %.4f\n", end_result$a, end_result$b))
  
  metrics_end <- compute_ae_metrics(end_result$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE = %.4f km², RMSE = %.4f km², R² = %.4f\n\n", 
              metrics_end$mae, metrics_end$rmse, metrics_end$r_squared))
  
  results$endpoint <- end_result
  results$endpoint$metrics <- metrics_end
  all_metrics$endpoint <- metrics_end
  
  # ---- 方法 3：混合校准 ----
  cat(strrep("-", 50), "\n")
  cat("Method 3: HYBRID CALIBRATION\n")
  cat(strrep("-", 50), "\n")
  
  hyb_result <- calibration_hybrid(reg_result, end_result, Z_pred_raw)
  
  cat(sprintf("  a = %.4f, b = %.4f\n", hyb_result$a, hyb_result$b))
  cat(sprintf("  Weight w = %.2f (based on R² = %.4f)\n", hyb_result$w, hyb_result$r_squared))
  
  metrics_hyb <- compute_ae_metrics(hyb_result$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE = %.4f km², RMSE = %.4f km², R² = %.4f\n\n", 
              metrics_hyb$mae, metrics_hyb$rmse, metrics_hyb$r_squared))
  
  results$hybrid <- hyb_result
  results$hybrid$metrics <- metrics_hyb
  all_metrics$hybrid <- metrics_hyb
  
  # ---- 方法 4：分段校准 ----
  cat(strrep("-", 50), "\n")
  cat("Method 4: PIECEWISE CALIBRATION (4 segments)\n")
  cat(strrep("-", 50), "\n")
  
  pw_result <- calibration_piecewise(Z_pred_raw, cell_area_km2, true_elev, true_area, 
                                      n_segments = 4)
  
  cat("  Segment parameters:\n")
  for (i in 1:length(pw_result$params)) {
    p <- pw_result$params[[i]]
    cat(sprintf("    Segment %d: a=%.4f, b=%.4f (area %.2f-%.2f km²)\n",
                i, p$a, p$b, p$area_range[1], p$area_range[2]))
  }
  
  metrics_pw <- compute_ae_metrics(pw_result$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE = %.4f km², RMSE = %.4f km², R² = %.4f\n\n", 
              metrics_pw$mae, metrics_pw$rmse, metrics_pw$r_squared))
  
  results$piecewise <- pw_result
  results$piecewise$metrics <- metrics_pw
  all_metrics$piecewise <- metrics_pw
  
  # ---- 模型选择 ----
  cat(strrep("=", 70), "\n")
  cat("  MODEL SELECTION & PERFORMANCE METRICS\n")
  cat(strrep("=", 70), "\n\n")
  
  method_names <- c("Regression", "Endpoint-Constrained", "Hybrid", "Piecewise")
  mae_values <- c(metrics_reg$mae, metrics_end$mae, metrics_hyb$mae, metrics_pw$mae)
  rmse_values <- c(metrics_reg$rmse, metrics_end$rmse, metrics_hyb$rmse, metrics_pw$rmse)
  r2_values <- c(metrics_reg$r_squared, metrics_end$r_squared, 
                  metrics_hyb$r_squared, metrics_pw$r_squared)
  
  cat("Performance Comparison:\n")
  cat(sprintf("  %-20s  %10s  %10s  %10s\n", "Method", "MAE(km²)", "RMSE(km²)", "R²"))
  cat(sprintf("  %s\n", strrep("-", 55)))
  for (i in 1:4) {
    cat(sprintf("  %-20s  %10.4f  %10.4f  %10.4f\n", 
                method_names[i], mae_values[i], rmse_values[i], r2_values[i]))
  }
  
  best_idx <- which.min(mae_values)
  best_method <- method_names[best_idx]
  best_metrics <- list(all_metrics$regression, all_metrics$endpoint, 
                        all_metrics$hybrid, all_metrics$piecewise)[[best_idx]]
  
  cat(sprintf("\n*** BEST METHOD: %s ***\n", best_method))
  cat(sprintf("    MAE  = %.4f km²\n", best_metrics$mae))
  cat(sprintf("    RMSE = %.4f km²\n", best_metrics$rmse))
  cat(sprintf("    R²   = %.4f\n", best_metrics$r_squared))
  cat(sprintf("    MAPE = %.2f%%\n", best_metrics$mape))
  cat(sprintf("    NSE  = %.4f\n\n", best_metrics$nse))
  
  # 选择最佳校准结果
  if (best_idx == 1) {
    best_result <- results$regression
  } else if (best_idx == 2) {
    best_result <- results$endpoint
  } else if (best_idx == 3) {
    best_result <- results$hybrid
  } else {
    best_result <- results$piecewise
  }
  
  # ---- 输出最终结果 ----
  results$best <- list(
    method = best_method,
    method_idx = best_idx,
    mae = best_metrics$mae,
    rmse = best_metrics$rmse,
    r_squared = best_metrics$r_squared,
    mape = best_metrics$mape,
    nse = best_metrics$nse,
    metrics = best_metrics,
    Z_calibrated = best_result$Z_calibrated,
    params = best_result
  )
  
  results$all_metrics <- data.frame(
    method = method_names,
    mae = mae_values,
    rmse = rmse_values,
    r_squared = r2_values
  )
  
  results$all_mae <- data.frame(
    method = method_names,
    mae = mae_values
  )
  
  results$true_ae <- data.frame(
    elevation = true_elev,
    area_km2 = true_area
  )
  
  cat("Calibration framework complete!\n")
  cat(sprintf("Final calibration: %s\n", best_method))
  cat(sprintf("Final MAE: %.4f km²\n\n", best_metrics$mae))
  
  return(results)
}


#' 生成校准比较图
#' Generate comparison plot of all calibration methods
#'
#' @param calib_results run_calibration_framework() 的返回结果
#' @param data_list load_and_prep_data() 返回的列表
#' @param save_path 保存路径
#' @return ggplot 对象
plot_calibration_comparison <- function(calib_results, data_list, save_path = NULL) {
  
  library(ggplot2)
  
  cell_area_km2 <- data_list$cell_area / 1e6
  true_ae <- calib_results$true_ae
  
  # 计算每种方法的 AE 曲线
  compute_ae_from_Z <- function(Z_cal, cell_area_km2) {
    Z_valid <- Z_cal[!is.na(Z_cal)]
    Z_sorted <- sort(Z_valid)
    cum_area <- (1:length(Z_sorted)) * cell_area_km2
    data.frame(elevation = Z_sorted, area_km2 = cum_area)
  }
  
  ae_reg <- compute_ae_from_Z(calib_results$regression$Z_calibrated, cell_area_km2)
  ae_end <- compute_ae_from_Z(calib_results$endpoint$Z_calibrated, cell_area_km2)
  ae_hyb <- compute_ae_from_Z(calib_results$hybrid$Z_calibrated, cell_area_km2)
  ae_pw <- compute_ae_from_Z(calib_results$piecewise$Z_calibrated, cell_area_km2)
  
  # 构建方法标签（使用 metrics 中的 MAE）
  label_true <- "True"
  label_reg <- sprintf("Regression (MAE=%.2f, R²=%.3f)", 
                        calib_results$regression$metrics$mae, 
                        calib_results$regression$metrics$r_squared)
  label_end <- sprintf("Endpoint (MAE=%.2f, R²=%.3f)", 
                        calib_results$endpoint$metrics$mae,
                        calib_results$endpoint$metrics$r_squared)
  label_hyb <- sprintf("Hybrid (MAE=%.2f, R²=%.3f)", 
                        calib_results$hybrid$metrics$mae,
                        calib_results$hybrid$metrics$r_squared)
  label_pw <- sprintf("Piecewise (MAE=%.2f, R²=%.3f)", 
                       calib_results$piecewise$metrics$mae,
                       calib_results$piecewise$metrics$r_squared)
  
  # 合并数据
  plot_data <- rbind(
    data.frame(true_ae, method = label_true),
    data.frame(ae_reg, method = label_reg),
    data.frame(ae_end, method = label_end),
    data.frame(ae_hyb, method = label_hyb),
    data.frame(ae_pw, method = label_pw)
  )
  
  # 构建颜色和线型映射
  color_values <- c("black", "blue", "red", "purple", "darkgreen")
  names(color_values) <- c(label_true, label_reg, label_end, label_hyb, label_pw)
  
  linetype_values <- c("solid", "dashed", "dotted", "dotdash", "longdash")
  names(linetype_values) <- c(label_true, label_reg, label_end, label_hyb, label_pw)
  
  # 绘图
  p <- ggplot(plot_data, aes(x = elevation, y = area_km2, color = method, linetype = method)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    labs(
      title = sprintf("Calibration Method Comparison: %s", data_list$lake_name),
      subtitle = sprintf("Best Method: %s (MAE = %.4f km²)", 
                         calib_results$best$method, calib_results$best$mae),
      x = "Elevation (m)",
      y = expression("Area (km"^2*")"),
      color = "Method",
      linetype = "Method"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 12, height = 8, dpi = 300)
    cat(sprintf("Calibration comparison plot saved: %s\n", save_path))
  }
  
  print(p)
  return(p)
}


#' 生成校准报告
#' Generate calibration report for documentation
#'
#' @param calib_results run_calibration_framework() 的返回结果
#' @return 字符串报告
generate_calibration_report <- function(calib_results) {
  
  report <- paste0(
    "## Calibration Results Summary\n\n",
    "### MAE Comparison\n\n",
    "| Method | MAE (km²) | Status |\n",
    "|--------|-----------|--------|\n"
  )
  
  for (i in 1:nrow(calib_results$all_mae)) {
    status <- ifelse(calib_results$all_mae$method[i] == calib_results$best$method, 
                     "**Selected**", "")
    report <- paste0(report, 
                     sprintf("| %s | %.4f | %s |\n",
                             calib_results$all_mae$method[i],
                             calib_results$all_mae$mae[i],
                             status))
  }
  
  report <- paste0(report, "\n### Selected Method Details\n\n")
  report <- paste0(report, sprintf("**Method**: %s\n\n", calib_results$best$method))
  report <- paste0(report, sprintf("**Final MAE**: %.4f km²\n\n", calib_results$best$mae))
  
  if (calib_results$best$method != "Piecewise") {
    report <- paste0(report, sprintf("**Parameters**: a = %.4f, b = %.4f\n\n",
                                     calib_results$best$params$a,
                                     calib_results$best$params$b))
  } else {
    report <- paste0(report, "**Piecewise Parameters**:\n\n")
    for (p in calib_results$best$params$params) {
      report <- paste0(report, sprintf("- Segment %d: a=%.4f, b=%.4f\n",
                                       p$segment, p$a, p$b))
    }
  }
  
  return(report)
}
