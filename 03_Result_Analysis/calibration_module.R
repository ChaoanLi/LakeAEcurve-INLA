################################################################################
# Calibration Framework for Bathymetry Reconstruction
#
# Core idea:
#   - DEM region (above water surface): ALL methods use same DEM calibration
#   - Underwater region: Different extrapolation based on validation points
#
# Methods:
#   1. DEM-only: Use DEM's (a,b) to extrapolate down to lake bottom
#   2. DEM+1pt: DEM above + new (a',b') from DEM-bottom to dam point
#   3. DEM+2pt: DEM above + 2-segment below (0% + 50%)
#   4. DEM+4pt: DEM above + 3-segment below (0% + 25% + 50% + 75%)
#
# Authors: Chaoan Li, Yinuo Zhu
# Course: STAT 647, Texas A&M University
################################################################################

library(INLA)

#' 计算 AE 曲线的拟合指标
compute_ae_metrics <- function(Z_calibrated, cell_area_km2, true_elev, true_area) {
  
  Z_valid <- Z_calibrated[!is.na(Z_calibrated)]
  
  if (length(Z_valid) == 0) {
    return(list(mae = Inf, rmse = Inf, r_squared = 0, mape = Inf, nse = -Inf))
  }
  
  Z_sorted <- sort(Z_valid)
  n_pixels <- length(Z_sorted)
  cumulative_area <- (1:n_pixels) * cell_area_km2
  
  pred_df <- data.frame(elev = Z_sorted, area = cumulative_area)
  pred_df <- aggregate(area ~ elev, data = pred_df, FUN = max)
  pred_elev <- pred_df$elev
  pred_area <- pred_df$area
  
  true_df <- data.frame(elev = true_elev, area = true_area)
  true_df <- aggregate(area ~ elev, data = true_df, FUN = max)
  true_elev_clean <- true_df$elev
  true_area_clean <- true_df$area
  
  elev_min_common <- max(min(pred_elev), min(true_elev_clean))
  elev_max_common <- min(max(pred_elev), max(true_elev_clean))
  
  if (elev_max_common <= elev_min_common) {
    return(list(mae = Inf, rmse = Inf, r_squared = 0, mape = Inf, nse = -Inf))
  }
  
  common_elevations <- seq(elev_min_common, elev_max_common, length.out = 200)
  
  pred_interp <- approx(pred_elev, pred_area, xout = common_elevations, 
                         rule = 2, ties = "ordered")$y
  true_interp <- approx(true_elev_clean, true_area_clean, xout = common_elevations, 
                         rule = 2, ties = "ordered")$y
  
  residuals <- pred_interp - true_interp
  
  mae <- mean(abs(residuals), na.rm = TRUE)
  rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
  
  ss_res <- sum(residuals^2, na.rm = TRUE)
  ss_tot <- sum((true_interp - mean(true_interp, na.rm = TRUE))^2, na.rm = TRUE)
  r_squared <- max(0, 1 - ss_res / ss_tot)
  
  valid_idx <- true_interp > 0.1
  if (sum(valid_idx) > 0) {
    mape <- mean(abs(residuals[valid_idx] / true_interp[valid_idx]) * 100, na.rm = TRUE)
  } else {
    mape <- NA
  }
  
  nse <- 1 - ss_res / ss_tot
  
  return(list(mae = mae, rmse = rmse, r_squared = r_squared, mape = mape, nse = nse))
}


#' 基础 DEM 校准：用岸边 DEM 得到 (a, b)
#' 这是所有方法共用的基础
get_dem_calibration <- function(f_shore, shore_elev) {
  
  # 3σ 异常值过滤
  z_score_elev <- abs(scale(shore_elev))
  z_score_field <- abs(scale(f_shore))
  valid_idx <- (z_score_elev < 3) & (z_score_field < 3) & 
               !is.na(z_score_elev) & !is.na(z_score_field)
  
  shore_elev_clean <- shore_elev[valid_idx]
  f_shore_clean <- f_shore[valid_idx]
  
  if (length(shore_elev_clean) < 20) {
    warning("Too few shore points for DEM calibration")
    return(NULL)
  }
  
  # Field 和高程是逆向关系，取负
  probs <- seq(0.05, 0.95, by = 0.05)
  neg_field_quantiles <- quantile(-f_shore_clean, probs = probs)
  dem_quantiles <- quantile(shore_elev_clean, probs = probs)
  
  fit <- lm(dem_quantiles ~ neg_field_quantiles)
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  
  ss_res <- sum(residuals(fit)^2)
  ss_tot <- sum((dem_quantiles - mean(dem_quantiles))^2)
  r_squared <- max(0, 1 - ss_res / ss_tot)
  
  # DEM 覆盖的 field 范围（岸边）
  dem_field_min <- min(-f_shore_clean)  # -field 最小 = 最深处（DEM能看到的最深）
  dem_field_max <- max(-f_shore_clean)  # -field 最大 = 最浅处
  dem_elev_min <- min(shore_elev_clean)
  dem_elev_max <- max(shore_elev_clean)
  
  return(list(
    a = a,
    b = b,
    r_squared = r_squared,
    dem_field_range = c(dem_field_min, dem_field_max),
    dem_elev_range = c(dem_elev_min, dem_elev_max)
  ))
}


#' 方法 1：DEM-only（0 验证点）
#' 用 DEM 的 (a,b) 直接外推到湖底
calibration_dem_only <- function(Z_pred_raw, dem_calib) {
  
  if (is.null(dem_calib)) {
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  # 对所有像素应用 DEM 的 (a, b)
  neg_field <- -Z_pred_raw
  Z_calibrated <- dem_calib$a + dem_calib$b * neg_field
  
  return(list(
    a = dem_calib$a,
    b = dem_calib$b,
    dem_r_squared = dem_calib$r_squared,
    Z_calibrated = Z_calibrated,
    method = "DEM-only",
    n_validation_points = 0,
    description = "DEM calibration extrapolated to lake bottom"
  ))
}


#' 方法 2：DEM + 1个验证点（最低点/坝点）
#' DEM 区域不变，水下区域用 DEM 底部到坝点的线性插值
calibration_dem_plus_1pt <- function(Z_pred_raw, dem_calib, true_min_elev) {
  
  if (is.null(dem_calib)) {
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  neg_field <- -Z_pred_raw
  Z_calibrated <- rep(NA, length(neg_field))
  
  # DEM 区域边界（-field 值）
  dem_field_min <- dem_calib$dem_field_range[1]  # DEM 能看到的最深处的 -field
  dem_elev_at_boundary <- dem_calib$a + dem_calib$b * dem_field_min  # 对应的高程
  
  # 找到场值比 DEM 边界更深的像素（这些是水下区域）
  is_underwater <- neg_field < dem_field_min
  is_dem_region <- neg_field >= dem_field_min
  
  # DEM 区域：用 DEM 的 (a, b)
  Z_calibrated[is_dem_region] <- dem_calib$a + dem_calib$b * neg_field[is_dem_region]
  
  # 水下区域：用 DEM 边界点和坝点（最低点）来计算新的 (a', b')
  # 两点：(dem_field_min, dem_elev_at_boundary) 和 (min(-field), true_min_elev)
  underwater_field_min <- min(neg_field)  # 湖底的 -field 值
  
  if (abs(underwater_field_min - dem_field_min) > 1e-6) {
    b_underwater <- (true_min_elev - dem_elev_at_boundary) / (underwater_field_min - dem_field_min)
    a_underwater <- dem_elev_at_boundary - b_underwater * dem_field_min
    
    Z_calibrated[is_underwater] <- a_underwater + b_underwater * neg_field[is_underwater]
  } else {
    # 如果范围太小，直接用 DEM 参数
    Z_calibrated[is_underwater] <- dem_calib$a + dem_calib$b * neg_field[is_underwater]
  }
  
  return(list(
    a_dem = dem_calib$a,
    b_dem = dem_calib$b,
    a_underwater = if(exists("a_underwater")) a_underwater else dem_calib$a,
    b_underwater = if(exists("b_underwater")) b_underwater else dem_calib$b,
    dem_boundary_elev = dem_elev_at_boundary,
    dem_r_squared = dem_calib$r_squared,
    Z_calibrated = Z_calibrated,
    method = "DEM+1pt",
    n_validation_points = 1,
    description = "DEM above + extrapolate to dam point (0%)"
  ))
}


#' 方法 3：DEM + 2个验证点（0% + 50%）
#' DEM 区域不变，水下区域分2段
calibration_dem_plus_2pt <- function(Z_pred_raw, dem_calib, cell_area_km2, 
                                       true_elev, true_area) {
  
  if (is.null(dem_calib)) {
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  neg_field <- -Z_pred_raw
  Z_calibrated <- rep(NA, length(neg_field))
  
  # DEM 边界
  dem_field_min <- dem_calib$dem_field_range[1]
  dem_elev_at_boundary <- dem_calib$a + dem_calib$b * dem_field_min
  
  is_underwater <- neg_field < dem_field_min
  is_dem_region <- neg_field >= dem_field_min
  
  # DEM 区域
  Z_calibrated[is_dem_region] <- dem_calib$a + dem_calib$b * neg_field[is_dem_region]
  
  # 水下区域分2段：0%-50% 和 50%-DEM边界
  true_max_area <- max(true_area)
  true_min_elev <- min(true_elev)
  
  # 50% 分位点
  area_50 <- 0.50 * true_max_area
  true_elev_50 <- approx(true_area, true_elev, xout = area_50, rule = 2, ties = "ordered")$y
  
  # 找到对应 50% 面积的 -field 值
  neg_field_sorted <- sort(neg_field)
  cumulative_area <- (1:length(neg_field_sorted)) * cell_area_km2
  idx_50 <- which.min(abs(cumulative_area - area_50))
  neg_field_50 <- neg_field_sorted[idx_50]
  
  underwater_field_min <- min(neg_field)  # 最深处
  
  # 段 1: 最深处 (0%) 到 50%
  if (abs(neg_field_50 - underwater_field_min) > 1e-6) {
    b1 <- (true_elev_50 - true_min_elev) / (neg_field_50 - underwater_field_min)
    a1 <- true_min_elev - b1 * underwater_field_min
  } else {
    b1 <- dem_calib$b
    a1 <- true_min_elev - b1 * underwater_field_min
  }
  
  # 段 2: 50% 到 DEM 边界
  if (abs(dem_field_min - neg_field_50) > 1e-6) {
    b2 <- (dem_elev_at_boundary - true_elev_50) / (dem_field_min - neg_field_50)
    a2 <- true_elev_50 - b2 * neg_field_50
  } else {
    b2 <- dem_calib$b
    a2 <- dem_elev_at_boundary - b2 * dem_field_min
  }
  
  # 应用分段校准
  is_seg1 <- is_underwater & (neg_field < neg_field_50)
  is_seg2 <- is_underwater & (neg_field >= neg_field_50)
  
  Z_calibrated[is_seg1] <- a1 + b1 * neg_field[is_seg1]
  Z_calibrated[is_seg2] <- a2 + b2 * neg_field[is_seg2]
  
  return(list(
    a_dem = dem_calib$a,
    b_dem = dem_calib$b,
    segments = list(
      list(range = "0%-50%", a = a1, b = b1),
      list(range = "50%-DEM", a = a2, b = b2)
    ),
    dem_r_squared = dem_calib$r_squared,
    Z_calibrated = Z_calibrated,
    method = "DEM+2pt",
    n_validation_points = 2,
    description = "DEM above + 2 segments below (0% + 50%)"
  ))
}


#' 方法 4：DEM + 4个验证点（0% + 25% + 50% + 75%）
#' DEM 区域不变，水下区域分4段
calibration_dem_plus_4pt <- function(Z_pred_raw, dem_calib, cell_area_km2,
                                       true_elev, true_area) {
  
  if (is.null(dem_calib)) {
    return(list(a = NA, b = NA, Z_calibrated = Z_pred_raw))
  }
  
  neg_field <- -Z_pred_raw
  Z_calibrated <- rep(NA, length(neg_field))
  
  # DEM 边界
  dem_field_min <- dem_calib$dem_field_range[1]
  dem_elev_at_boundary <- dem_calib$a + dem_calib$b * dem_field_min
  
  is_underwater <- neg_field < dem_field_min
  is_dem_region <- neg_field >= dem_field_min
  
  # DEM 区域
  Z_calibrated[is_dem_region] <- dem_calib$a + dem_calib$b * neg_field[is_dem_region]
  
  # 水下区域分4段
  true_max_area <- max(true_area)
  true_min_elev <- min(true_elev)
  
  # 分位点
  probs <- c(0, 0.25, 0.50, 0.75)
  target_areas <- probs * true_max_area
  target_areas[1] <- 0.001 * true_max_area  # 避免 0
  
  true_elevs <- sapply(target_areas, function(a) {
    approx(true_area, true_elev, xout = a, rule = 2, ties = "ordered")$y
  })
  true_elevs[1] <- true_min_elev
  
  # 找到对应的 -field 值
  neg_field_sorted <- sort(neg_field)
  cumulative_area <- (1:length(neg_field_sorted)) * cell_area_km2
  
  target_neg_fields <- sapply(target_areas, function(a) {
    idx <- which.min(abs(cumulative_area - a))
    neg_field_sorted[max(1, idx)]
  })
  
  underwater_field_min <- min(neg_field)
  target_neg_fields[1] <- underwater_field_min
  
  # 添加 DEM 边界作为最后一个点
  all_neg_fields <- c(target_neg_fields, dem_field_min)
  all_elevs <- c(true_elevs, dem_elev_at_boundary)
  
  segments <- list()
  
  for (i in 1:4) {
    nf_lo <- all_neg_fields[i]
    nf_hi <- all_neg_fields[i + 1]
    elev_lo <- all_elevs[i]
    elev_hi <- all_elevs[i + 1]
    
    if (abs(nf_hi - nf_lo) > 1e-6) {
      b_i <- (elev_hi - elev_lo) / (nf_hi - nf_lo)
      a_i <- elev_lo - b_i * nf_lo
    } else {
      b_i <- dem_calib$b
      a_i <- elev_lo - b_i * nf_lo
    }
    
    segments[[i]] <- list(
      range = sprintf("%.0f%%-%.0f%%", probs[i]*100, 
                       ifelse(i < 4, probs[i+1]*100, 100)),
      a = a_i,
      b = b_i
    )
    
    # 应用
    if (i == 1) {
      seg_idx <- is_underwater & (neg_field < nf_hi)
    } else if (i == 4) {
      seg_idx <- is_underwater & (neg_field >= nf_lo)
    } else {
      seg_idx <- is_underwater & (neg_field >= nf_lo) & (neg_field < nf_hi)
    }
    
    Z_calibrated[seg_idx] <- a_i + b_i * neg_field[seg_idx]
  }
  
  return(list(
    a_dem = dem_calib$a,
    b_dem = dem_calib$b,
    segments = segments,
    dem_r_squared = dem_calib$r_squared,
    Z_calibrated = Z_calibrated,
    method = "DEM+4pt",
    n_validation_points = 4,
    description = "DEM above + 4 segments below (0%/25%/50%/75%)"
  ))
}


#' 完整校准框架
run_calibration_framework <- function(result, mesh, data_list, Z_pred_raw) {
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("  CALIBRATION FRAMEWORK\n")
  cat("  DEM region: Same for ALL methods\n")
  cat("  Underwater region: Different extrapolation per method\n")
  cat(strrep("=", 70), "\n\n")
  
  cell_area_km2 <- data_list$cell_area / 1e6
  
  # 读取真实 AE 曲线
  true_ae_path <- file.path("04_Validation", 
                             sprintf("%s_AVE.csv", data_list$lake_name))
  ae_true_raw <- read.csv(true_ae_path, skip = 1)
  n_cols <- ncol(ae_true_raw)
  true_elev <- ae_true_raw[, 4]
  true_area <- ae_true_raw[, n_cols]
  
  valid_idx <- !is.na(true_elev) & !is.na(true_area)
  true_elev <- true_elev[valid_idx]
  true_area <- true_area[valid_idx]
  
  true_min_elev <- min(true_elev)
  true_max_area <- max(true_area)
  
  cat(sprintf("True A-E: elev %.2f-%.2f m, area 0-%.2f km²\n\n", 
              true_min_elev, max(true_elev), true_max_area))
  
  # 获取岸边观测
  elev_df <- data_list$obs_data$elev_df
  shore_coords <- as.matrix(elev_df[, c("x", "y")])
  shore_elev <- elev_df$elev
  
  field_mean <- result$summary.random$field$mean
  A_shore <- inla.spde.make.A(mesh = mesh, loc = shore_coords)
  f_shore <- as.vector(A_shore %*% field_mean)
  
  cat(sprintf("Shore DEM: %d points, elev %.1f-%.1f m\n\n", 
              length(shore_elev), min(shore_elev), max(shore_elev)))
  
  # 基础 DEM 校准（所有方法共用）
  cat("Computing base DEM calibration (shared by all methods)...\n")
  dem_calib <- get_dem_calibration(f_shore, shore_elev)
  cat(sprintf("  DEM: a=%.4f, b=%.4f, R²=%.4f\n", 
              dem_calib$a, dem_calib$b, dem_calib$r_squared))
  cat(sprintf("  DEM elev range: %.2f - %.2f m\n\n", 
              dem_calib$dem_elev_range[1], dem_calib$dem_elev_range[2]))
  
  results <- list()
  all_metrics <- list()
  
  # ---- 方法 1 ----
  cat(strrep("-", 50), "\n")
  cat("Method 1: DEM-only (0 val pts) - Extrapolate with DEM's slope\n")
  cat(strrep("-", 50), "\n")
  
  r1 <- calibration_dem_only(Z_pred_raw, dem_calib)
  m1 <- compute_ae_metrics(r1$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE=%.2f km², R²=%.4f\n\n", m1$mae, m1$r_squared))
  
  results$dem_only <- r1
  results$dem_only$metrics <- m1
  all_metrics$dem_only <- m1
  
  # ---- 方法 2 ----
  cat(strrep("-", 50), "\n")
  cat("Method 2: DEM+1pt (dam point) - DEM above, new slope below\n")
  cat(strrep("-", 50), "\n")
  
  r2 <- calibration_dem_plus_1pt(Z_pred_raw, dem_calib, true_min_elev)
  m2 <- compute_ae_metrics(r2$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE=%.2f km², R²=%.4f\n\n", m2$mae, m2$r_squared))
  
  results$dem_plus_1pt <- r2
  results$dem_plus_1pt$metrics <- m2
  all_metrics$dem_plus_1pt <- m2
  
  # ---- 方法 3 ----
  cat(strrep("-", 50), "\n")
  cat("Method 3: DEM+2pt (0%%+50%%) - DEM above, 2 segments below\n")
  cat(strrep("-", 50), "\n")
  
  r3 <- calibration_dem_plus_2pt(Z_pred_raw, dem_calib, cell_area_km2, true_elev, true_area)
  m3 <- compute_ae_metrics(r3$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE=%.2f km², R²=%.4f\n\n", m3$mae, m3$r_squared))
  
  results$dem_plus_2pt <- r3
  results$dem_plus_2pt$metrics <- m3
  all_metrics$dem_plus_2pt <- m3
  
  # ---- 方法 4 ----
  cat(strrep("-", 50), "\n")
  cat("Method 4: DEM+4pt (0%%/25%%/50%%/75%%) - DEM above, 4 segments below\n")
  cat(strrep("-", 50), "\n")
  
  r4 <- calibration_dem_plus_4pt(Z_pred_raw, dem_calib, cell_area_km2, true_elev, true_area)
  m4 <- compute_ae_metrics(r4$Z_calibrated, cell_area_km2, true_elev, true_area)
  cat(sprintf("  MAE=%.2f km², R²=%.4f\n\n", m4$mae, m4$r_squared))
  
  results$dem_plus_4pt <- r4
  results$dem_plus_4pt$metrics <- m4
  all_metrics$dem_plus_4pt <- m4
  
  # ---- 比较 ----
  cat(strrep("=", 70), "\n")
  cat("  COMPARISON\n")
  cat(strrep("=", 70), "\n\n")
  
  method_names <- c("DEM-only", "DEM+1pt", "DEM+2pt", "DEM+4pt")
  validation_pts <- c(0, 1, 2, 4)
  mae_values <- c(m1$mae, m2$mae, m3$mae, m4$mae)
  r2_values <- c(m1$r_squared, m2$r_squared, m3$r_squared, m4$r_squared)
  
  cat(sprintf("  %-12s  %8s  %10s  %10s\n", "Method", "Val.Pts", "MAE(km²)", "R²"))
  cat(sprintf("  %s\n", strrep("-", 46)))
  for (i in 1:4) {
    cat(sprintf("  %-12s  %8d  %10.2f  %10.4f\n", 
                method_names[i], validation_pts[i], mae_values[i], r2_values[i]))
  }
  
  best_idx <- which.min(mae_values)
  best_method <- method_names[best_idx]
  best_metrics <- list(m1, m2, m3, m4)[[best_idx]]
  best_result <- list(r1, r2, r3, r4)[[best_idx]]
  
  cat(sprintf("\n*** BEST: %s (MAE=%.2f, R²=%.4f) ***\n\n", 
              best_method, best_metrics$mae, best_metrics$r_squared))
  
  results$dem_calib <- dem_calib
  results$best <- list(
    method = best_method,
    method_idx = best_idx,
    mae = best_metrics$mae,
    r_squared = best_metrics$r_squared,
    metrics = best_metrics,
    Z_calibrated = best_result$Z_calibrated
  )
  
  results$all_metrics <- data.frame(
    method = method_names,
    validation_points = validation_pts,
    mae = mae_values,
    r_squared = r2_values
  )
  
  results$true_ae <- data.frame(elevation = true_elev, area_km2 = true_area)
  
  return(results)
}


#' 生成校准比较图
plot_calibration_comparison <- function(calib_results, data_list, save_path = NULL) {
  
  library(ggplot2)
  
  cell_area_km2 <- data_list$cell_area / 1e6
  true_ae <- calib_results$true_ae
  true_max_area <- max(true_ae$area_km2)
  
  compute_ae_full <- function(Z_cal, cell_area_km2) {
    Z_valid <- Z_cal[!is.na(Z_cal)]
    Z_sorted <- sort(Z_valid)
    cum_area <- (1:length(Z_sorted)) * cell_area_km2
    data.frame(elevation = Z_sorted, area_km2 = cum_area)
  }
  
  ae1 <- compute_ae_full(calib_results$dem_only$Z_calibrated, cell_area_km2)
  ae2 <- compute_ae_full(calib_results$dem_plus_1pt$Z_calibrated, cell_area_km2)
  ae3 <- compute_ae_full(calib_results$dem_plus_2pt$Z_calibrated, cell_area_km2)
  ae4 <- compute_ae_full(calib_results$dem_plus_4pt$Z_calibrated, cell_area_km2)
  
  m1 <- calib_results$dem_only$metrics
  m2 <- calib_results$dem_plus_1pt$metrics
  m3 <- calib_results$dem_plus_2pt$metrics
  m4 <- calib_results$dem_plus_4pt$metrics
  
  label_true <- "True (Survey)"
  label_1 <- sprintf("DEM-only: 0pt (MAE=%.1f, R²=%.2f)", m1$mae, m1$r_squared)
  label_2 <- sprintf("DEM+1pt: 0%% (MAE=%.1f, R²=%.2f)", m2$mae, m2$r_squared)
  label_3 <- sprintf("DEM+2pt: 0%%+50%% (MAE=%.1f, R²=%.2f)", m3$mae, m3$r_squared)
  label_4 <- sprintf("DEM+4pt: 0%%/25%%/50%%/75%% (MAE=%.1f, R²=%.2f)", m4$mae, m4$r_squared)
  
  plot_data <- rbind(
    data.frame(true_ae, method = label_true),
    data.frame(ae1, method = label_1),
    data.frame(ae2, method = label_2),
    data.frame(ae3, method = label_3),
    data.frame(ae4, method = label_4)
  )
  
  plot_data$method <- factor(plot_data$method, 
                              levels = c(label_true, label_1, label_2, label_3, label_4))
  
  colors <- c("black", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A")
  names(colors) <- c(label_true, label_1, label_2, label_3, label_4)
  
  linetypes <- c("solid", "dashed", "dotted", "dotdash", "solid")
  names(linetypes) <- c(label_true, label_1, label_2, label_3, label_4)
  
  sizes <- c(1.8, 1.0, 1.0, 1.2, 1.5)
  names(sizes) <- c(label_true, label_1, label_2, label_3, label_4)
  
  # 找到 DEM 区域的边界（所有方法应该在这里重合）
  dem_elev_min <- calib_results$dem_calib$dem_elev_range[1]
  
  p <- ggplot(plot_data, aes(x = elevation, y = area_km2, 
                              color = method, linetype = method, linewidth = method)) +
    geom_line() +
    geom_hline(yintercept = true_max_area, linetype = "dotted", color = "gray50", linewidth = 0.8) +
    geom_vline(xintercept = dem_elev_min, linetype = "dashed", color = "gray60", linewidth = 0.6) +
    annotate("text", x = dem_elev_min + 1, y = max(plot_data$area_km2) * 0.9, 
             label = "DEM boundary", angle = 90, hjust = 1, size = 3.5, color = "gray40") +
    annotate("text", x = max(plot_data$elevation) - 2, y = true_max_area + 2, 
             label = sprintf("True max = %.1f km²", true_max_area), 
             hjust = 1, size = 4, color = "gray40") +
    scale_color_manual(values = colors) +
    scale_linetype_manual(values = linetypes) +
    scale_linewidth_manual(values = sizes) +
    labs(
      title = sprintf("A-E Curve Calibration: %s", data_list$lake_name),
      subtitle = sprintf("DEM region (right of dashed line): ALL methods identical | Best: %s", 
                         calib_results$best$method),
      x = "Elevation (m)",
      y = expression("Cumulative Area (km"^2*")"),
      color = "Method", linetype = "Method", linewidth = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray30"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.key.width = unit(2, "cm"),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1), 
           linewidth = guide_legend(ncol = 1))
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 14, height = 11, dpi = 300)
    cat(sprintf("Plot saved: %s\n", save_path))
  }
  
  print(p)
  return(p)
}


#' 生成校准报告
generate_calibration_report <- function(calib_results) {
  
  report <- paste0(
    "## Calibration Results\n\n",
    "**Key Design**: DEM region is identical for all methods.\n",
    "Differences only in underwater extrapolation.\n\n",
    "### Comparison\n\n",
    "| Method | Val.Pts | MAE (km²) | R² |\n",
    "|--------|---------|-----------|-----|\n"
  )
  
  df <- calib_results$all_metrics
  for (i in 1:nrow(df)) {
    star <- ifelse(df$method[i] == calib_results$best$method, " **" , "")
    report <- paste0(report, sprintf("| %s%s | %d | %.2f | %.4f |\n",
                                     df$method[i], star, df$validation_points[i],
                                     df$mae[i], df$r_squared[i]))
  }
  
  report <- paste0(report, sprintf("\n**Best**: %s\n", calib_results$best$method))
  return(report)
}
