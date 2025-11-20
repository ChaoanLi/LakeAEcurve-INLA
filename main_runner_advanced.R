################################################################################
# 主运行脚本（高级版）：使用配置文件
# Main Runner Script (Advanced): Using Configuration File
#
# 这个版本从 config_advanced.R 读取所有参数
# This version reads all parameters from config_advanced.R
################################################################################

# ============================
# 0. 初始化
# ============================

cat("\n")
cat("================================================================================\n")
cat("  Lake Bathymetry Reconstruction using INLA-SPDE (Advanced)\n")
cat("  湖泊测深图重建与 Area-Elevation 曲线估计（高级版）\n")
cat("================================================================================\n")
cat("\n")

# 清理工作空间
rm(list = ls())
gc()

# ============================
# 1. 加载配置文件
# ============================

cat("1. 加载配置文件...\n")

project_dir <- getwd()
config_file <- file.path(project_dir, "config_advanced.R")

if (!file.exists(config_file)) {
  stop("配置文件不存在: ", config_file)
}

source(config_file)

# 打印配置摘要
print_config_summary(CONFIG)

# 验证配置
if (!validate_config(CONFIG)) {
  stop("配置验证失败，请检查配置文件")
}

# ============================
# 2. 加载必需的包
# ============================

cat("2. 检查并加载 R 包...\n")

required_packages <- c("terra", "sf", "ggplot2", "viridis", "patchwork", "INLA")


for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("   安装缺失的包: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# INLA 特殊处理
if (!require("INLA", quietly = TRUE)) {
  cat("安装 INLA 包...\n")
  install.packages("INLA", 
                   repos = c(getOption("repos"), 
                             INLA = "https://inla.r-inla-download.org/R/stable"), 
                   dep = TRUE)
  library(INLA)
}

cat("✓ 所有必需包已加载!\n\n")

# 设置 INLA 线程数
if (!is.null(CONFIG$performance$num_threads)) {
  inla.setOption(num.threads = CONFIG$performance$num_threads)
  cat(sprintf("INLA 线程数设置为: %d\n", CONFIG$performance$num_threads))
}

# ============================
# 3. 加载自定义函数模块
# ============================

cat("\n3. 加载自定义函数模块...\n")

source(file.path(CONFIG$project_dir, "01_Data_Prep", "data_generation.R"))
source(file.path(CONFIG$project_dir, "01_Data_Prep", "mesh_setup.R"))
source(file.path(CONFIG$project_dir, "02_Model_Implementation", "spde_definition.R"))
source(file.path(CONFIG$project_dir, "02_Model_Implementation", "fit_inlabru_model.R"))
source(file.path(CONFIG$project_dir, "03_Result_Analysis", "reconstruct_map.R"))
source(file.path(CONFIG$project_dir, "03_Result_Analysis", "ae_curve_ppd.R"))

cat("✓ 所有模块已加载!\n\n")

# ============================
# 4. 创建输出目录
# ============================

if (!dir.exists(CONFIG$output_dir)) {
  dir.create(CONFIG$output_dir, recursive = TRUE)
  cat(sprintf("创建输出目录: %s\n", CONFIG$output_dir))
}

# ============================
# 5. 确定要处理的湖泊列表
# ============================

cat("\n5. 确定要处理的湖泊...\n")

lakes_to_process <- list()
for (lake_id in names(CONFIG$lakes)) {
  lake <- CONFIG$lakes[[lake_id]]
  if (lake$enabled) {
    lakes_to_process[[lake_id]] <- lake
    cat(sprintf("   ✓ %s\n", lake$name))
  }
}

if (length(lakes_to_process) == 0) {
  stop("没有启用的湖泊！请在配置文件中设置 enabled = TRUE")
}

cat(sprintf("\n将处理 %d 个湖泊\n\n", length(lakes_to_process)))

# ============================
# 6. 主循环：处理每个湖泊
# ============================

cat("================================================================================\n")
cat("6. 开始处理湖泊\n")
cat("================================================================================\n\n")

all_results <- list()

for (lake_id in names(lakes_to_process)) {
  
  lake_info <- lakes_to_process[[lake_id]]
  lake_name <- lake_info$name
  
  cat("\n")
  cat("################################################################################\n")
  cat(sprintf("# Processing: %s\n", lake_name))
  cat("################################################################################\n")
  cat("\n")
  
  # ---- 6.1 加载和预处理数据 ----
  data_list <- load_and_prep_data(
    dem_path = lake_info$dem_path,
    water_freq_path = lake_info$water_freq_path,
    perm_water_path = lake_info$perm_water_path,
    lake_name = lake_name
  )
  
  # ---- 6.2 构建观测数据 ----
  obs_data <- build_observation_data(
    data_list = data_list,
    N_trials = CONFIG$data_prep$N_trials,
    use_shore_elev = CONFIG$data_prep$use_shore_elev
  )
  
  # ---- 6.3 可视化原始数据 ----
  if (CONFIG$visualization$save_all_plots) {
    save_path <- file.path(CONFIG$output_dir, 
                           sprintf("%s_original_data.%s", 
                                   lake_name, CONFIG$visualization$output_format))
  } else {
    save_path <- NULL
  }
  
  plot_original_data(data_list = data_list, save_path = save_path)
  
  # ---- 6.4 构建 SPDE mesh ----
  mesh <- build_mesh(
    data_list = data_list,
    max_edge = CONFIG$mesh$max_edge,
    cutoff = CONFIG$mesh$cutoff,
    offset = CONFIG$mesh$offset
  )
  
  # 可视化 mesh
  if (CONFIG$mesh$save_mesh_plot) {
    save_path <- file.path(CONFIG$output_dir, 
                           sprintf("%s_mesh.%s", 
                                   lake_name, CONFIG$visualization$output_format))
  } else {
    save_path <- NULL
  }
  
  plot_mesh(mesh = mesh, data_list = data_list, save_path = save_path)
  
  # ---- 6.5 定义 SPDE 模型 ----
  spde <- define_spde(
    mesh = mesh,
    prior_range = CONFIG$spde_priors$prior_range,
    prior_sigma = CONFIG$spde_priors$prior_sigma
  )
  
  # ---- 6.6 构建投影矩阵 ----
  proj_matrices <- build_projection_matrices(
    mesh = mesh,
    obs_data = obs_data
  )
  
  # ---- 6.7 构建 INLA stacks ----
  stack_list <- build_inla_stacks(
    spde = spde,
    obs_data = obs_data,
    proj_matrices = proj_matrices
  )
  
  # ---- 6.8 拟合 INLA 模型 ----
  result <- fit_inla_model(
    stack_list = stack_list,
    spde = spde,
    obs_data = obs_data,
    use_elev = CONFIG$data_prep$use_shore_elev,
    beta_prior = c(CONFIG$inla$beta_prior_mean, 
                   CONFIG$inla$beta_prior_precision)
  )
  
  # ---- 6.9 提取超参数 ----
  hyperpar <- extract_spde_hyperpar(result = result, spde = spde)
  
  # ---- 6.10 重建湖底 DEM ----
  bathy_result <- reconstruct_bathymetry(
    result = result,
    mesh = mesh,
    data_list = data_list,
    template_raster = CONFIG$prediction$template_raster
  )
  
  # 可视化重建的 DEM
  if (CONFIG$visualization$save_all_plots) {
    save_path <- file.path(CONFIG$output_dir, 
                           sprintf("%s_bathymetry.%s", 
                                   lake_name, CONFIG$visualization$output_format))
  } else {
    save_path <- NULL
  }
  
  plot_bathymetry(bathy_result = bathy_result, 
                  data_list = data_list, 
                  save_path = save_path)
  
  # ---- 6.11 计算 Area-Elevation 曲线 ----
  ae_df <- compute_ae_curve(
    bathy_result = bathy_result,
    data_list = data_list,
    elevation_step = CONFIG$ae_curve$elevation_step
  )
  
  # 绘制 A-E 曲线
  if (CONFIG$visualization$save_all_plots) {
    save_path <- file.path(CONFIG$output_dir, 
                           sprintf("%s_ae_curve.%s", 
                                   lake_name, CONFIG$visualization$output_format))
  } else {
    save_path <- NULL
  }
  
  plot_ae_curve(ae_df = ae_df, data_list = data_list, save_path = save_path)
  
  # ---- 6.12（可选）计算带不确定性的 A-E 曲线 ----
  ae_unc <- NULL
  if (CONFIG$ae_curve$compute_uncertainty) {
    save_path <- file.path(CONFIG$output_dir, 
                           sprintf("%s_ae_curve_uncertainty.%s", 
                                   lake_name, CONFIG$visualization$output_format))
    
    ae_unc <- compute_ae_curve_with_uncertainty(
      result = result,
      mesh = mesh,
      data_list = data_list,
      bathy_result = bathy_result,
      n_samples = CONFIG$ae_curve$n_posterior_samples,
      elevation_step = CONFIG$ae_curve$elevation_step_uncertainty,
      save_path = save_path
    )
  }
  
  # ---- 6.13 保存结果 ----
  cat("\n保存结果...\n")
  
  # 确定文件名前缀
  prefix <- ifelse(is.null(CONFIG$output$file_prefix), 
                   lake_name, 
                   sprintf("%s_%s", CONFIG$output$file_prefix, lake_name))
  
  # 保存 RData
  if (CONFIG$output$save_rdata) {
    lake_results <- list(
      lake_name = lake_name,
      data_list = data_list,
      obs_data = obs_data,
      mesh = mesh,
      spde = spde,
      result = if (CONFIG$performance$save_full_result) result else NULL,
      hyperpar = hyperpar,
      bathy_result = bathy_result,
      ae_df = ae_df,
      ae_unc = ae_unc,
      config = CONFIG  # 保存配置以便重现
    )
    
    save(lake_results, 
         file = file.path(CONFIG$output_dir, sprintf("%s_results.RData", prefix)))
    cat(sprintf("✓ RData 已保存: %s_results.RData\n", prefix))
  }
  
  # 保存 CSV
  if (CONFIG$output$save_csv) {
    write.csv(ae_df, 
              file = file.path(CONFIG$output_dir, sprintf("%s_ae_curve.csv", prefix)),
              row.names = FALSE)
    cat(sprintf("✓ CSV 已保存: %s_ae_curve.csv\n", prefix))
  }
  
  # 保存栅格
  if (CONFIG$output$save_raster) {
    writeRaster(bathy_result$mean, 
                file.path(CONFIG$output_dir, sprintf("%s_bathymetry_mean.tif", prefix)),
                overwrite = TRUE)
    writeRaster(bathy_result$sd, 
                file.path(CONFIG$output_dir, sprintf("%s_bathymetry_sd.tif", prefix)),
                overwrite = TRUE)
    cat(sprintf("✓ 栅格已保存: %s_bathymetry_*.tif\n", prefix))
  }
  
  # 清理内存（如果需要）
  if (CONFIG$performance$cleanup_after_fit) {
    gc()
  }
  
  # 添加到总结果列表
  if (CONFIG$output$save_rdata) {
    all_results[[lake_name]] <- lake_results
  }
  
  cat("\n")
  cat(sprintf("✓✓✓ %s 处理完成! ✓✓✓\n", lake_name))
  cat("\n")
}

# ============================
# 7. 最终总结
# ============================

cat("\n")
cat("================================================================================\n")
cat("7. 最终总结\n")
cat("================================================================================\n\n")

cat(sprintf("成功处理了 %d 个湖泊:\n", length(all_results)))
for (lake_name in names(all_results)) {
  cat(sprintf("  - %s\n", lake_name))
}

cat(sprintf("\n所有输出文件保存在: %s\n", CONFIG$output_dir))

# 打印关键结果摘要
cat("\n")
cat("================================================================================\n")
cat("关键结果摘要\n")
cat("================================================================================\n\n")

for (lake_name in names(all_results)) {
  lake_res <- all_results[[lake_name]]
  
  cat(sprintf("\n--- %s ---\n", lake_name))
  
  # 模型参数
  if (!is.null(lake_res$result)) {
    fixed_summary <- lake_res$result$summary.fixed
    
    if ("intercept_water" %in% rownames(fixed_summary)) {
      alpha <- fixed_summary["intercept_water", "mean"]
      alpha_sd <- fixed_summary["intercept_water", "sd"]
      cat(sprintf("  Alpha (截距): %.4f ± %.4f\n", alpha, alpha_sd))
    }
    
    if ("neg_field_water" %in% rownames(fixed_summary)) {
      beta <- fixed_summary["neg_field_water", "mean"]
      beta_sd <- fixed_summary["neg_field_water", "sd"]
      cat(sprintf("  Beta (高程系数): %.4f ± %.4f\n", beta, beta_sd))
      
      if (beta > 0) {
        cat("    ✓ Beta > 0: 高程越高 → 水频率越低（符合预期）\n")
      } else {
        cat("    ⚠ Beta < 0: 不符合预期，请检查数据\n")
      }
    }
  }
  
  # 空间场参数
  cat(sprintf("  空间 Range: %.2f m (%.2f - %.2f)\n", 
              lake_res$hyperpar$range$mean,
              lake_res$hyperpar$range$quant0.025,
              lake_res$hyperpar$range$quant0.975))
  cat(sprintf("  空间 Sigma: %.2f m (%.2f - %.2f)\n", 
              lake_res$hyperpar$sigma$mean,
              lake_res$hyperpar$sigma$quant0.025,
              lake_res$hyperpar$sigma$quant0.975))
  
  # A-E 曲线范围
  ae_df <- lake_res$ae_df
  cat(sprintf("  A-E 曲线范围:\n"))
  cat(sprintf("    高程: %.2f 到 %.2f m\n", 
              min(ae_df$elevation), max(ae_df$elevation)))
  cat(sprintf("    面积: %.4f 到 %.4f km²\n", 
              min(ae_df$area_km2), max(ae_df$area_km2)))
  
  # 湖底高程统计
  cat(sprintf("  重建的湖底高程:\n"))
  cat(sprintf("    均值范围: %.2f 到 %.2f m\n", 
              min(lake_res$bathy_result$Z_pred_mean), 
              max(lake_res$bathy_result$Z_pred_mean)))
  cat(sprintf("    平均不确定性: %.2f m\n", 
              mean(lake_res$bathy_result$Z_pred_sd)))
}

cat("\n")
cat("================================================================================\n")
cat("  分析完成！\n")
cat("  Analysis Complete!\n")
cat("================================================================================\n")
cat("\n")

cat("程序运行完毕！\n")
cat("Program finished successfully!\n\n")

