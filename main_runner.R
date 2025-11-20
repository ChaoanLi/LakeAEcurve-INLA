################################################################################
# 主运行脚本：湖泊 Area-Elevation 曲线估计
# Main Runner Script: Lake Area-Elevation Curve Estimation
#
# 项目：使用 INLA + SPDE 从多源栅格数据重建湖底测深图并计算 A-E 曲线
# 
# 作者：Cursor AI + User
# 日期：2025
################################################################################

# ============================
# 0. 清理环境和加载包
# ============================

cat("\n")
cat("================================================================================\n")
cat("  Lake Bathymetry Reconstruction using INLA-SPDE\n")
cat("  湖泊测深图重建与 Area-Elevation 曲线估计\n")
cat("================================================================================\n")
cat("\n")

# 清理工作空间
rm(list = ls())
gc()

# 检查并加载必需的包
required_packages <- c("terra", "sf", "ggplot2", "viridis", "patchwork", "INLA")


cat("检查并加载 R 包...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("   安装缺失的包: %s\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("✓ 所有必需包已加载!\n\n")

# 如果 INLA 未安装，需要从特定源安装
if (!require("INLA", quietly = TRUE)) {
  cat("安装 INLA 包...\n")
  install.packages("INLA", 
                   repos = c(getOption("repos"), 
                             INLA = "https://inla.r-inla-download.org/R/stable"), 
                   dep = TRUE)
  library(INLA)
}

# ============================
# 1. 设置路径和参数
# ============================

cat("================================================================================\n")
cat("1. 配置文件路径和参数\n")
cat("================================================================================\n\n")

# 项目根目录（如果在项目根目录运行，则为 "."）
project_dir <- getwd()
cat(sprintf("项目目录: %s\n\n", project_dir))

# 源代码目录
source(file.path(project_dir, "01_Data_Prep", "data_generation.R"))
source(file.path(project_dir, "01_Data_Prep", "mesh_setup.R"))
source(file.path(project_dir, "02_Model_Implementation", "spde_definition.R"))
source(file.path(project_dir, "02_Model_Implementation", "fit_inlabru_model.R"))
source(file.path(project_dir, "03_Result_Analysis", "reconstruct_map.R"))
source(file.path(project_dir, "03_Result_Analysis", "ae_curve_ppd.R"))

cat("✓ 所有模块已加载!\n\n")

# 数据目录
data_dir <- file.path(project_dir, "01_Data_Prep", "data")

# 输出目录（创建如果不存在）
output_dir <- file.path(project_dir, "outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================
# 2. 选择要分析的湖泊
# ============================

cat("================================================================================\n")
cat("2. 选择要分析的湖泊\n")
cat("================================================================================\n\n")

cat("可用湖泊:\n")
cat("  1. Belton Lake\n")
cat("  2. E.V. Spence Reservoir\n")
cat("  3. Both (两个湖泊)\n\n")

# 用户选择（可以修改这里）
lake_choice <- 1  # 1 = Belton, 2 = E.V. Spence, 3 = Both

# 根据选择设置湖泊列表
if (lake_choice == 1) {
  lakes_to_process <- list(
    list(
      name = "Belton",
      dem_path = file.path(data_dir, "DEM", "Belton_DEM_100m_Buffer.tif"),
      water_freq_path = file.path(data_dir, "occurrence", 
                                   "Belton_GSW_Occurrence_Buffer100m.tif"),
      perm_water_path = file.path(data_dir, "permanent_water", 
                                   "Belton_Min_Permanent_Water.shp"),
      dam_point_path = file.path(data_dir, "DAM_point", "Belton_dam.shp"),
      true_ae_path = "04_Validation/Belton_AVE.csv"
    )
  )
} else if (lake_choice == 2) {
  lakes_to_process <- list(
    list(
      name = "EVSpence",
      dem_path = file.path(data_dir, "DEM", "EVSpence_DEM_200m_Buffer.tif"),
      water_freq_path = file.path(data_dir, "occurrence", 
                                   "EVSpence_GSW_Occurrence_Buffer200m.tif"),
      perm_water_path = file.path(data_dir, "permanent_water", 
                                   "EVSpence_Reservoir_Min_Permanent_Water.shp"),
      dam_point_path = file.path(data_dir, "DAM_point", "EVSpence_dam.shp"),
      true_ae_path = "04_Validation/EVspence_AVE.csv"
    )
  )
} else {
  lakes_to_process <- list(
    list(
      name = "Belton",
      dem_path = file.path(data_dir, "DEM", "Belton_DEM_100m_Buffer.tif"),
      water_freq_path = file.path(data_dir, "occurrence", 
                                   "Belton_GSW_Occurrence_Buffer100m.tif"),
      perm_water_path = file.path(data_dir, "permanent_water", 
                                   "Belton_Min_Permanent_Water.shp"),
      dam_point_path = file.path(data_dir, "DAM_point", "Belton_dam.shp"),
      true_ae_path = "04_Validation/Belton_AVE.csv"
    ),
    list(
      name = "EVSpence",
      dem_path = file.path(data_dir, "DEM", "EVSpence_DEM_200m_Buffer.tif"),
      water_freq_path = file.path(data_dir, "occurrence", 
                                   "EVSpence_GSW_Occurrence_Buffer200m.tif"),
      perm_water_path = file.path(data_dir, "permanent_water", 
                                   "EVSpence_Reservoir_Min_Permanent_Water.shp"),
      dam_point_path = file.path(data_dir, "DAM_point", "EVSpence_dam.shp"),
      true_ae_path = "04_Validation/EVspence_AVE.csv"
    )
  )
}

cat(sprintf("将处理 %d 个湖泊\n\n", length(lakes_to_process)))

# ============================
# 3. 主循环：处理每个湖泊
# ============================

cat("================================================================================\n")
cat("3. 开始处理湖泊\n")
cat("================================================================================\n\n")

# 存储结果
all_results <- list()

for (lake_idx in seq_along(lakes_to_process)) {
  
  lake_info <- lakes_to_process[[lake_idx]]
  lake_name <- lake_info$name
  
  cat("\n")
  cat("################################################################################\n")
  cat(sprintf("# Processing Lake: %s (%d/%d)\n", lake_name, lake_idx, 
              length(lakes_to_process)))
  cat("################################################################################\n")
  cat("\n")
  
  # ---- 3.1 加载和预处理数据 ----
  data_list <- load_and_prep_data(
    dem_path = lake_info$dem_path,
    water_freq_path = lake_info$water_freq_path,
    perm_water_path = lake_info$perm_water_path,
    lake_name = lake_name
  )
  
  # ---- 3.2 构建观测数据 ----
  obs_data <- build_observation_data(
    data_list = data_list,
    N_trials = 100,
    use_shore_elev = TRUE,
    dam_point_path = lake_info$dam_point_path,  # 添加大坝点约束
    dam_point_weight = 50  # 大坝点重复50次以增强约束
  )
  
  # 将 obs_data 添加到 data_list 以便后处理使用
  data_list$obs_data <- obs_data
  
  # ---- 3.3 可视化原始数据 ----
  plot_original_data(
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_original_data.png", lake_name))
  )
  
  # ---- 3.4 构建 SPDE mesh ----
  mesh <- build_mesh(
    data_list = data_list,
    max_edge = c(100, 500),
    cutoff = 50,
    offset = c(100, 500)
  )
  
  # 可视化 mesh
  plot_mesh(
    mesh = mesh,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_mesh.png", lake_name))
  )
  
  # ---- 3.5 定义 SPDE 模型 ----
  spde <- define_spde(
    mesh = mesh,
    prior_range = c(500, 0.5),
    prior_sigma = c(1, 0.01)
  )
  
  # ---- 3.6 构建投影矩阵 ----
  proj_matrices <- build_projection_matrices(
    mesh = mesh,
    obs_data = obs_data
  )
  
  # ---- 3.7 构建 INLA stacks ----
  stack_list <- build_inla_stacks(
    spde = spde,
    obs_data = obs_data,
    proj_matrices = proj_matrices
  )
  
  # ---- 3.8 拟合 INLA 模型 ----
  result <- fit_inla_model(
    stack_list = stack_list,
    spde = spde,
    obs_data = obs_data,
    use_elev = TRUE,
    beta_prior = c(0, 0.01)
  )
  
  # ---- 3.9 提取超参数 ----
  hyperpar <- extract_spde_hyperpar(
    result = result,
    spde = spde
  )
  
  # ---- 3.10 重建湖底 DEM ----
  bathy_result <- reconstruct_bathymetry(
    result = result,
    mesh = mesh,
    data_list = data_list,
    template_raster = data_list$dem
  )
  
  # 可视化重建的 DEM
  plot_bathymetry(
    bathy_result = bathy_result,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_bathymetry.png", lake_name))
  )
  
  # ---- 3.11 计算 Area-Elevation 曲线 ----
  ae_df <- compute_ae_curve(
    bathy_result = bathy_result,
    data_list = data_list,
    elevation_step = 0.1
  )
  
  # 绘制 A-E 曲线（包含真实A-E曲线对比）
  plot_ae_curve(
    ae_df = ae_df,
    data_list = data_list,
    save_path = file.path(output_dir, sprintf("%s_ae_curve.png", lake_name)),
    true_ae_path = lake_info$true_ae_path  # 添加真实A-E曲线用于验证对比
  )
  
  # ---- 3.12（可选）计算带不确定性的 A-E 曲线 ----
  # 注意：这一步会比较耗时，如果数据量大可以跳过或减少样本数
  if (FALSE) {  # 设置为 TRUE 以启用
    ae_unc <- compute_ae_curve_with_uncertainty(
      result = result,
      mesh = mesh,
      data_list = data_list,
      bathy_result = bathy_result,
      n_samples = 50,
      elevation_step = 0.5,
      save_path = file.path(output_dir, 
                            sprintf("%s_ae_curve_uncertainty.png", lake_name))
    )
  }
  
  # ---- 3.13 保存结果 ----
  cat("\n")
  cat("Saving results to RData file...\n")
  
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
  
  cat(sprintf("✓ Results saved to: %s_results.RData\n", lake_name))
  
  # 保存 A-E 曲线为 CSV
  write.csv(ae_df, 
            file = file.path(output_dir, sprintf("%s_ae_curve.csv", lake_name)),
            row.names = FALSE)
  
  cat(sprintf("✓ A-E curve saved to: %s_ae_curve.csv\n", lake_name))
  
  # 添加到总结果列表
  all_results[[lake_name]] <- lake_results
  
  cat("\n")
  cat(sprintf("✓✓✓ Lake %s processing complete! ✓✓✓\n", lake_name))
  cat("\n")
}

# ============================
# 4. 最终总结
# ============================

cat("\n")
cat("================================================================================\n")
cat("4. 最终总结\n")
cat("================================================================================\n\n")

cat(sprintf("成功处理了 %d 个湖泊:\n", length(all_results)))
for (lake_name in names(all_results)) {
  cat(sprintf("  - %s\n", lake_name))
}

cat("\n所有输出文件保存在: %s\n", output_dir)
cat("\n输出文件包括:\n")
cat("  - *_original_data.png: 原始数据可视化\n")
cat("  - *_mesh.png: SPDE mesh 可视化\n")
cat("  - *_bathymetry_mean.png: 重建的湖底高程（后验均值）\n")
cat("  - *_bathymetry_sd.png: 湖底高程的不确定性（后验标准差）\n")
cat("  - *_ae_curve.png: Area-Elevation 曲线\n")
cat("  - *_ae_curve.csv: A-E 曲线数据（CSV 格式）\n")
cat("  - *_results.RData: 完整结果（R 对象）\n")

cat("\n")
cat("================================================================================\n")
cat("  分析完成！\n")
cat("  Analysis Complete!\n")
cat("================================================================================\n")
cat("\n")

# ============================
# 5. 打印一些关键结果
# ============================

cat("\n")
cat("================================================================================\n")
cat("5. 关键结果摘要\n")
cat("================================================================================\n\n")

for (lake_name in names(all_results)) {
  lake_res <- all_results[[lake_name]]
  
  cat(sprintf("\n--- %s ---\n", lake_name))
  
  # 模型参数
  fixed_summary <- lake_res$result$summary.fixed
  
  if ("intercept_water" %in% rownames(fixed_summary)) {
    alpha <- fixed_summary["intercept_water", "mean"]
    cat(sprintf("  Alpha (intercept): %.4f\n", alpha))
  }
  
  if ("neg_field_water" %in% rownames(fixed_summary)) {
    beta <- fixed_summary["neg_field_water", "mean"]
    cat(sprintf("  Beta (elevation coef): %.4f\n", beta))
  }
  
  # 空间场参数
  cat(sprintf("  Spatial range: %.2f m\n", lake_res$hyperpar$range$mean))
  cat(sprintf("  Spatial sigma: %.2f m\n", lake_res$hyperpar$sigma$mean))
  
  # A-E 曲线范围
  ae_df <- lake_res$ae_df
  cat(sprintf("  A-E curve: %.2f to %.2f km² (at elev %.2f to %.2f m)\n",
              min(ae_df$area_km2), max(ae_df$area_km2),
              min(ae_df$elevation), max(ae_df$elevation)))
}

cat("\n")
cat("================================================================================\n")
cat("\n程序运行完毕！\n")
cat("Program finished successfully!\n")
cat("\n")
