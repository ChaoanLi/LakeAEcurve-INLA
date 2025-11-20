################################################################################
# 示例脚本：单个湖泊分析
# Example Script: Single Lake Analysis
#
# 这个脚本展示如何一步步地运行分析流程
# This script demonstrates step-by-step analysis workflow
################################################################################

# 清理环境
rm(list = ls())
gc()

# 加载必需的包
library(terra)
library(sf)
library(INLA)
library(ggplot2)
library(viridis)

# 加载自定义函数
project_dir <- getwd()
source(file.path(project_dir, "01_Data_Prep", "data_generation.R"))
source(file.path(project_dir, "01_Data_Prep", "mesh_setup.R"))
source(file.path(project_dir, "02_Model_Implementation", "spde_definition.R"))
source(file.path(project_dir, "02_Model_Implementation", "fit_inlabru_model.R"))
source(file.path(project_dir, "03_Result_Analysis", "reconstruct_map.R"))
source(file.path(project_dir, "03_Result_Analysis", "ae_curve_ppd.R"))

cat("\n")
cat("================================================================\n")
cat("  示例：Belton Lake 分析\n")
cat("  Example: Belton Lake Analysis\n")
cat("================================================================\n")
cat("\n")

# ============================
# 步骤 1：设置文件路径
# ============================

cat("步骤 1：设置文件路径...\n")

data_dir <- file.path(project_dir, "01_Data_Prep", "data")

# Belton Lake 数据路径
dem_path <- file.path(data_dir, "DEM", "Belton_DEM_100m_Buffer.tif")
water_freq_path <- file.path(data_dir, "occurrence", 
                              "Belton_GSW_Occurrence_Buffer100m.tif")
perm_water_path <- file.path(data_dir, "permanent_water", 
                              "Belton_Min_Permanent_Water.shp")

# 检查文件是否存在
if (!file.exists(dem_path)) {
  stop("DEM 文件不存在: ", dem_path)
}
if (!file.exists(water_freq_path)) {
  stop("水频率文件不存在: ", water_freq_path)
}
if (!file.exists(perm_water_path)) {
  stop("永久水域文件不存在: ", perm_water_path)
}

cat("✓ 所有数据文件存在\n\n")

# ============================
# 步骤 2：加载和预处理数据
# ============================

cat("步骤 2：加载和预处理数据...\n")

data_list <- load_and_prep_data(
  dem_path = dem_path,
  water_freq_path = water_freq_path,
  perm_water_path = perm_water_path,
  lake_name = "Belton Lake"
)

# 查看数据结构
cat("\n数据列表内容:\n")
cat(sprintf("  - DEM: %s\n", class(data_list$dem)))
cat(sprintf("  - Water frequency: %s\n", class(data_list$water_freq)))
cat(sprintf("  - Permanent water: %s\n", class(data_list$perm_water)))
cat(sprintf("  - Lake boundary: %s\n", class(data_list$lake_boundary)))
cat(sprintf("  - Cell area: %.2f m²\n", data_list$cell_area))

# ============================
# 步骤 3：构建观测数据
# ============================

cat("\n步骤 3：构建观测数据...\n")

obs_data <- build_observation_data(
  data_list = data_list,
  N_trials = 100,              # Binomial 试验次数
  use_shore_elev = TRUE,       # 使用岸边高程
  shore_buffer_m = 200         # 岸边缓冲区
)

# 查看观测数据
cat("\n观测数据摘要:\n")
cat(sprintf("  - 高程观测点数: %d\n", nrow(obs_data$elev_df)))
if (nrow(obs_data$elev_df) > 0) {
  cat(sprintf("  - 高程范围: %.2f 到 %.2f m\n", 
              min(obs_data$elev_df$elev), max(obs_data$elev_df$elev)))
}
cat(sprintf("  - 水频率观测点数: %d\n", nrow(obs_data$water_df)))
cat(sprintf("  - 水频率范围: %.3f 到 %.3f\n", 
            min(obs_data$water_df$p), max(obs_data$water_df$p)))

# ============================
# 步骤 4：可视化原始数据
# ============================

cat("\n步骤 4：可视化原始数据...\n")

plot_original_data(data_list = data_list)

# ============================
# 步骤 5：构建 SPDE mesh
# ============================

cat("\n步骤 5：构建 SPDE mesh...\n")

# 可以调整这些参数来控制 mesh 密度
mesh <- build_mesh(
  data_list = data_list,
  max_edge = c(100, 500),      # 内部/外部最大边长（米）
  cutoff = 50,                 # 最小节点间距（米）
  offset = c(100, 500)         # 边界扩展（米）
)

# 可视化 mesh
plot_mesh(mesh = mesh, data_list = data_list)

# ============================
# 步骤 6：定义 SPDE 模型
# ============================

cat("\n步骤 6：定义 SPDE 模型...\n")

# PC priors：
# - P(spatial range < 500m) = 0.5
# - P(spatial sigma > 1m) = 0.01
spde <- define_spde(
  mesh = mesh,
  prior_range = c(500, 0.5),
  prior_sigma = c(1, 0.01)
)

# ============================
# 步骤 7：构建投影矩阵
# ============================

cat("\n步骤 7：构建投影矩阵...\n")

proj_matrices <- build_projection_matrices(
  mesh = mesh,
  obs_data = obs_data
)

# ============================
# 步骤 8：构建 INLA stacks
# ============================

cat("\n步骤 8：构建 INLA stacks...\n")

stack_list <- build_inla_stacks(
  spde = spde,
  obs_data = obs_data,
  proj_matrices = proj_matrices
)

# ============================
# 步骤 9：拟合 INLA 模型
# ============================

cat("\n步骤 9：拟合 INLA 模型...\n")
cat("（这可能需要几分钟时间...）\n\n")

result <- fit_inla_model(
  stack_list = stack_list,
  spde = spde,
  obs_data = obs_data,
  use_elev = TRUE,
  beta_prior = c(0, 0.01)      # Beta ~ N(0, precision=0.01)
)

# ============================
# 步骤 10：提取超参数
# ============================

cat("\n步骤 10：提取空间场超参数...\n")

hyperpar <- extract_spde_hyperpar(
  result = result,
  spde = spde
)

# ============================
# 步骤 11：重建湖底 DEM
# ============================

cat("\n步骤 11：重建湖底 DEM...\n")

bathy_result <- reconstruct_bathymetry(
  result = result,
  mesh = mesh,
  data_list = data_list,
  template_raster = data_list$dem
)

# 可视化重建的 DEM
plot_bathymetry(
  bathy_result = bathy_result,
  data_list = data_list
)

# ============================
# 步骤 12：计算 A-E 曲线
# ============================

cat("\n步骤 12：计算 Area-Elevation 曲线...\n")

ae_df <- compute_ae_curve(
  bathy_result = bathy_result,
  data_list = data_list,
  elevation_step = 0.1         # 高程步长（米）
)

# 绘制 A-E 曲线
plot_ae_curve(
  ae_df = ae_df,
  data_list = data_list
)

# 查看 A-E 曲线前几行
cat("\nA-E 曲线数据（前 10 行）:\n")
print(head(ae_df, 10))

# ============================
# 步骤 13：保存结果（可选）
# ============================

cat("\n步骤 13：保存结果...\n")

output_dir <- file.path(project_dir, "outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存完整结果
save(
  data_list,
  obs_data,
  mesh,
  spde,
  result,
  hyperpar,
  bathy_result,
  ae_df,
  file = file.path(output_dir, "belton_example_results.RData")
)

# 保存 A-E 曲线
write.csv(
  ae_df,
  file = file.path(output_dir, "belton_ae_curve.csv"),
  row.names = FALSE
)

cat("\n✓ 结果已保存到 outputs/ 文件夹\n")

# ============================
# 步骤 14：输出关键统计量
# ============================

cat("\n")
cat("================================================================\n")
cat("  关键结果摘要\n")
cat("  Key Results Summary\n")
cat("================================================================\n\n")

# 固定效应
fixed_summary <- result$summary.fixed
cat("固定效应 (Fixed Effects):\n")
print(fixed_summary)
cat("\n")

# 空间场参数
cat("空间场参数 (Spatial Field Parameters):\n")
cat(sprintf("  Range (空间相关距离): %.2f m\n", hyperpar$range$mean))
cat(sprintf("  Sigma (标准差): %.2f m\n", hyperpar$sigma$mean))
cat("\n")

# A-E 曲线范围
cat("A-E 曲线 (Area-Elevation Curve):\n")
cat(sprintf("  高程范围: %.2f 到 %.2f m\n", 
            min(ae_df$elevation), max(ae_df$elevation)))
cat(sprintf("  面积范围: %.4f 到 %.4f km²\n", 
            min(ae_df$area_km2), max(ae_df$area_km2)))
cat("\n")

# 湖底高程统计
cat("重建的湖底高程 (Reconstructed Bathymetry):\n")
cat(sprintf("  均值范围: %.2f 到 %.2f m\n", 
            min(bathy_result$Z_pred_mean), max(bathy_result$Z_pred_mean)))
cat(sprintf("  平均不确定性: %.2f m\n", mean(bathy_result$Z_pred_sd)))
cat(sprintf("  最大不确定性: %.2f m\n", max(bathy_result$Z_pred_sd)))
cat("\n")

cat("================================================================\n")
cat("  示例分析完成！\n")
cat("  Example analysis complete!\n")
cat("================================================================\n\n")

cat("提示：您可以修改脚本中的参数来测试不同的模型设置。\n")
cat("Tip: You can modify parameters in the script to test different model settings.\n\n")

