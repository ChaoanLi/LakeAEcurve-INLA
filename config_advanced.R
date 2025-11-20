################################################################################
# 高级配置文件
# Advanced Configuration File
#
# 在这个文件中自定义所有模型参数和设置
# Customize all model parameters and settings in this file
################################################################################

# ============================
# 1. 项目路径设置
# ============================

CONFIG <- list()

# 项目根目录（通常是当前工作目录）
CONFIG$project_dir <- getwd()

# 数据目录
CONFIG$data_dir <- file.path(CONFIG$project_dir, "01_Data_Prep", "data")

# 输出目录
CONFIG$output_dir <- file.path(CONFIG$project_dir, "outputs")

# 是否创建输出子目录（按日期时间）
CONFIG$use_timestamp_subdir <- FALSE
if (CONFIG$use_timestamp_subdir) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  CONFIG$output_dir <- file.path(CONFIG$output_dir, timestamp)
}

# ============================
# 2. 湖泊数据配置
# ============================

# 定义要分析的湖泊
# 格式：list(name, dem_path, water_freq_path, perm_water_path)

CONFIG$lakes <- list(
  
  # Belton Lake
  belton = list(
    name = "Belton",
    dem_path = file.path(CONFIG$data_dir, "DEM", "Belton_DEM_100m_Buffer.tif"),
    water_freq_path = file.path(CONFIG$data_dir, "occurrence", 
                                 "Belton_GSW_Occurrence_Buffer100m.tif"),
    perm_water_path = file.path(CONFIG$data_dir, "permanent_water", 
                                 "Belton_Min_Permanent_Water.shp"),
    enabled = TRUE               # 是否分析这个湖泊
  ),
  
  # E.V. Spence Reservoir
  evspence = list(
    name = "EVSpence",
    dem_path = file.path(CONFIG$data_dir, "DEM", "EVSpence_DEM_200m_Buffer.tif"),
    water_freq_path = file.path(CONFIG$data_dir, "occurrence", 
                                 "EVSpence_GSW_Occurrence_Buffer200m.tif"),
    perm_water_path = file.path(CONFIG$data_dir, "permanent_water", 
                                 "EVSpence_Reservoir_Min_Permanent_Water.shp"),
    enabled = TRUE
  )
)

# ============================
# 3. 数据预处理参数
# ============================

CONFIG$data_prep <- list(
  
  # 湖区边界缓冲区距离（米）
  # 用于从永久水域创建分析区域边界
  boundary_buffer = 500,
  
  # Binomial 模型的试验次数
  # 水频率 p 会被转换为 y = round(p * N_trials)
  N_trials = 100,
  
  # 是否使用岸边高程作为观测
  use_shore_elev = TRUE,
  
  # 岸边缓冲区宽度（米）
  # 用于定义"岸边"区域（有 DEM 但非永久水域的区域）
  shore_buffer = 200,
  
  # 是否在预处理时下采样数据（减少计算量）
  downsample = FALSE,
  downsample_factor = 2        # 下采样因子（2 = 每隔一个像元采样）
)

# ============================
# 4. SPDE Mesh 参数
# ============================

CONFIG$mesh <- list(
  
  # 三角形最大边长（米）
  # [内部区域, 外部扩展区域]
  # 较小的值 = 更密集的 mesh = 更精细但计算更慢
  max_edge = c(100, 500),
  
  # 最小节点间距（米）
  # 防止节点过于密集
  cutoff = 50,
  
  # 边界外扩展距离（米）
  # [内部扩展, 外部扩展]
  # 用于减少边界效应
  offset = c(100, 500),
  
  # 是否保存 mesh 可视化图
  save_mesh_plot = TRUE
)

# ============================
# 5. SPDE 模型先验
# ============================

CONFIG$spde_priors <- list(
  
  # Range 参数的 PC prior
  # c(range0, prob): P(range < range0) = prob
  # 
  # 解释：我们认为空间相关距离小于 range0 的概率是 prob
  # 例如：c(500, 0.5) 表示我们认为空间相关距离的中位数约为 500 米
  prior_range = c(500, 0.5),
  
  # Sigma（标准差）参数的 PC prior
  # c(sigma0, prob): P(sigma > sigma0) = prob
  # 
  # 解释：我们认为空间场标准差大于 sigma0 的概率是 prob
  # 例如：c(1, 0.01) 表示我们认为标准差大于 1 米的概率很小（1%）
  prior_sigma = c(1, 0.01)
)

# ============================
# 6. INLA 模型参数
# ============================

CONFIG$inla <- list(
  
  # Beta 系数（高程对水频率的影响）的先验
  # c(mean, precision)
  # Beta ~ N(mean, 1/precision)
  # 
  # 注意：Beta 应该是正的（高程越高 → 水频率越低）
  beta_prior_mean = 0,
  beta_prior_precision = 0.01,  # 弱信息先验
  
  # Alpha（截距）的先验
  alpha_prior_mean = 0,
  alpha_prior_precision = 0.001,
  
  # 高程观测精度的先验（如果使用高程数据）
  elev_prior_precision = 0.001,
  
  # INLA 控制参数
  control_compute = list(
    dic = TRUE,                 # 计算 DIC
    waic = TRUE,                # 计算 WAIC
    config = TRUE               # 保存配置（用于后验样本）
  ),
  
  # INLA 拟合选项
  verbose = TRUE,               # 输出详细信息
  
  # 数值稳定性设置
  control_inla = list(
    strategy = "gaussian",      # 或 "simplified.laplace"
    int.strategy = "eb"         # 或 "grid"
  )
)

# ============================
# 7. 后验预测参数
# ============================

CONFIG$prediction <- list(
  
  # 是否使用模板栅格（NULL = 使用 DEM 作为模板）
  template_raster = NULL,
  
  # 预测网格分辨率（如果重新采样）
  # NULL = 使用原始分辨率
  pred_resolution = NULL
)

# ============================
# 8. A-E 曲线计算参数
# ============================

CONFIG$ae_curve <- list(
  
  # 高程步长（米）
  # 较小的值 = 更精细的曲线但计算更慢
  elevation_step = 0.1,
  
  # 是否计算带不确定性的 A-E 曲线（基于后验样本）
  compute_uncertainty = FALSE,
  
  # 后验样本数量（如果计算不确定性）
  # 注意：这会显著增加计算时间
  n_posterior_samples = 100,
  
  # 不确定性计算的高程步长（可以比主曲线更粗）
  elevation_step_uncertainty = 0.5
)

# ============================
# 9. 可视化参数
# ============================

CONFIG$visualization <- list(
  
  # 图形输出格式
  output_format = "png",        # "png", "pdf", "tiff"
  
  # 图形分辨率（DPI）
  dpi = 300,
  
  # 图形尺寸（英寸）
  fig_width = 10,
  fig_height = 8,
  
  # 配色方案
  color_palette_elevation = "viridis",  # viridis, terrain, cividis
  color_palette_uncertainty = "magma",
  color_palette_water = "plasma",
  
  # 是否保存所有中间图形
  save_all_plots = TRUE,
  
  # 是否在屏幕上显示图形
  show_plots = TRUE
)

# ============================
# 10. 计算性能参数
# ============================

CONFIG$performance <- list(
  
  # 并行计算设置（INLA 支持）
  num_threads = 4,              # 使用的线程数（NULL = 自动）
  
  # 是否在拟合后清理大对象以节省内存
  cleanup_after_fit = FALSE,
  
  # 是否保存完整的 INLA 结果对象
  # （包含所有后验信息，文件可能很大）
  save_full_result = TRUE
)

# ============================
# 11. 输出控制
# ============================

CONFIG$output <- list(
  
  # 是否保存 RData 文件
  save_rdata = TRUE,
  
  # 是否保存 CSV 文件（A-E 曲线）
  save_csv = TRUE,
  
  # 是否保存重建的栅格（GeoTIFF）
  save_raster = TRUE,
  
  # 文件名前缀（NULL = 使用湖泊名称）
  file_prefix = NULL,
  
  # 是否输出详细日志
  verbose = TRUE
)

# ============================
# 12. 模型诊断参数
# ============================

CONFIG$diagnostics <- list(
  
  # 是否进行残差诊断
  check_residuals = FALSE,
  
  # 是否绘制后验边际分布
  plot_marginals = FALSE,
  
  # 是否进行交叉验证
  cross_validation = FALSE,
  cv_folds = 5
)

# ============================
# 辅助函数：打印配置摘要
# ============================

print_config_summary <- function(config = CONFIG) {
  cat("\n")
  cat("================================================================\n")
  cat("  配置摘要 | Configuration Summary\n")
  cat("================================================================\n\n")
  
  cat("1. 项目路径:\n")
  cat(sprintf("   - 项目目录: %s\n", config$project_dir))
  cat(sprintf("   - 输出目录: %s\n", config$output_dir))
  cat("\n")
  
  cat("2. 要分析的湖泊:\n")
  for (lake_id in names(config$lakes)) {
    lake <- config$lakes[[lake_id]]
    if (lake$enabled) {
      cat(sprintf("   - %s ✓\n", lake$name))
    } else {
      cat(sprintf("   - %s (disabled)\n", lake$name))
    }
  }
  cat("\n")
  
  cat("3. Mesh 参数:\n")
  cat(sprintf("   - max_edge: [%.1f, %.1f] m\n", 
              config$mesh$max_edge[1], config$mesh$max_edge[2]))
  cat(sprintf("   - cutoff: %.1f m\n", config$mesh$cutoff))
  cat("\n")
  
  cat("4. SPDE 先验:\n")
  cat(sprintf("   - Range prior: P(range < %.1f) = %.2f\n", 
              config$spde_priors$prior_range[1], 
              config$spde_priors$prior_range[2]))
  cat(sprintf("   - Sigma prior: P(sigma > %.2f) = %.2f\n", 
              config$spde_priors$prior_sigma[1], 
              config$spde_priors$prior_sigma[2]))
  cat("\n")
  
  cat("5. A-E 曲线参数:\n")
  cat(sprintf("   - 高程步长: %.2f m\n", config$ae_curve$elevation_step))
  cat(sprintf("   - 计算不确定性: %s\n", 
              ifelse(config$ae_curve$compute_uncertainty, "是", "否")))
  cat("\n")
  
  cat("================================================================\n\n")
}

# ============================
# 辅助函数：验证配置
# ============================

validate_config <- function(config = CONFIG) {
  
  errors <- character(0)
  
  # 检查数据文件是否存在
  for (lake_id in names(config$lakes)) {
    lake <- config$lakes[[lake_id]]
    if (lake$enabled) {
      if (!file.exists(lake$dem_path)) {
        errors <- c(errors, sprintf("DEM 文件不存在: %s", lake$dem_path))
      }
      if (!file.exists(lake$water_freq_path)) {
        errors <- c(errors, sprintf("水频率文件不存在: %s", lake$water_freq_path))
      }
      if (!file.exists(lake$perm_water_path)) {
        errors <- c(errors, sprintf("永久水域文件不存在: %s", lake$perm_water_path))
      }
    }
  }
  
  # 检查参数合理性
  if (config$mesh$cutoff <= 0) {
    errors <- c(errors, "mesh$cutoff 必须大于 0")
  }
  
  if (config$data_prep$N_trials <= 0) {
    errors <- c(errors, "N_trials 必须大于 0")
  }
  
  if (config$ae_curve$elevation_step <= 0) {
    errors <- c(errors, "elevation_step 必须大于 0")
  }
  
  # 报告错误
  if (length(errors) > 0) {
    cat("\n配置验证失败！\n")
    cat("发现以下错误:\n")
    for (err in errors) {
      cat(sprintf("  - %s\n", err))
    }
    return(FALSE)
  } else {
    cat("\n✓ 配置验证通过！\n\n")
    return(TRUE)
  }
}

# ============================
# 导出配置对象
# ============================

# 在主脚本中使用:
# source("config_advanced.R")
# print_config_summary(CONFIG)
# if (validate_config(CONFIG)) {
#   # 运行分析...
# }

cat("\n高级配置文件已加载\n")
cat("使用 print_config_summary(CONFIG) 查看配置摘要\n")
cat("使用 validate_config(CONFIG) 验证配置\n\n")

