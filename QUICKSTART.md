# 快速入门指南 | Quick Start Guide

## 湖泊测深图重建 - 5 分钟上手

---

## 方法 1：最简单（推荐初学者）

### 步骤 1：准备数据

确保您的数据文件位于正确位置：

```
01_Data_Prep/data/
├── DEM/
│   ├── Belton_DEM_100m_Buffer.tif
│   └── EVSpence_DEM_200m_Buffer.tif
├── occurrence/
│   ├── Belton_GSW_Occurrence_Buffer100m.tif
│   └── EVSpence_GSW_Occurrence_Buffer200m.tif
└── permanent_water/
    ├── Belton_Min_Permanent_Water.shp
    └── EVSpence_Reservoir_Min_Permanent_Water.shp
```

### 步骤 2：安装依赖

在 R 中运行：

```r
# 安装基础包
install.packages(c("terra", "sf", "ggplot2", "viridis", "patchwork"))

# 安装 INLA
install.packages("INLA", 
                 repos = c(getOption("repos"), 
                           INLA = "https://inla.r-inla-download.org/R/stable"), 
                 dep = TRUE)
```

### 步骤 3：运行主脚本

```r
# 设置工作目录到项目根目录
setwd("path/to/LakeAEcurve-INLA")

# 运行主脚本
source("main_runner.R")
```

**就这样！** 脚本会自动：
- 加载数据
- 构建 SPDE mesh
- 拟合 INLA 模型
- 重建湖底 DEM
- 计算 A-E 曲线
- 保存所有结果到 `outputs/` 文件夹

---

## 方法 2：单个湖泊示例

如果您只想测试一个湖泊：

```r
setwd("path/to/LakeAEcurve-INLA")
source("example_single_lake.R")
```

这个脚本会：
- 逐步展示每个分析步骤
- 提供详细的进度信息
- 在屏幕上显示所有图形
- 适合学习和调试

---

## 方法 3：自定义参数（高级用户）

### 步骤 1：修改配置文件

打开 `config_advanced.R`，根据需要修改参数：

```r
# 例如：调整 mesh 密度
CONFIG$mesh <- list(
  max_edge = c(50, 300),    # 更密集的 mesh
  cutoff = 25,
  offset = c(100, 500)
)

# 调整 A-E 曲线精度
CONFIG$ae_curve <- list(
  elevation_step = 0.05,     # 更精细的步长
  compute_uncertainty = TRUE, # 启用不确定性分析
  n_posterior_samples = 50
)
```

### 步骤 2：运行高级脚本

```r
source("main_runner_advanced.R")
```

---

## 常见问题 | FAQ

### Q1: 脚本运行需要多长时间？

**A:** 取决于数据大小和计算机性能：
- Belton Lake: 约 5-15 分钟
- E.V. Spence: 约 10-30 分钟
- 带不确定性分析: 额外 10-30 分钟

### Q2: 如何只分析其中一个湖泊？

**A:** 在 `main_runner.R` 中修改：

```r
# 第 2 部分：选择要分析的湖泊
lake_choice <- 1  # 1 = 仅 Belton, 2 = 仅 E.V. Spence, 3 = 两者都分析
```

### Q3: 内存不足怎么办？

**A:** 减小 mesh 密度：

```r
# 在 main_runner.R 的第 3.4 部分
mesh <- build_mesh(
  data_list = data_list,
  max_edge = c(200, 800),    # 增大这些值
  cutoff = 100               # 增大这个值
)
```

### Q4: 如何查看模型拟合质量？

**A:** 查看输出的模型摘要：

```r
# 加载保存的结果
load("outputs/Belton_results.RData")

# 查看固定效应
print(lake_results$result$summary.fixed)

# 查看超参数
print(lake_results$result$summary.hyperpar)

# 查看 DIC/WAIC
print(lake_results$result$dic$dic)
print(lake_results$result$waic$waic)
```

### Q5: 如何修改先验？

**A:** 在 `main_runner.R` 中找到第 3.5 部分：

```r
spde <- define_spde(
  mesh = mesh,
  prior_range = c(500, 0.5),   # 修改这里：P(range < 500m) = 0.5
  prior_sigma = c(1, 0.01)      # 修改这里：P(sigma > 1m) = 0.01
)
```

### Q6: 输出文件在哪里？

**A:** 所有输出保存在 `outputs/` 文件夹：

- `*_original_data.png`: 原始数据可视化
- `*_mesh.png`: SPDE mesh
- `*_bathymetry_mean.png`: 重建的湖底高程
- `*_bathymetry_sd.png`: 不确定性图
- `*_ae_curve.png`: A-E 曲线
- `*_ae_curve.csv`: A-E 曲线数据
- `*_results.RData`: 完整 R 对象

### Q7: 如何解释结果？

**A:** 重点关注：

1. **Intercept**（在模型摘要中）：
   - 水频率的基线 logit 值
   - 典型值：-2 到 5

2. **空间 Range**：
   - 空间相关距离
   - 典型值：100-1000 米

3. **空间 Sigma**：
   - 空间场的变异性（相对深度）
   - 典型值：0.5-5

4. **Calibration Offset**：
   - 将相对深度转换为绝对高程的校准值
   - 应接近湖面平均水位

5. **A-E 曲线形状**：
   - 平缓段 = 盆地平坦
   - 陡峭段 = 湖岸陡峭

---

## 故障排除 | Troubleshooting

### 错误: "CRS mismatch"

**原因**: 栅格数据坐标系不一致

**解决**: 代码会自动重投影，但如果仍有问题，请手动统一坐标系：

```r
library(terra)
dem <- rast("path/to/dem.tif")
water <- rast("path/to/water.tif")

# 重投影 water 到 dem 的坐标系
water <- project(water, dem)
```

### 错误: "INLA not found"

**解决**: 手动安装 INLA：

```r
install.packages("INLA", 
                 repos = "https://inla.r-inla-download.org/R/stable")
```

如果网络问题，从 GitHub 安装：

```r
remotes::install_github("inbo/INLA")
```

### 模型运行时间太长

**解决**: 
1. 减小 mesh 密度（增大 `max_edge`, `cutoff`）
2. 降低栅格分辨率（重采样）
3. 增加 INLA 线程数：

```r
library(INLA)
inla.setOption(num.threads = 4)  # 使用 4 个线程
```

---

## 下一步 | Next Steps

- 查看完整文档：`README_Implementation.md`
- 理解统计模型：`tutorial.md`
- 自定义参数：`config_advanced.R`
- 单步调试：`example_single_lake.R`

---

## 示例输出 | Example Output

运行成功后，您应该看到：

```
================================================================================
  Lake Bathymetry Reconstruction using INLA-SPDE
  湖泊测深图重建与 Area-Elevation 曲线估计
================================================================================

✓ 所有必需包已加载!
✓ 所有模块已加载!

################################################################################
# Processing Lake: Belton (1/1)
################################################################################

====================================
Loading data for Belton...
====================================
   DEM loaded: ...
   Water frequency loaded: ...
   Permanent water mask loaded: ...
   
✓ Data loading and preprocessing complete!

====================================
Building observation data frames...
====================================
   Elevation observations from shoreline: 5423 points
   Water frequency observations: 12847 points
   
✓ Observation data frames built!

====================================
Building SPDE mesh...
====================================
   Number of vertices: 1247
   Number of triangles: 2341
   
✓ Mesh built successfully!

====================================
Fitting INLA model...
====================================
   (This may take several minutes...)
   
✓ Model fitting complete!

====================================
Model Summary
====================================

Fixed Effects:
                       mean     sd
intercept           2.3456  0.1234

Hyperparameters:
                              mean     sd
Range for field            450.23  78.45
Stdev for field              1.23   0.34

Calibration:
Calibration offset (median): 102.456 m

✓✓✓ Lake Belton processing complete! ✓✓✓

================================================================================
  分析完成！
  Analysis Complete!
================================================================================
```

---

**准备好了吗？开始您的分析吧！**

```r
setwd("path/to/LakeAEcurve-INLA")
source("main_runner.R")
```

**祝您分析顺利！ Good luck!**

