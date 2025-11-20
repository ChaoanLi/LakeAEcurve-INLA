# 湖泊测深图重建与 Area-Elevation 曲线估计

## Lake Bathymetry Reconstruction and Area-Elevation Curve Estimation using INLA-SPDE

---

## 项目概述 | Project Overview

本项目实现了一个完整的 R 管道，使用 **INLA（Integrated Nested Laplace Approximation）** 和 **SPDE（Stochastic Partial Differential Equation）** 方法，从多源栅格数据重建湖底测深图（bathymetry），并计算湖泊的 **Area-Elevation (A-E) 曲线**。

### 输入数据 | Input Data

项目需要三类栅格数据：

1. **岸边 DEM（Shoreline Digital Elevation Model）**
   - 提供陆地和湖岸区域的高程观测
   - 作为 Gaussian likelihood 的观测数据

2. **水出现频率栅格（Water Occurrence Frequency）**
   - 每个像元在观测期内被水覆盖的频率（0-1 或 0-100）
   - 转换为 Binomial 观测数据

3. **永久水域 mask（Permanent Water Mask）**
   - 标识始终有水的区域（湖中最深部分）
   - 用于强约束模型

### 统计模型 | Statistical Model

采用层次空间模型框架：

**潜在场（Latent Field）**：
- 相对深度场 Z(s) ~ GRF（高斯随机场）
- 使用 Matérn 协方差结构（通过 SPDE-GMRF 近似）

**观测模型（Observation Model）**：

**水频率观测（Binomial + logit）**：
```
Y_water(s) ~ Binomial(N, p(s))
logit(p(s)) = α + f(s)
```
- 含义：空间场 f(s) 表征相对深度（越深越容易被水淹没）

**校准步骤（Calibration）**：
- 使用岸边 DEM 观测值校准相对深度场，得到绝对高程：
```
Elevation(s) = f(s) + offset
offset = median(DEM_shore - f_predicted_shore)
```

---

## 代码结构 | Code Structure

```
LakeAEcurve-INLA/
│
├── 01_Data_Prep/              # 数据预处理模块
│   ├── data_generation.R      # 加载栅格、构建观测数据
│   ├── mesh_setup.R           # SPDE mesh 构建与可视化
│   └── data/                  # 原始数据文件夹
│       ├── DEM/
│       ├── occurrence/
│       └── permanent_water/
│
├── 02_Model_Implementation/   # 模型实现模块
│   ├── spde_definition.R      # SPDE 模型定义、投影矩阵
│   └── fit_inlabru_model.R    # INLA 模型拟合、超参数提取
│
├── 03_Result_Analysis/        # 结果分析模块
│   ├── reconstruct_map.R      # 重建湖底 DEM、可视化
│   └── ae_curve_ppd.R         # A-E 曲线计算与不确定性分析
│
├── main_runner.R              # 主运行脚本（整合所有模块）
├── README_Implementation.md   # 本文件
└── outputs/                   # 输出文件夹（自动创建）
```

---

## 安装依赖 | Installation

### 必需的 R 包 | Required R Packages

```r
# 基础空间数据处理
install.packages(c("terra", "sf", "sp"))

# 可视化
install.packages(c("ggplot2", "viridis", "patchwork"))

# INLA（需要从特定源安装）
install.packages("INLA", 
                 repos = c(getOption("repos"), 
                           INLA = "https://inla.r-inla-download.org/R/stable"), 
                 dep = TRUE)
```

### 检查安装 | Verify Installation

```r
library(terra)
library(sf)
library(INLA)
library(ggplot2)

cat("All packages loaded successfully!\n")
```

---

## 使用方法 | Usage

### 快速开始 | Quick Start

1. **准备数据**：确保数据文件放在 `01_Data_Prep/data/` 目录下

2. **运行主脚本**：
   ```r
   source("main_runner.R")
   ```

3. **查看结果**：所有输出保存在 `outputs/` 文件夹

### 自定义参数 | Customization

在 `main_runner.R` 中可以修改：

#### 选择湖泊
```r
# 第 2 部分：选择要分析的湖泊
lake_choice <- 1  # 1 = Belton, 2 = E.V. Spence, 3 = Both
```

#### 调整 Mesh 参数
```r
# 第 3.4 部分：构建 SPDE mesh
mesh <- build_mesh(
  data_list = data_list,
  max_edge = c(100, 500),    # 内部/外部三角形最大边长
  cutoff = 50,               # 最小节点间距
  offset = c(100, 500)       # 边界扩展距离
)
```

#### 修改先验
```r
# 第 3.5 部分：定义 SPDE 模型
spde <- define_spde(
  mesh = mesh,
  prior_range = c(500, 0.5),   # P(range < 500) = 0.5
  prior_sigma = c(1, 0.01)      # P(sigma > 1) = 0.01
)
```

#### 设置 Beta 系数先验
```r
# 第 3.8 部分：拟合 INLA 模型
result <- fit_inla_model(
  stack_list = stack_list,
  spde = spde,
  obs_data = obs_data,
  use_elev = TRUE,
  beta_prior = c(0, 0.01)      # Beta ~ N(mean=0, precision=0.01)
)
```

#### A-E 曲线精度
```r
# 第 3.11 部分：计算 Area-Elevation 曲线
ae_df <- compute_ae_curve(
  bathy_result = bathy_result,
  data_list = data_list,
  elevation_step = 0.1         # 高程步长（米）
)
```

---

## 主要函数说明 | Main Functions

### 1. 数据预处理 (`data_generation.R`)

#### `load_and_prep_data()`
- **功能**：加载和预处理三类栅格数据
- **输入**：DEM、水频率、永久水域 shapefile 路径
- **输出**：包含所有预处理数据的列表
- **检查**：CRS 一致性、分辨率对齐、值域规范化

#### `build_observation_data()`
- **功能**：构建 INLA 观测数据框
- **输入**：预处理后的数据列表
- **输出**：高程观测 (`elev_df`) 和水频率观测 (`water_df`)

### 2. Mesh 构建 (`mesh_setup.R`)

#### `build_mesh()`
- **功能**：构建 SPDE 三角网格
- **参数**：
  - `max_edge`：三角形最大边长 [内部, 外部]
  - `cutoff`：最小节点间距
  - `offset`：边界扩展 [内部, 外部]
- **输出**：INLA mesh 对象

#### `plot_mesh()`
- **功能**：可视化 mesh 结构

### 3. SPDE 定义 (`spde_definition.R`)

#### `define_spde()`
- **功能**：定义 SPDE 模型（使用 PC priors）
- **参数**：
  - `prior_range`：空间相关距离先验
  - `prior_sigma`：边际标准差先验

#### `build_projection_matrices()`
- **功能**：构建观测位置到 mesh 的投影矩阵（A 矩阵）

### 4. 模型拟合 (`fit_inlabru_model.R`)

#### `build_inla_stacks()`
- **功能**：构建 INLA 数据栈（使用 water 数据）
- **策略**：单 likelihood 模型，岸边高程用于后处理校准

#### `fit_inla_model()`
- **功能**：拟合 Binomial-SPDE 模型
- **输出**：后验参数估计、空间场分布
- **解释**：空间场表示相对深度（非绝对高程）

#### `extract_spde_hyperpar()`
- **功能**：提取空间场超参数（range, sigma）

### 5. 结果分析 (`reconstruct_map.R`, `ae_curve_ppd.R`)

#### `reconstruct_bathymetry()`
- **功能**：从后验场重建栅格 DEM
- **关键步骤**：使用岸边 DEM 观测校准相对深度场，得到绝对高程
- **输出**：校准后的后验均值和标准差栅格

#### `plot_bathymetry()`
- **功能**：可视化重建的湖底高程和不确定性

#### `compute_ae_curve()`
- **功能**：计算 Area-Elevation 曲线
- **方法**：对每个水位 h，计算 Z(s) ≤ h 的面积

#### `compute_ae_curve_with_uncertainty()`
- **功能**：基于后验样本计算 A-E 曲线的 95% 置信区间

---

## 输出文件 | Output Files

运行完成后，`outputs/` 文件夹包含：

### 可视化图片

| 文件名 | 内容 |
|--------|------|
| `*_original_data.png` | 原始三类栅格数据可视化 |
| `*_mesh.png` | SPDE mesh 结构与湖区边界 |
| `*_bathymetry_mean.png` | 重建的湖底高程（后验均值） |
| `*_bathymetry_sd.png` | 湖底高程的不确定性（后验标准差） |
| `*_ae_curve.png` | Area-Elevation 曲线 |
| `*_ae_curve_uncertainty.png` | 带 95% 置信区间的 A-E 曲线（可选） |

### 数据文件

| 文件名 | 格式 | 内容 |
|--------|------|------|
| `*_ae_curve.csv` | CSV | A-E 曲线数据表 |
| `*_results.RData` | R 对象 | 完整分析结果（模型、后验、栅格等） |

---

## 模型诊断与解释 | Model Diagnostics & Interpretation

### 检查模型拟合质量

1. **固定效应**：
   - `Intercept`：水频率的基线 logit 值

2. **空间场超参数**：
   - `Range`：空间相关距离（米）
     - 典型值：100-1000 米（取决于湖泊大小）
   - `Sigma`：空间场的边际标准差（相对深度单位）
     - 典型值：0.5-5（取决于深度变异性）

3. **校准质量**：
   - 检查校准偏移量（calibration offset）是否合理
   - 查看岸边点的拟合残差范围
   - 典型偏移量应接近湖面平均水位

4. **不确定性评估**：
   - 查看 `*_bathymetry_sd.png`
   - 高不确定性区域通常在：
     - 远离岸边的深水区
     - 水频率信息稀疏的区域

### A-E 曲线解释

- **曲线形状**：
  - 平缓段：盆地平坦区域（单位高程变化对应大面积变化）
  - 陡峭段：湖岸陡壁（单位高程变化对应小面积变化）

- **应用**：
  - 水库容量估算
  - 不同水位下的湖泊面积预测
  - 洪水/干旱情景模拟

---

## 故障排除 | Troubleshooting

### 常见问题

1. **INLA 安装失败**
   ```r
   # 尝试手动安装
   install.packages("INLA", 
                    repos = "https://inla.r-inla-download.org/R/stable")
   ```

2. **CRS 不匹配错误**
   - 确保所有栅格和矢量数据在同一坐标系
   - 代码会自动重投影，但建议预先统一

3. **内存不足**
   - 减小 mesh 密度（增大 `max_edge`，增大 `cutoff`）
   - 降低栅格分辨率（重采样）
   - 减少后验样本数（`n_samples`）

4. **模型不收敛**
   - 调整先验参数
   - 检查数据质量（是否有异常值）
   - 增加观测数据（特别是岸边高程点）

5. **Beta 系数为负**
   - 可能是数据问题（水频率与高程关系颠倒）
   - 检查 DEM 是否正确（应该是绝对高程）
   - 可以在模型中强制 Beta > 0（通过 log 变换）

---

## 扩展与改进 | Extensions & Improvements

### 可能的改进方向

1. **非平稳空间场**：
   - 使用 `inla.barrier.pcmatern` 处理陆地障碍

2. **时间动态模型**：
   - 如果有多时相数据，可加入时间维度

3. **多湖泊联合建模**：
   - 共享超参数，提高小湖泊估计精度

4. **贝叶斯模型选择**：
   - 比较不同协方差函数（Matérn vs Exponential）
   - 使用 DIC、WAIC 选择最优模型

5. **集成其他数据源**：
   - 卫星测高数据（Altimetry）
   - 实测水深点（In-situ measurements）

---

## 引用与参考 | Citation & References

### 本项目基于以下方法论：

1. **INLA 方法**：
   - Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B*, 71(2), 319-392.

2. **SPDE 方法**：
   - Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *Journal of the Royal Statistical Society: Series B*, 73(4), 423-498.

3. **空间数据融合**：
   - Banerjee, S., Carlin, B. P., & Gelfand, A. E. (2014). *Hierarchical modeling and analysis for spatial data*. CRC press.

### R 包引用

```r
citation("INLA")
citation("terra")
citation("sf")
```

---

## 许可证 | License

本项目代码采用 MIT License。

数据来源请查看原始数据提供方的许可协议。

---

## 联系与支持 | Contact & Support

如有问题或建议，请联系项目维护者。

**项目创建日期**：2025

**最后更新**：2025

---

## 致谢 | Acknowledgments

感谢：
- INLA 开发团队
- R-spatial 社区
- STAT 647 课程团队

---

**祝您分析顺利！Good luck with your analysis!**

