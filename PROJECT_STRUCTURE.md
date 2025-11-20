# 项目文件结构说明 | Project File Structure

本文档列出了项目中所有代码文件及其用途。

---

## 核心模块 | Core Modules

### 1. 数据预处理模块（01_Data_Prep/）

#### `data_generation.R`
**功能**: 数据加载与预处理

**主要函数**:
- `load_and_prep_data()`: 
  - 加载三类栅格数据（DEM、水频率、永久水域）
  - 检查 CRS 和分辨率一致性
  - 创建湖区边界和 mask
  - 输出预处理后的数据列表

- `build_observation_data()`:
  - 构建高程观测数据框（Gaussian likelihood）
  - 构建水频率观测数据框（Binomial likelihood）
  - 将水频率转换为 Binomial 计数
  - 处理永久水域点

**输入**: 栅格文件路径
**输出**: 预处理数据列表、观测数据框

---

#### `mesh_setup.R`
**功能**: SPDE mesh 构建与可视化

**主要函数**:
- `build_mesh()`:
  - 基于湖区边界构建 INLA SPDE mesh
  - 使用 `inla.mesh.2d()` 创建三角网格
  - 支持自定义 mesh 参数（max_edge, cutoff, offset）

- `plot_mesh()`:
  - 可视化 mesh 结构
  - 叠加湖区边界和永久水域
  - 显示 mesh 统计信息

**输入**: 预处理数据列表、mesh 参数
**输出**: INLA mesh 对象、mesh 可视化图

---

### 2. 模型实现模块（02_Model_Implementation/）

#### `spde_definition.R`
**功能**: SPDE 模型定义与投影矩阵构建

**主要函数**:
- `define_spde()`:
  - 使用 PC priors 定义 SPDE 模型
  - 调用 `inla.spde2.pcmatern()`
  - 设置空间场先验（range, sigma）

- `build_projection_matrices()`:
  - 构建从 mesh 到观测位置的投影矩阵（A 矩阵）
  - 分别为高程观测和水频率观测创建 A 矩阵

**输入**: mesh 对象、观测数据、先验参数
**输出**: SPDE 对象、投影矩阵列表

---

#### `fit_inlabru_model.R`
**功能**: INLA 模型拟合与结果提取

**主要函数**:
- `build_inla_stacks()`:
  - 构建 INLA 数据栈（多 likelihood 情况）
  - 创建高程 stack（Gaussian）
  - 创建水频率 stack（Binomial）
  - 合并 stacks

- `fit_inla_model()`:
  - 拟合联合 INLA 模型
  - 设置模型公式和先验
  - 使用 `inla()` 进行贝叶斯推断
  - 输出后验参数估计

- `extract_spde_hyperpar()`:
  - 提取空间场超参数（range, sigma）
  - 使用 `inla.spde2.result()` 转换 theta 参数
  - 计算后验统计量

**输入**: SPDE 对象、数据栈、先验设置
**输出**: INLA 拟合结果、超参数后验分布

---

### 3. 结果分析模块（03_Result_Analysis/）

#### `reconstruct_map.R`
**功能**: 湖底 DEM 重建与可视化

**主要函数**:
- `reconstruct_bathymetry()`:
  - 从 INLA 后验场重建栅格 DEM
  - 构建预测网格的 A 矩阵
  - 投影后验均值和标准差到栅格
  - 生成测深图（bathymetry）

- `plot_bathymetry()`:
  - 可视化重建的湖底高程（后验均值）
  - 可视化不确定性（后验标准差）
  - 叠加湖区边界和永久水域

- `plot_original_data()`:
  - 可视化三类原始栅格数据
  - 使用 patchwork 创建组合图

**输入**: INLA 结果、mesh、预处理数据
**输出**: 重建的 DEM 栅格、可视化图

---

#### `ae_curve_ppd.R`
**功能**: Area-Elevation 曲线计算与不确定性分析

**主要函数**:
- `compute_ae_curve()`:
  - 计算 Area-Elevation 曲线
  - 对每个水位高度 h，计算 Z(s) ≤ h 的面积
  - 生成 A-E 曲线数据框

- `plot_ae_curve()`:
  - 绘制 A-E 曲线
  - 美化图形（坐标轴标签、单位转换等）

- `compute_ae_curve_with_uncertainty()`:
  - 基于后验样本计算 A-E 曲线的不确定性
  - 使用 `inla.posterior.sample()` 抽样
  - 计算 95% 置信区间
  - 绘制带置信带的 A-E 曲线

**输入**: 重建的 DEM、预处理数据、INLA 结果
**输出**: A-E 曲线数据框、可视化图

---

## 主运行脚本 | Main Scripts

### `main_runner.R`
**类型**: 主运行脚本（标准版）

**功能**:
- 整合所有模块的完整分析流程
- 自动处理选定的湖泊（Belton, E.V. Spence, 或两者）
- 逐步执行：数据加载 → mesh 构建 → 模型拟合 → 结果分析
- 保存所有输出到 `outputs/` 文件夹

**使用方法**:
```r
source("main_runner.R")
```

**适用场景**: 
- 默认参数分析
- 快速运行整个流程
- 初学者推荐

---

### `main_runner_advanced.R`
**类型**: 主运行脚本（高级版）

**功能**:
- 从配置文件 `config_advanced.R` 读取所有参数
- 支持完全自定义的分析流程
- 包含配置验证和错误检查
- 更灵活的输出控制

**使用方法**:
```r
source("main_runner_advanced.R")
```

**适用场景**:
- 需要自定义参数
- 批量运行不同配置
- 高级用户

---

### `example_single_lake.R`
**类型**: 示例脚本

**功能**:
- 单个湖泊（Belton）的完整分析示例
- 逐步展示每个分析步骤，带详细注释
- 适合学习和调试
- 在屏幕上显示所有中间结果

**使用方法**:
```r
source("example_single_lake.R")
```

**适用场景**:
- 学习代码结构
- 理解每个步骤的输出
- 单步调试
- 教学演示

---

## 配置文件 | Configuration Files

### `config_advanced.R`
**类型**: 高级配置文件

**功能**:
- 集中管理所有模型参数和设置
- 包含详细的参数说明
- 提供配置验证函数
- 支持多湖泊、多场景配置

**主要配置项**:
1. **项目路径**: 数据目录、输出目录
2. **湖泊数据**: 文件路径、启用/禁用
3. **数据预处理**: 缓冲区、试验次数、下采样
4. **Mesh 参数**: max_edge, cutoff, offset
5. **SPDE 先验**: range, sigma 的 PC priors
6. **INLA 参数**: beta/alpha 先验、数值选项
7. **A-E 曲线**: 高程步长、不确定性计算
8. **可视化**: 图形格式、配色方案、DPI
9. **性能**: 线程数、内存管理
10. **输出控制**: 保存选项、文件格式

**辅助函数**:
- `print_config_summary()`: 打印配置摘要
- `validate_config()`: 验证配置有效性

---

## 文档 | Documentation

### `README_Implementation.md`
**类型**: 完整实现文档

**内容**:
- 项目概述和背景
- 统计模型详细说明
- 代码结构和模块说明
- 主要函数 API 文档
- 使用方法和自定义指南
- 输出文件说明
- 模型诊断与解释
- 故障排除
- 扩展与改进建议
- 引用和参考文献

**适用人群**: 
- 需要深入理解项目的用户
- 希望修改或扩展代码的开发者
- 写论文需要引用方法的研究者

---

### `QUICKSTART.md`
**类型**: 快速入门指南

**内容**:
- 三种使用方法（简单、示例、高级）
- 常见问题 FAQ
- 故障排除
- 示例输出展示

**适用人群**:
- 初学者
- 希望快速上手的用户
- 需要快速参考的用户

---

### `tutorial.md`
**类型**: 原始规范文档

**内容**:
- 研究目标与总体思路
- 数据结构与预处理假设
- 潜在空间过程与观测模型（统计建模规范）
- INLA 实现结构要求
- 结果输出与作业可解释性要求

**用途**:
- 项目的"设计文档"
- 理解统计模型的理论基础
- 代码实现的参考规范

---

### `PROJECT_STRUCTURE.md`
**类型**: 本文件

**内容**: 项目文件结构和功能说明

---

## 辅助文件 | Auxiliary Files

### `.gitignore`
**功能**: Git 版本控制忽略文件

**内容**:
- R 临时文件（.Rhistory, .RData）
- 输出文件夹（outputs/）
- 图形文件（*.png, *.pdf）
- 系统文件（.DS_Store, Thumbs.db）
- IDE 配置文件

---

## 数据文件夹 | Data Folders

### `01_Data_Prep/data/`
**结构**:
```
data/
├── DEM/                    # 岸边高程栅格
│   ├── Belton_DEM_100m_Buffer.tif
│   └── EVSpence_DEM_200m_Buffer.tif
├── occurrence/             # 水出现频率栅格
│   ├── Belton_GSW_Occurrence_Buffer100m.tif
│   └── EVSpence_GSW_Occurrence_Buffer200m.tif
└── permanent_water/        # 永久水域矢量
    ├── Belton_Min_Permanent_Water.shp (+ .dbf, .shx, .prj, etc.)
    └── EVSpence_Reservoir_Min_Permanent_Water.shp (+ ...)
```

**说明**: 
- 这些是原始输入数据（未被 git 追踪）
- 确保文件名与配置文件中的路径一致

---

## 输出文件夹 | Output Folder

### `outputs/`
**自动生成的输出**:

#### 可视化图片
- `{lake}_original_data.png`: 原始三类栅格数据
- `{lake}_mesh.png`: SPDE mesh 结构
- `{lake}_bathymetry_mean.png`: 重建的湖底高程
- `{lake}_bathymetry_sd.png`: 高程不确定性
- `{lake}_ae_curve.png`: Area-Elevation 曲线
- `{lake}_ae_curve_uncertainty.png`: 带置信区间的 A-E 曲线（可选）

#### 数据文件
- `{lake}_ae_curve.csv`: A-E 曲线数据表
- `{lake}_bathymetry_mean.tif`: 重建的湖底 DEM（GeoTIFF）
- `{lake}_bathymetry_sd.tif`: 不确定性栅格（GeoTIFF）
- `{lake}_results.RData`: 完整 R 对象（包含模型、数据、结果）

---

## 依赖关系图 | Dependency Graph

```
main_runner.R  或  main_runner_advanced.R
    │
    ├─→ config_advanced.R (仅 advanced 版本)
    │
    ├─→ 01_Data_Prep/data_generation.R
    │   └─→ 包: terra, sf
    │
    ├─→ 01_Data_Prep/mesh_setup.R
    │   └─→ 包: INLA, sf, ggplot2
    │
    ├─→ 02_Model_Implementation/spde_definition.R
    │   └─→ 包: INLA
    │
    ├─→ 02_Model_Implementation/fit_inlabru_model.R
    │   └─→ 包: INLA
    │
    ├─→ 03_Result_Analysis/reconstruct_map.R
    │   └─→ 包: INLA, terra, sf, ggplot2, viridis, patchwork
    │
    └─→ 03_Result_Analysis/ae_curve_ppd.R
        └─→ 包: INLA, terra, ggplot2
```

---

## 典型工作流 | Typical Workflow

### 场景 1：快速分析（默认参数）

```r
# 1. 设置工作目录
setwd("path/to/LakeAEcurve-INLA")

# 2. 运行主脚本
source("main_runner.R")

# 3. 查看结果
# 输出在 outputs/ 文件夹
```

---

### 场景 2：自定义参数分析

```r
# 1. 设置工作目录
setwd("path/to/LakeAEcurve-INLA")

# 2. 修改配置文件
# 打开 config_advanced.R，根据需要修改参数

# 3. 运行高级脚本
source("main_runner_advanced.R")

# 4. 查看结果
list.files("outputs/")
```

---

### 场景 3：学习和调试

```r
# 1. 设置工作目录
setwd("path/to/LakeAEcurve-INLA")

# 2. 运行示例脚本（逐步执行）
source("example_single_lake.R")

# 3. 单步调试（在 RStudio 中）
# 打开 example_single_lake.R
# 逐行运行，观察中间结果
```

---

### 场景 4：加载已保存的结果

```r
# 加载完整结果
load("outputs/Belton_results.RData")

# 查看对象
names(lake_results)
# [1] "lake_name"    "data_list"    "obs_data"     
# [4] "mesh"         "spde"         "result"       
# [7] "hyperpar"     "bathy_result" "ae_df"

# 提取 A-E 曲线
ae_df <- lake_results$ae_df
head(ae_df)

# 重新绘制图形
library(ggplot2)
ggplot(ae_df, aes(x = area_km2, y = elevation)) +
  geom_line() +
  labs(title = "Area-Elevation Curve")
```

---

## 文件大小估计 | File Size Estimates

| 文件类型 | 典型大小 | 说明 |
|---------|---------|------|
| R 脚本 (.R) | 1-10 KB | 代码文件 |
| 文档 (.md) | 10-100 KB | Markdown 文档 |
| 输入栅格 (.tif) | 1-50 MB | 取决于分辨率和范围 |
| 输入矢量 (.shp) | 10-100 KB | Shapefile |
| 输出图片 (.png) | 100-500 KB | 300 DPI |
| 输出栅格 (.tif) | 1-10 MB | 重建的 DEM |
| 结果对象 (.RData) | 10-100 MB | 包含完整 INLA 结果 |

---

## 总代码量统计 | Code Statistics

| 模块 | 文件数 | 函数数 | 总行数（约） |
|------|--------|--------|-------------|
| 数据预处理 | 2 | 4 | 400 |
| 模型实现 | 2 | 5 | 600 |
| 结果分析 | 2 | 7 | 800 |
| 主脚本 | 3 | - | 1500 |
| 配置和文档 | 5 | 2 | 2000 |
| **总计** | **14** | **18** | **~5300** |

---

## 维护和更新 | Maintenance

### 添加新湖泊

1. 准备三类数据文件（DEM、水频率、永久水域）
2. 放置在 `01_Data_Prep/data/` 相应子文件夹
3. 在 `config_advanced.R` 或 `main_runner.R` 中添加湖泊配置：

```r
CONFIG$lakes$new_lake <- list(
  name = "NewLake",
  dem_path = "...",
  water_freq_path = "...",
  perm_water_path = "...",
  enabled = TRUE
)
```

### 修改模型

- 修改先验：编辑 `config_advanced.R` 或相应模块函数
- 改变 likelihood：修改 `fit_inlabru_model.R` 中的模型公式
- 添加协变量：在 `build_inla_stacks()` 中添加固定效应

### 扩展功能

- 添加新的可视化：在 `03_Result_Analysis/` 创建新函数
- 集成其他数据源：修改 `data_generation.R`
- 实现交叉验证：参考 `config_advanced.R` 中的 diagnostics 部分

---

**文档更新日期**: 2025

**项目版本**: 1.0

**维护者**: Cursor AI + User

