# 模型改进已集成到主流程

## 改进内容

已将大坝点约束和真实DEM验证功能直接集成到主程序中，无需额外脚本。

---

## 修改的文件

### 1. `01_Data_Prep/data_generation.R`

**`build_observation_data()` 函数**

新增参数：
- `dam_point_path` - 大坝点shapefile路径（可选）
- `dam_point_weight` - 大坝点重复次数（默认50），用于增强约束

**功能：**
- 自动读取大坝点
- 从DEM提取大坝位置的高程
- 将大坝点重复多次添加到高程观测中
- 作为最低点约束，改进深水区估计

---

### 2. `03_Result_Analysis/ae_curve_ppd.R`

**`plot_ae_curve()` 函数**

新增参数：
- `true_ae_path` - 真实A-E曲线CSV路径（来自04_Validation文件夹，可选）

**功能：**
- 读取验证数据CSV文件（格式：Belton_AVE.csv, EVspence_AVE.csv）
- 自动提取Elevation (m)和Area (km2)列
- 在同一图中绘制预测和真实曲线
  - **黑色实线** = 真实曲线（from validation data）
  - **蓝色虚线** = 预测曲线（INLA-SPDE）
- 自动计算并输出误差统计（MAE, RMSE, MAPE）

---

### 3. `main_runner.R`

**湖泊信息扩展**

每个湖泊配置新增：
- `dam_point_path` - 大坝点路径
- `true_ae_path` - 真实A-E曲线CSV路径（用于验证）

**调用更新：**

```r
# 构建观测数据（自动包含大坝点）
obs_data <- build_observation_data(
  data_list = data_list,
  N_trials = 100,
  use_shore_elev = TRUE,
  dam_point_path = lake_info$dam_point_path,  # 新增
  dam_point_weight = 50                       # 新增
)

# 绘制A-E曲线（自动对比真实验证数据）
plot_ae_curve(
  ae_df = ae_df,
  data_list = data_list,
  save_path = file.path(output_dir, sprintf("%s_ae_curve.png", lake_name)),
  true_ae_path = lake_info$true_ae_path  # 新增：真实A-E曲线验证数据
)
```

---

## 使用方法

### 直接运行主程序

```r
source("main_runner.R")
```

**程序会自动：**
1. ✅ 读取大坝点并提取高程
2. ✅ 将大坝点作为强约束添加到模型中
3. ✅ 使用约束后的模型拟合
4. ✅ 从真实DEM计算真实A-E曲线
5. ✅ 在同一图中绘制预测和真实曲线
6. ✅ 输出误差统计

---

## 输出结果

### 控制台输出示例

```
====================================
Building observation data frames...
====================================

1. Building elevation observation data...
   Elevation observations from shoreline: 5423 points
   Elevation range: 180.23 to 238.51 m

   Adding dam point as lowest elevation constraint...
   Dam point location: (-97.456789, 31.123456)
   Dam elevation from DEM: 176.45 m
   Replicating 50 times to increase constraint weight
   ✓ Total elevation observations: 5473 (including 50 dam replicates)

...

Plotting Area-Elevation curve...
   Loading true A-E curve from validation data...
   Loaded 1134 true A-E data points
   True elevation range: 146.74 to 238.51 m
   Validation errors (in common elevation range):
     MAE  = 0.1234 km²
     RMSE = 0.2345 km²
     MAPE = 2.34%
   A-E curve plot saved to: outputs/Belton_ae_curve.png
```

### 图片输出

`outputs/Belton_ae_curve.png` 包含：
- **黑色实线** - 真实A-E曲线（参考标准）
- **蓝色虚线** - 预测A-E曲线（模型结果）
- 图例区分两条曲线
- 标题和坐标轴标签清晰

---

## 改进效果

### 预期改进

通过添加大坝点约束：
- ✅ 最低点（大坝位置）高程更准确
- ✅ 深水区估计改善
- ✅ A-E曲线整体更接近真实值
- ✅ 误差指标（MAE, RMSE, MAPE）降低

### 验证方法

查看控制台输出的误差统计：
- **MAE** < 0.5 km² → 很好
- **RMSE** < 1.0 km² → 很好
- **MAPE** < 5% → 很好

---

## 参数调整

### 增加大坝点约束权重

如果想让大坝点约束更强：

在 `main_runner.R` 中修改：

```r
obs_data <- build_observation_data(
  data_list = data_list,
  N_trials = 100,
  use_shore_elev = TRUE,
  dam_point_path = lake_info$dam_point_path,
  dam_point_weight = 100  # 从50改为100
)
```

### 禁用大坝点约束

如果想看没有约束的效果（用于对比）：

```r
obs_data <- build_observation_data(
  data_list = data_list,
  N_trials = 100,
  use_shore_elev = TRUE,
  dam_point_path = NULL  # 设为NULL禁用
)
```

---

## 扩展到其他湖泊

代码已自动配置好两个湖泊：
- Belton Lake
- E.V. Spence Reservoir

在 `main_runner.R` 中设置：

```r
lake_choice <- 1  # 1 = Belton
lake_choice <- 2  # 2 = E.V. Spence
lake_choice <- 3  # 3 = 两个都运行
```

---

## 技术细节

### 大坝点约束原理

1. **提取高程**：从DEM提取大坝位置的真实高程
2. **重复观测**：将该点重复50次（相当于增加权重）
3. **空间传播**：SPDE的平滑性让约束影响周围区域
4. **锚定作用**：将相对深度场"锚定"到正确的高程范围

### 真实A-E曲线对比原理

1. **验证数据加载**：从CSV文件（04_Validation/）读取真实A-E曲线
2. **数据格式**：自动提取Elevation (m)和Area (km2)列
3. **高程范围对齐**：在预测和真实曲线的重叠高程范围内计算误差
4. **插值对比**：插值到相同的高程点进行点对点对比
5. **误差统计**：计算MAE、RMSE、MAPE指标
6. **可视化对比**：在同一图中展示两条曲线

---

## 故障排除

### 问题1：找不到大坝点文件

**错误：** 大坝点文件不存在

**解决：** 确保文件路径正确

```r
file.exists("01_Data_Prep/data/DAM_point/Belton_dam.shp")  # 应返回TRUE
```

### 问题2：大坝点高程为NA

**错误：** 大坝点不在DEM范围内

**解决：** 检查CRS是否一致

```r
library(sf)
library(terra)

dam <- st_read("01_Data_Prep/data/DAM_point/Belton_dam.shp")
dem <- rast("01_Data_Prep/data/DEM/Belton_DEM_100m_Buffer.tif")

print(st_crs(dam))
print(crs(dem))
```

### 问题3：缓存数据问题

如果修改了参数但结果没变化，可能是使用了旧缓存。

**解决：** 删除观测数据缓存

```r
file.remove("cache/Belton_obs_data.RData")
```

然后重新运行 `source("main_runner.R")`

---

## 重要提示

⚠️ **已清除观测数据缓存**

由于修改了 `build_observation_data()` 函数，首次运行会重新生成观测数据（包含大坝点约束）。这是正常的，不会影响其他缓存（DEM、water frequency等）。

---

**现在直接运行 `main_runner.R` 即可获得改进的结果！** 🚀

```r
source("main_runner.R")
```

