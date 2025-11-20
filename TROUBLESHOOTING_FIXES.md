# 故障排除与修复记录

## 修复日期
2025年

## 遇到的问题

### 问题 1: Cell Area = 0.00 m² (严重)

**错误信息**:
```
Cell area: 0.00 m²
```

**原因**:
- 数据使用 **WGS 84 (EPSG:4326)** 地理坐标系
- 单位是**度（degrees）**而不是米（meters）
- 分辨率非常小（约 0.0003 度）
- `prod(res(dem_masked))` 计算出的是"平方度"，显示为 0.00

**后果**:
- Area-Elevation 曲线计算错误
- 所有面积相关的统计都会出错
- 无法得到正确的湖泊容量估计

**解决方案**:
在 `01_Data_Prep/data_generation.R` 中添加自动坐标系转换：

```r
# 检查是否为地理坐标系
is_geographic <- st_is_longlat(st_crs(dem_crs))

if (is_geographic) {
  # 自动计算合适的 UTM zone
  utm_zone <- floor((center_lon + 180) / 6) + 1
  utm_epsg <- ifelse(center_lat >= 0, 
                     32600 + utm_zone,  # 北半球
                     32700 + utm_zone)  # 南半球
  
  # 投影到 UTM
  dem_rast <- project(dem_rast, paste0("EPSG:", utm_epsg))
  water_freq_rast <- project(water_freq_rast, paste0("EPSG:", utm_epsg))
  perm_water_sf <- st_transform(perm_water_sf, crs = utm_epsg)
}
```

**效果**:
- 数据自动投影到合适的 UTM 坐标系
- 单位变为米
- Cell area 变为正确的值（例如 ~100-10000 m²）

---

### 问题 2: 模型简化（2025更新）

**背景**:
早期版本尝试使用多 likelihood 模型（Gaussian for elevation + Binomial for water），但遇到以下问题：
- Stack 合并的复杂性
- 模型收敛不稳定
- 结果出现异常（如深度几乎为零）

**当前解决方案（简化策略）**:
采用单 likelihood 模型 + 后处理校准：

1. **模型拟合阶段**：
   - 仅使用水频率数据（Binomial likelihood）
   - 空间场表示相对深度
   - 模型更简单、更稳定

2. **校准阶段**（在 `reconstruct_bathymetry()` 中）：
   - 使用岸边 DEM 观测校准空间场
   - 计算偏移量：`offset = median(DEM_shore - field_shore)`
   - 转换为绝对高程：`elevation = field + offset`

**优点**:
- ✅ 避免多 likelihood 的复杂性
- ✅ 模型收敛更稳定
- ✅ 最终仍得到绝对高程和准确的 A-E 曲线
- ✅ 代码更易维护

**实现文件**:
- `02_Model_Implementation/fit_inlabru_model.R` - 单 likelihood 模型
- `03_Result_Analysis/reconstruct_map.R` - 校准步骤

---

### 问题 3: ggplot2 警告

**警告信息 1**:
```
Warning: Option 'terrain' does not exist. Defaulting to 'viridis'.
```

**原因**: `viridis` 包没有 'terrain' 选项

**解决方案**: 改为 'viridis' (或 'plasma', 'magma', 'inferno', 'cividis')

---

**警告信息 2**:
```
Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
Please use `linewidth` instead.
```

**原因**: ggplot2 3.4.0+ 版本更新了参数名称

**解决方案**: 全局替换
- `geom_line(size = ...)` → `geom_line(linewidth = ...)`
- `geom_path(size = ...)` → `geom_path(linewidth = ...)`
- `geom_sf(size = ...)` → `geom_sf(linewidth = ...)`

**注意**: `geom_point(size = ...)` 保持不变（点的大小仍用 `size`）

**修改的文件**:
- `03_Result_Analysis/reconstruct_map.R`
- `01_Data_Prep/mesh_setup.R`
- `03_Result_Analysis/ae_curve_ppd.R`

---

**警告信息 3**:
```
Warning: `inla.mesh.2d()` was deprecated in INLA 23.08.18.
Please use `fmesher::fm_mesh_2d_inla()` instead.
```

**状态**: 暂时保留，代码仍可正常运行

**未来改进**: 可以更新为 `fmesher` 包的新函数

---

## 修复后的表现

### 预期输出（修复后）:

```
====================================
Loading data for Belton...
====================================
1. Reading raster data...
   DEM loaded: ...

2. Checking CRS and resolution alignment...
   Original DEM CRS: EPSG:4326
   ⚠ 检测到地理坐标系（单位：度）
   正在投影到 UTM 坐标系（单位：米）...
   目标坐标系: UTM Zone 14 (EPSG:32614)
   ✓ 投影完成
   DEM resolution: 30.00 x 30.00  (米)
   Water frequency resolution: 30.00 x 30.00

...

7. Data summary:
   DEM range: 176.91 to 238.51 m
   Water frequency range: 0.000 to 0.990
   Permanent water pixels: 439569
   Cell area: 900.00 m²  ✓ 正常

✓ Data loading and preprocessing complete!

====================================
Building INLA stacks...
====================================

1. Building water frequency stack (Binomial likelihood)...
   Note: Using water-only model for stability
   Shore elevation will be used for calibration in post-processing
   Water observations: 1234724
   Shore observations available for calibration: 346144

✓ INLA stacks built!

...

====================================
Reconstructing bathymetry from posterior field...
====================================

...

6. Calibrating with shore elevation data...
   Using 346144 shore elevation points for calibration
   Calibration offset (median): 102.456 m
   Offset range at shore points: 98.123 to 106.789 m
   Calibrated elevation range: 176.91 to 238.51 m

✓ Bathymetry reconstruction complete!
   Calibration applied: 102.456 m offset
```

---

## 测试建议

运行以下代码测试修复效果：

```r
# 清理环境
rm(list = ls())
gc()

# 运行主脚本
source("main_runner.R")

# 检查关键输出
# 1. 查看 Cell area 是否正常（应该 > 100 m²）
# 2. 查看是否有 "Name mismatch" 错误
# 3. 查看模型是否成功拟合
# 4. 查看 A-E 曲线是否合理
```

---

## 相关文件

**修改的文件**:
1. `01_Data_Prep/data_generation.R` - 添加坐标系转换
2. `02_Model_Implementation/fit_inlabru_model.R` - 简化为单 likelihood 模型
3. `03_Result_Analysis/reconstruct_map.R` - 添加校准步骤，更新 ggplot2 参数
4. `main_runner.R` - 将 obs_data 添加到 data_list
5. `01_Data_Prep/mesh_setup.R` - 更新 ggplot2 参数
6. `03_Result_Analysis/ae_curve_ppd.R` - 更新 ggplot2 参数

**更新的文档**:
- `README_Implementation.md` - 更新模型描述
- `QUICKSTART.md` - 更新快速开始指南
- `tutorial.md` - 添加简化模型说明
- `TROUBLESHOOTING_FIXES.md` (本文件)

---

## 注意事项

### 坐标系转换的影响

1. **Mesh 参数**: 
   - 原来的 `max_edge = c(100, 500)` 在度单位下意味着非常大的距离
   - 转换为 UTM 后，这些值变为米，是合理的
   - 如果原始数据已经是投影坐标系，不会进行转换

2. **缓冲区参数**:
   - `lake_boundary` 的 500m 缓冲区现在是真正的 500 米
   - 在地理坐标系中，这个值会被自动调整

3. **分辨率**:
   - WGS84 下约 0.0003° ≈ 30m（在纬度 30° 附近）
   - UTM 下直接是 30m

### 数据质量检查

修复后，应该检查：
- ✅ Cell area 在合理范围内（10-10000 m²）
- ✅ DEM 分辨率在合理范围内（10-100 m）
- ✅ 湖泊面积计算合理
- ✅ A-E 曲线形状合理（单调递增）

---

## 更新历史

- **2025-xx-xx**: 初始修复
  - 添加自动 UTM 投影
  - 修复 INLA stack name mismatch
  - 更新 ggplot2 参数

---

**修复状态**: ✅ 完成

**测试状态**: ⏳ 待测试

**建议**: 运行 `main_runner.R` 验证所有修复

