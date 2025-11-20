# 缓存功能使用说明

## 什么是缓存？

为了加速重复运行，我们为以下慢速步骤添加了**自动缓存**功能：

1. **数据预处理** (`load_and_prep_data`) - 包括：
   - 读取栅格
   - CRS 转换和投影
   - 栅格重采样
   - 裁剪和 mask
   
2. **观测数据构建** (`build_observation_data`) - 包括：
   - 提取岸边高程点
   - 构建水频率观测点
   - Binomial 计数转换

---

## 如何使用

### 自动模式（推荐）✅

**默认已启用**！无需任何修改：

```r
source("main_runner.R")
```

**第一次运行**：
- 正常处理数据（较慢）
- 自动保存缓存到 `cache/` 文件夹

**第二次及以后**：
- 自动检测并加载缓存
- **跳过耗时的预处理步骤**
- 极速启动！⚡

---

## 输出示例

### 第一次运行（创建缓存）

```
====================================
Loading data for Belton...
====================================
1. Reading raster data...
   DEM loaded: ...
   
2. Checking CRS and resolution alignment...
   ⚠ 检测到地理坐标系（单位：度）
   正在投影到 UTM 坐标系（单位：米）...
   ✓ 投影完成
   
... (完整的预处理步骤)

✓ Data loading and preprocessing complete!

Saving preprocessed data to cache...
✓ Cache saved: cache/Belton_prep_data.RData

====================================
Building observation data frames...
====================================
... (完整的观测数据构建)

✓ Observation data frames built!

Saving observation data to cache...
✓ Cache saved: cache/Belton_obs_data.RData
```

### 第二次运行（使用缓存）⚡

```
====================================
Loading cached data for Belton...
====================================
   Cache file: cache/Belton_prep_data.RData
✓ Cached data loaded successfully!
   Skipping preprocessing (saved time!)

====================================
Loading cached observation data...
====================================
   Cache file: cache/Belton_obs_data.RData
✓ Cached observation data loaded successfully!
   Skipping observation data construction (saved time!)
```

**速度提升**: 从几分钟 → 几秒钟！

---

## 缓存文件

缓存保存在项目根目录的 `cache/` 文件夹：

```
cache/
├── Belton_prep_data.RData         # Belton 元数据
├── Belton_dem_cache.tif           # Belton DEM 栅格
├── Belton_water_freq_cache.tif    # Belton 水频率栅格
├── Belton_perm_water_cache.tif    # Belton 永久水域栅格
├── Belton_obs_data.RData          # Belton 观测数据
├── EVSpence_prep_data.RData       # EV Spence 元数据
├── EVSpence_dem_cache.tif         # EV Spence DEM 栅格
├── EVSpence_water_freq_cache.tif  # EV Spence 水频率栅格
├── EVSpence_perm_water_cache.tif  # EV Spence 永久水域栅格
└── EVSpence_obs_data.RData        # EV Spence 观测数据
```

**文件大小**: 通常 10-100 MB（取决于数据量）

**注意**: 栅格数据保存为独立的 `.tif` 文件以避免外部指针失效问题。

---

## 何时需要清除缓存

如果您修改了：
- ✅ 原始数据文件（DEM、水频率、永久水域）
- ✅ 预处理参数（如分辨率、缓冲区大小）
- ✅ 观测数据参数（如 `N_trials`）

**清除方法**:

### 方法 1: 删除缓存文件夹

```r
# 在 R 中
unlink("cache", recursive = TRUE)
```

或者手动删除 `cache/` 文件夹

### 方法 2: 禁用缓存并重新运行

```r
# 在函数调用时
data_list <- load_and_prep_data(
  dem_path = ...,
  water_freq_path = ...,
  perm_water_path = ...,
  lake_name = "Belton",
  use_cache = FALSE  # 禁用缓存
)
```

### 方法 3: 删除特定湖泊的缓存

```r
# 仅删除 Belton 的缓存
file.remove("cache/Belton_prep_data.RData")
file.remove("cache/Belton_obs_data.RData")
```

---

## 高级用法

### 自定义缓存目录

```r
data_list <- load_and_prep_data(
  dem_path = ...,
  water_freq_path = ...,
  perm_water_path = ...,
  lake_name = "Belton",
  cache_dir = "my_custom_cache"  # 自定义缓存文件夹
)
```

### 检查缓存是否存在

```r
# 检查 Belton 的缓存
cache_file <- "cache/Belton_prep_data.RData"
if (file.exists(cache_file)) {
  cat("缓存存在\n")
  # 查看文件大小
  cat(sprintf("大小: %.2f MB\n", file.size(cache_file) / 1024^2))
  # 查看修改时间
  cat(sprintf("修改时间: %s\n", file.mtime(cache_file)))
} else {
  cat("缓存不存在，将在下次运行时创建\n")
}
```

### 在主脚本中完全禁用缓存

编辑 `main_runner.R`：

```r
# ---- 3.1 加载和预处理数据 ----
data_list <- load_and_prep_data(
  dem_path = lake_info$dem_path,
  water_freq_path = lake_info$water_freq_path,
  perm_water_path = lake_info$perm_water_path,
  lake_name = lake_name,
  use_cache = FALSE  # 添加这一行
)

# ---- 3.2 构建观测数据 ----
obs_data <- build_observation_data(
  data_list = data_list,
  N_trials = 100,
  use_shore_elev = TRUE,
  use_cache = FALSE  # 添加这一行
)
```

---

## 故障排除

### 问题 1: "external pointer is not valid" ⚠️

**症状**: 
```
错误于eval(ei, envir): external pointer is not valid
```

**原因**: 旧版本的缓存直接保存了 terra 栅格对象,其外部指针在重新加载后失效

**解决**: 删除旧缓存文件

```r
# 删除所有预处理缓存
file.remove(list.files("cache", pattern = "_prep_data.RData$", full.names = TRUE))

# 或删除所有缓存
unlink("cache", recursive = TRUE)
```

重新运行后将使用新的缓存格式(栅格保存为 .tif 文件)。

**状态**: ✅ 已在新版本中修复

### 问题 2: "Error in load(...)"

**原因**: 缓存文件损坏

**解决**: 删除该缓存文件并重新运行

```r
file.remove("cache/Belton_prep_data.RData")
```

### 问题 3: 数据更新后结果没变

**原因**: 仍在使用旧缓存

**解决**: 清除所有缓存

```r
unlink("cache", recursive = TRUE)
```

### 问题 4: 磁盘空间不足

**原因**: 缓存文件占用空间(新版本会保存额外的 .tif 文件)

**解决**: 定期清理不需要的缓存

```r
# 查看缓存文件大小
list.files("cache", full.names = TRUE, recursive = TRUE) |>
  sapply(file.size) |>
  sum() / 1024^2  # MB

# 删除所有缓存
unlink("cache", recursive = TRUE)
```

---

## 缓存的优势

✅ **速度**: 第二次运行快 10-100 倍  
✅ **自动**: 无需手动操作  
✅ **智能**: 自动检测是否存在  
✅ **灵活**: 随时可以禁用或清除  
✅ **安全**: 不影响原始数据

---

## 注意事项

⚠️ **缓存不包含**:
- INLA 模型拟合结果（每次都要重新拟合）
- Mesh 构建（较快，每次重新构建）
- 可视化图形

⚠️ **何时不使用缓存**:
- 正在调试数据预处理代码
- 需要测试不同的预处理参数
- 原始数据已更新

⚠️ **注意事项**:
- 缓存文件可能较大（10-100 MB）
- 确保有足够的磁盘空间
- `.gitignore` 已配置忽略 `cache/` 文件夹

---

## 总结

- ✅ **默认启用**: 自动加速，无需配置
- ✅ **第一次慢**: 创建缓存
- ✅ **之后快**: 使用缓存
- ✅ **数据变化**: 清除缓存
- ✅ **灵活控制**: 可随时禁用

**享受极速分析体验！** ⚡


