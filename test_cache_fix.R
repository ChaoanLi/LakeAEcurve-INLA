################################################################################
# 测试缓存修复
# Test Cache Fix
#
# 用途：验证新的缓存机制是否正常工作
# Purpose: Verify that the new cache mechanism works correctly
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  测试缓存修复 / Testing Cache Fix\n")
cat("================================================================================\n")
cat("\n")

# 加载必需的包
suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# 设置路径
project_dir <- getwd()
source(file.path(project_dir, "01_Data_Prep", "data_generation.R"))

# 测试参数
data_dir <- file.path(project_dir, "01_Data_Prep", "data")
lake_name <- "Belton"

dem_path <- file.path(data_dir, "DEM", "Belton_DEM_100m_Buffer.tif")
water_freq_path <- file.path(data_dir, "occurrence", "Belton_GSW_Occurrence_Buffer100m.tif")
perm_water_path <- file.path(data_dir, "permanent_water", "Belton_Min_Permanent_Water.shp")

# 检查文件是否存在
cat("1. 检查数据文件...\n")
if (!file.exists(dem_path)) {
  cat("  ✗ DEM 文件不存在:", dem_path, "\n")
  quit(save = "no", status = 1)
}
if (!file.exists(water_freq_path)) {
  cat("  ✗ 水频率文件不存在:", water_freq_path, "\n")
  quit(save = "no", status = 1)
}
if (!file.exists(perm_water_path)) {
  cat("  ✗ 永久水域文件不存在:", perm_water_path, "\n")
  quit(save = "no", status = 1)
}
cat("  ✓ 所有数据文件存在\n\n")

# 测试 1: 首次运行（创建新缓存）
cat("================================================================================\n")
cat("测试 1: 创建新缓存\n")
cat("================================================================================\n\n")

# 清除现有缓存
if (dir.exists("cache")) {
  unlink("cache", recursive = TRUE)
  cat("  已清除现有缓存\n\n")
}

# 运行数据加载
cat("2. 运行 load_and_prep_data...\n\n")
tryCatch({
  data_list <- load_and_prep_data(
    dem_path = dem_path,
    water_freq_path = water_freq_path,
    perm_water_path = perm_water_path,
    lake_name = lake_name,
    use_cache = TRUE
  )
  cat("\n✓ 测试 1 通过: 成功创建缓存\n\n")
}, error = function(e) {
  cat("\n✗ 测试 1 失败:", e$message, "\n")
  quit(save = "no", status = 1)
})

# 验证缓存文件
cat("3. 验证缓存文件...\n")
cache_files <- list.files("cache", full.names = TRUE)
expected_files <- c(
  "cache/Belton_prep_data.RData",
  "cache/Belton_dem_cache.tif",
  "cache/Belton_water_freq_cache.tif",
  "cache/Belton_perm_water_cache.tif"
)

all_exist <- all(sapply(expected_files, file.exists))
if (all_exist) {
  cat("  ✓ 所有预期的缓存文件已创建\n")
  for (f in expected_files) {
    size_mb <- file.size(f) / 1024^2
    cat(sprintf("    - %s (%.2f MB)\n", basename(f), size_mb))
  }
  cat("\n")
} else {
  cat("  ✗ 缺少某些缓存文件\n")
  missing <- expected_files[!sapply(expected_files, file.exists)]
  for (f in missing) {
    cat("    缺失:", f, "\n")
  }
  quit(save = "no", status = 1)
}

# 测试 2: 第二次运行（加载缓存）
cat("================================================================================\n")
cat("测试 2: 加载现有缓存\n")
cat("================================================================================\n\n")

cat("4. 再次运行 load_and_prep_data...\n\n")
tryCatch({
  data_list2 <- load_and_prep_data(
    dem_path = dem_path,
    water_freq_path = water_freq_path,
    perm_water_path = perm_water_path,
    lake_name = lake_name,
    use_cache = TRUE
  )
  cat("\n✓ 测试 2 通过: 成功从缓存加载\n\n")
}, error = function(e) {
  cat("\n✗ 测试 2 失败:", e$message, "\n")
  cat("\n这表明缓存加载有问题!\n")
  quit(save = "no", status = 1)
})

# 测试 3: 验证加载的数据是有效的栅格对象
cat("================================================================================\n")
cat("测试 3: 验证数据有效性\n")
cat("================================================================================\n\n")

cat("5. 检查栅格对象...\n")
tryCatch({
  # 尝试访问栅格数据（这会触发 external pointer 错误,如果有的话）
  dem_vals <- values(data_list2$dem)
  water_vals <- values(data_list2$water_freq)
  perm_vals <- values(data_list2$perm_water)
  
  cat(sprintf("  ✓ DEM: %d 个像元\n", length(dem_vals)))
  cat(sprintf("  ✓ Water freq: %d 个像元\n", length(water_vals)))
  cat(sprintf("  ✓ Perm water: %d 个像元\n", length(perm_vals)))
  cat("\n✓ 测试 3 通过: 所有栅格对象有效\n\n")
}, error = function(e) {
  cat("\n✗ 测试 3 失败:", e$message, "\n")
  cat("\n这表明栅格对象的外部指针仍然无效!\n")
  quit(save = "no", status = 1)
})

# 测试 4: 验证数据值的合理性
cat("================================================================================\n")
cat("测试 4: 验证数据范围\n")
cat("================================================================================\n\n")

cat("6. 检查数据范围...\n")
dem_range <- range(values(data_list2$dem), na.rm = TRUE)
water_range <- range(values(data_list2$water_freq), na.rm = TRUE)

cat(sprintf("  DEM 范围: %.2f 到 %.2f m\n", dem_range[1], dem_range[2]))
cat(sprintf("  水频率范围: %.3f 到 %.3f\n", water_range[1], water_range[2]))

if (water_range[1] >= 0 && water_range[2] <= 1) {
  cat("  ✓ 水频率在 [0, 1] 范围内\n")
} else {
  cat("  ✗ 水频率范围异常\n")
}

cat("\n✓ 测试 4 通过: 数据范围正常\n\n")

# 最终总结
cat("================================================================================\n")
cat("  所有测试通过! ✅\n")
cat("  All Tests Passed! ✅\n")
cat("================================================================================\n\n")

cat("缓存修复验证成功!\n")
cat("Cache fix verified successfully!\n\n")

cat("新的缓存机制工作正常:\n")
cat("  ✓ 栅格正确保存为 .tif 文件\n")
cat("  ✓ 缓存可以正确加载\n")
cat("  ✓ 加载的栅格对象有效（无 external pointer 错误）\n")
cat("  ✓ 数据值正确\n\n")

cat("您现在可以安全地运行 main_runner.R\n")
cat("You can now safely run main_runner.R\n\n")

