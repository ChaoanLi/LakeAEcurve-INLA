################################################################################
# 清除缓存脚本
# Clear Cache Script
#
# 用途：删除所有缓存文件,强制重新处理数据
# Purpose: Remove all cache files to force data reprocessing
################################################################################

cat("\n")
cat("================================================================================\n")
cat("  清除缓存工具 / Cache Clearing Tool\n")
cat("================================================================================\n")
cat("\n")

cache_dir <- "cache"

if (!dir.exists(cache_dir)) {
  cat("✓ 缓存目录不存在,无需清理\n")
  cat("✓ Cache directory does not exist, nothing to clear\n")
  cat("\n")
  quit(save = "no")
}

# 列出所有缓存文件
cache_files <- list.files(cache_dir, full.names = TRUE, recursive = TRUE)

if (length(cache_files) == 0) {
  cat("✓ 缓存目录为空,无需清理\n")
  cat("✓ Cache directory is empty, nothing to clear\n")
  cat("\n")
  quit(save = "no")
}

cat(sprintf("发现 %d 个缓存文件:\n", length(cache_files)))
cat(sprintf("Found %d cache files:\n\n", length(cache_files)))

# 按类型分组显示
rdata_files <- grep("\\.RData$", cache_files, value = TRUE)
tif_files <- grep("\\.tif$", cache_files, value = TRUE)
other_files <- setdiff(cache_files, c(rdata_files, tif_files))

if (length(rdata_files) > 0) {
  cat(sprintf("  RData 文件 (%d):\n", length(rdata_files)))
  for (f in rdata_files) {
    size_mb <- file.size(f) / 1024^2
    cat(sprintf("    - %s (%.2f MB)\n", basename(f), size_mb))
  }
  cat("\n")
}

if (length(tif_files) > 0) {
  cat(sprintf("  TIF 栅格文件 (%d):\n", length(tif_files)))
  for (f in tif_files) {
    size_mb <- file.size(f) / 1024^2
    cat(sprintf("    - %s (%.2f MB)\n", basename(f), size_mb))
  }
  cat("\n")
}

if (length(other_files) > 0) {
  cat(sprintf("  其他文件 (%d):\n", length(other_files)))
  for (f in other_files) {
    size_mb <- file.size(f) / 1024^2
    cat(sprintf("    - %s (%.2f MB)\n", basename(f), size_mb))
  }
  cat("\n")
}

# 计算总大小
total_size_mb <- sum(sapply(cache_files, file.size)) / 1024^2
cat(sprintf("总大小 / Total size: %.2f MB\n\n", total_size_mb))

# 询问用户确认（如果是交互模式）
if (interactive()) {
  cat("是否删除所有缓存文件？\n")
  cat("Delete all cache files?\n")
  cat("  1 - 是/Yes (删除所有)\n")
  cat("  2 - 仅删除 RData 文件 (保留 .tif 栅格)\n")
  cat("  3 - 否/No (取消)\n")
  cat("\n")
  
  choice <- readline("请选择 / Choose (1-3): ")
  
  if (choice == "1") {
    cat("\n正在删除所有缓存文件...\n")
    cat("Deleting all cache files...\n")
    unlink(cache_dir, recursive = TRUE)
    cat("✓ 所有缓存已清除\n")
    cat("✓ All cache cleared\n")
  } else if (choice == "2") {
    cat("\n正在删除 RData 文件...\n")
    cat("Deleting RData files...\n")
    for (f in rdata_files) {
      file.remove(f)
      cat(sprintf("  ✓ 已删除: %s\n", basename(f)))
    }
    cat("\n✓ RData 缓存已清除(TIF 文件保留)\n")
    cat("✓ RData cache cleared (TIF files kept)\n")
  } else {
    cat("\n取消操作\n")
    cat("Operation cancelled\n")
  }
} else {
  # 非交互模式：删除所有
  cat("非交互模式：删除所有缓存\n")
  cat("Non-interactive mode: deleting all cache\n\n")
  unlink(cache_dir, recursive = TRUE)
  cat("✓ 所有缓存已清除\n")
  cat("✓ All cache cleared\n")
}

cat("\n")
cat("================================================================================\n")
cat("完成 / Done\n")
cat("================================================================================\n")
cat("\n")

cat("提示 / Tip:\n")
cat("  下次运行 main_runner.R 时将重新生成缓存\n")
cat("  Cache will be regenerated on next run of main_runner.R\n")
cat("\n")

