################################################################################
# 重新绘制 A-E 曲线（包含真实验证数据对比）
# Replot A-E curve with validation data comparison
################################################################################

cat("================================================================================\n")
cat("  重新绘制 A-E 曲线\n")
cat("  Replotting A-E Curve with Validation Data\n")
cat("================================================================================\n\n")

# 加载必要的库和函数
source("03_Result_Analysis/ae_curve_ppd.R")

# 选择湖泊
lake_name <- "Belton"  # 或 "EVSpence"

cat(sprintf("处理湖泊: %s\n\n", lake_name))

# 设置路径
results_file <- sprintf("outputs/%s_results.RData", lake_name)
ae_csv_file <- sprintf("outputs/%s_ae_curve.csv", lake_name)
true_ae_file <- sprintf("04_Validation/%s_AVE.csv", lake_name)

# 检查文件是否存在
if (!file.exists(results_file)) {
  stop(sprintf("找不到结果文件: %s\n请先运行 main_runner.R", results_file))
}

if (!file.exists(true_ae_file)) {
  stop(sprintf("找不到验证数据文件: %s", true_ae_file))
}

# 加载结果
cat("加载模型结果...\n")
load(results_file)

# 如果有保存的A-E曲线CSV，直接读取
if (file.exists(ae_csv_file)) {
  cat("读取已计算的A-E曲线...\n")
  ae_df <- read.csv(ae_csv_file)
} else {
  cat("从结果中计算A-E曲线...\n")
  source("03_Result_Analysis/ae_curve_ppd.R")
  
  ae_df <- compute_ae_curve(
    bathy_result = lake_results$bathy_result,
    data_list = lake_results$data_list,
    elevation_step = 0.1
  )
}

# 重新绘制A-E曲线（包含验证数据对比）
cat("\n绘制A-E曲线（包含验证数据对比）...\n")

plot_ae_curve(
  ae_df = ae_df,
  data_list = lake_results$data_list,
  save_path = sprintf("outputs/%s_ae_curve.png", lake_name),
  true_ae_path = true_ae_file
)

cat("\n================================================================================\n")
cat("  ✓ 完成！\n")
cat("  Complete!\n")
cat("================================================================================\n\n")
cat(sprintf("图片已保存: outputs/%s_ae_curve.png\n", lake_name))
cat("请查看图片中的预测曲线（蓝色虚线）和真实曲线（黑色实线）对比\n\n")

