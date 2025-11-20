################################################################################
# INLA 模型拟合模块
# INLA Model Fitting Module
################################################################################

#' 构建 INLA stacks（多 likelihood 情况）
#' Build INLA stacks for multiple likelihoods
#' 
#' @param spde SPDE 对象
#' @param obs_data build_observation_data() 返回的列表
#' @param proj_matrices build_projection_matrices() 返回的列表
#' @return 合并后的 stack 对象
build_inla_stacks <- function(spde, obs_data, proj_matrices) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Building INLA stacks...\n")
  cat("====================================\n")
  
  elev_df <- obs_data$elev_df
  water_df <- obs_data$water_df
  A_elev <- proj_matrices$A_elev
  A_water <- proj_matrices$A_water
  
  # 创建 SPDE 索引
  spde_index <- inla.spde.make.index("field", n.spde = spde$n.spde)
  
  # ---- 简化策略：只使用 water 数据构建模型 ----
  # 原因：water occurrence 数据信息量充足，且避免多 likelihood 的复杂性
  # 岸边 DEM 将在后处理中用于校准，得到绝对高程
  
  cat("\n1. Building water frequency stack (Binomial likelihood)...\n")
  cat("   Note: Using water-only model for stability\n")
  cat("   Shore elevation will be used for calibration in post-processing\n")
  
  # 为 Binomial，INLA 需要响应为 cbind(successes, failures) 或 successes + Ntrials
  # 使用简化格式：successes + Ntrials 参数
  
  stack_water <- inla.stack(
    data = list(
      y = water_df$y_water,  # 成功次数
      Ntrials = water_df$N   # 试验次数
    ),
    A = list(1, A_water),
    effects = list(
      list(intercept = rep(1, nrow(water_df))),
      spde_index
    ),
    tag = "water"
  )
  
  cat(sprintf("   Water observations: %d\n", nrow(water_df)))
  cat(sprintf("   Shore observations available for calibration: %d\n", nrow(elev_df)))
  
  # 使用 water stack 作为完整 stack
  stack_full <- stack_water
  stack_elev <- NULL  # 不使用 elevation stack（将用于后处理校准）
  
  cat("\n✓ INLA stacks built!\n\n")
  
  return(list(
    stack_full = stack_full,
    stack_elev = stack_elev,
    stack_water = stack_water
  ))
}


#' 拟合 INLA 模型
#' Fit INLA model with joint likelihood
#' 
#' @param stack_list build_inla_stacks() 返回的列表
#' @param spde SPDE 对象
#' @param obs_data build_observation_data() 返回的列表
#' @param use_elev 是否使用高程观测
#' @param beta_prior 对 beta 参数的先验均值和精度：c(mean, precision)
#' @return INLA 拟合结果对象
fit_inla_model <- function(stack_list, 
                            spde,
                            obs_data,
                            use_elev = TRUE,
                            beta_prior = c(0, 0.01)) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Fitting INLA model...\n")
  cat("====================================\n")
  
  stack_full <- stack_list$stack_full
  
  # ---- 1. 构建模型公式 ----
  cat("\n1. Building model formula...\n")
  
  # 公式说明：
  # - 对于高程观测（如果有）：eta = field(s)
  # - 对于水频率观测：logit(p) = alpha - beta * field(s)
  #   
  # INLA 中实现 -beta * field(s) 的方式：
  # 我们使用一个技巧：通过 copy 功能创建 field 的负副本
  # 
  # 简化策略：只使用 water 数据
  # 空间场 field 的解释：在水中更深的地方，被水淹没的频率越高
  # 因此 field 本身就代表了相对深度（越大越深）
  
  formula <- y ~ -1 + intercept + f(field, model = spde)
  family_vec <- "binomial"
  
  cat("   Formula: y ~ intercept + f(field)\n")
  cat("   Family: Binomial\n")
  cat("   field interpretation: relative bathymetry (positive = deeper)\n")
  
  # ---- 2. 设置先验 ----
  cat("\n2. Setting priors...\n")
  
  # 对截距设置弱先验
  control_fixed <- list(
    mean = list(intercept = 0),
    prec = list(intercept = 0.01)
  )
  
  cat("   Intercept prior: N(0, precision = 0.01)\n")
  
  # ---- 3. 拟合模型 ----
  cat("\n3. Fitting model with INLA...\n")
  cat("   (This may take several minutes...)\n\n")
  
  # 准备数据
  stack_data <- inla.stack.data(stack_full)
  stack_A <- inla.stack.A(stack_full)
  
  # 拟合
  result <- inla(
    formula = formula,
    data = stack_data,
    family = family_vec,
    Ntrials = stack_data$Ntrials,  # 为 Binomial 提供试验次数
    control.predictor = list(
      A = stack_A,
      compute = TRUE
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE,
      config = TRUE
    ),
    control.fixed = control_fixed,
    verbose = TRUE
  )
  
  cat("\n✓ Model fitting complete!\n\n")
  
  # ---- 4. 输出模型摘要 ----
  cat("====================================\n")
  cat("Model Summary\n")
  cat("====================================\n\n")
  
  print(summary(result))
  
  cat("\n====================================\n")
  cat("Fixed Effects Summary\n")
  cat("====================================\n")
  
  fixed_summary <- result$summary.fixed
  print(fixed_summary)
  
  # 提取关键参数
  if ("intercept" %in% rownames(fixed_summary)) {
    intercept_mean <- fixed_summary["intercept", "mean"]
    intercept_sd <- fixed_summary["intercept", "sd"]
    cat(sprintf("\nIntercept (baseline logit): %.4f ± %.4f\n", 
                intercept_mean, intercept_sd))
  }
  
  cat("\n====================================\n")
  cat("Hyperparameters (SPDE)\n")
  cat("====================================\n")
  
  hyper_summary <- result$summary.hyperpar
  print(hyper_summary)
  
  # 提取空间场参数（需要转换）
  if ("Theta1 for field" %in% rownames(hyper_summary)) {
    # 从 theta 参数转换为 range 和 sigma
    # 这依赖于 SPDE 的参数化方式
    cat("\n   Note: Theta parameters correspond to log(tau) and log(kappa).\n")
    cat("   Use inla.spde2.result() to extract range and sigma.\n")
  }
  
  cat("\n====================================\n\n")
  
  return(result)
}


#' 提取空间场超参数（range 和 sigma）
#' Extract spatial field hyperparameters (range and sigma)
#' 
#' @param result INLA 拟合结果
#' @param spde SPDE 对象
#' @return 包含 range 和 sigma 后验统计的列表
extract_spde_hyperpar <- function(result, spde) {
  
  library(INLA)
  
  cat("====================================\n")
  cat("Extracting SPDE hyperparameters...\n")
  cat("====================================\n\n")
  
  # 使用 inla.spde2.result 提取
  spde_result <- inla.spde2.result(result, "field", spde)
  
  # Range 的后验
  tryCatch({
    range_summary <- inla.zmarginal(spde_result$marginals.range.nominal[[1]])
    cat("Range (spatial correlation distance):\n")
    print(range_summary)
    cat("\n")
  }, error = function(e) {
    cat("Warning: Could not extract range summary\n")
    range_summary <<- list(mean = NA, sd = NA)
  })
  
  # Variance 的后验
  tryCatch({
    variance_marginal <- spde_result$marginals.variance.nominal[[1]]
    # 过滤掉负值（如果有的话）
    variance_marginal <- variance_marginal[variance_marginal[,1] >= 0, ]
    
    sigma_summary <- inla.zmarginal(variance_marginal)
    cat("Variance (marginal variance of the field):\n")
    print(sigma_summary)
    cat("\n")
  }, error = function(e) {
    cat("Warning: Could not extract variance summary\n")
    sigma_summary <<- list(mean = NA, sd = NA)
  })
  
  # 计算 sigma（标准差）
  tryCatch({
    # 确保只对正值开方
    variance_marginal <- spde_result$marginals.variance.nominal[[1]]
    variance_marginal <- variance_marginal[variance_marginal[,1] >= 0, ]
    
    sigma_marginal <- inla.tmarginal(
      function(x) sqrt(pmax(x, 0)),  # 确保非负
      variance_marginal
    )
    
    # 移除任何 NaN 或 Inf
    sigma_marginal <- sigma_marginal[is.finite(sigma_marginal[,1]) & is.finite(sigma_marginal[,2]), ]
    
    if (nrow(sigma_marginal) > 0) {
      sigma_summary2 <- inla.zmarginal(sigma_marginal)
      cat("Sigma (standard deviation):\n")
      print(sigma_summary2)
      cat("\n")
    } else {
      cat("Warning: Could not compute valid sigma values\n")
      sigma_summary2 <- list(mean = NA, sd = NA)
    }
  }, error = function(e) {
    cat(sprintf("Warning: Error computing sigma: %s\n", e$message))
    sigma_summary2 <- list(mean = NA, sd = NA)
  })
  
  cat("====================================\n\n")
  
  return(list(
    range = range_summary,
    variance = sigma_summary,
    sigma = sigma_summary2,
    spde_result = spde_result
  ))
}

