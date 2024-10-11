# 定义计算 BIC 的函数
estswm <- function(lambda_T, lambda_T1) {
  result.1 <- lasso_stage(lambda_T = lambda_T)
  result.2 <- adaptivelasso_stage(lambda_T = lambda_T,lambda_T1 = lambda_T1,result.1 = result.1)
  
  ksi_hat = result.2$ksi_hat
  A_hat = result.2$A_hat
  beta_hat = result.2$beta_hat
  delta_hat = result.2$delta_hat
  lasso_results_ad = result.2$lasso_results.ad
  
  S_hat <- sum(ksi_hat != 0) 
  BIC <- log(T ^ (-2) * norm(
    t(B_full) %*% y_full - t(B_full) %*% Z_full %*% ksi_hat - t(B_full) %*% Z_full %*% V_full %*% delta_hat -
      t(B_full) %*% X_beta_func(beta_hat) %*% as.vector(diag(1,N,N))) ^ 2) +
    abs(S_hat) * log(T)/T * log(log(2*N - 2))
  return(BIC)
}

# 创建参数组合
lambda_T_values <- 10^seq(log10(0.01), log10(0.02), length.out = 6)
lambda_T1_values <- 10^seq(0, log10(0.005), length.out = 6)
# 创建一个包含所有参数组合的网格
tuning_grid <- expand.grid(lambda_T = lambda_T_values, lambda_T1 = lambda_T1_values, stringsAsFactors = FALSE)

# 定义一个函数来计算BIC
evaluate_bic <- function(params) {
  lambda_T <- as.numeric(params["lambda_T"])
  lambda_T1 <- as.numeric(params["lambda_T1"])
  BIC <- estswm(lambda_T = lambda_T, lambda_T1 = lambda_T1)
  return(BIC)
}

# 计算每个参数组合的 BIC
bic_values <- sapply(1:nrow(tuning_grid), function(i) {
  evaluate_bic(tuning_grid[i, ])
})

# 找到最小的BIC以及对应的参数组合
min_bic_index <- which.min(bic_values)
optimal_params <- tuning_grid[min_bic_index, ]

print(optimal_params)
print(paste("Minimum BIC:", bic_values[min_bic_index]))
