#GLMNET_AdaptiveLassoStage
#========================== The Block Coordinate Descent Algorihrm (Adaptive LASSO Stage) ========================
adaptivelasso_stage <- function(lambda_T,lambda_T1,result.1) {
  
  r <- 1
  j <- 1
  max_eval <- 10000
  pred_rel <- 1e-2
  
  
  ksi_tilde = result.1$ksi_tilde
  A_tilde = result.1$A_tilde
  beta_tilde = result.1$beta_tilde
  delta_tilde = result.1$delta_tilde
  lasso_results = result.1$lasso_results
  # adaptive lasso stage notation
  res.ad <- list()
  eta_results.ad <- list()
  
  v_temp <- ksi_tilde
  v_temp[v_temp == 0] <- 1e-8
  u_temp <- delta_tilde
  u_temp[u_temp == 0] <- 1e-8
  
  u_full <- as.matrix(sapply(sapply(u_temp,abs),function(element) element^-1),1,M)
  v_full <- as.matrix(sapply(sapply(v_temp,abs),function(element) element^-1),1,N*N)
  
  
  # 实际上在写关于g_function的函数，ksi是自变量
  # 这段代码不运行，eq3.8
  #h_hat <- T^(-1/2) * N^(-1/2) * temp %*% solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*%
  #                                   t(X_full) %*% B_nu %*% t(B_nu) %*% (diag(1,T*N,T*N) - kronecker(diag(1,T,T),matrix(A_hat,N,N))) %*% y_nu  
  
  
  # 以下代码正常运行
  # 计算g^{(r-1)} #第二个typo
  g_function <- function(ksi) {
    temp <- kronecker(covariate_sequence[[1]], (IV_sequence[[1]] - IV_average) %*% Gamma)
    for (t in 2:T) {
      temp <- temp + kronecker(covariate_sequence[[t]], (IV_sequence[[t]] - IV_average) %*% Gamma)
    }
    g_r1 <- T^(-1/2) * N^(-1/2) * temp %*% solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full)  %*%
      t(X_full) %*% B_nu %*% t(B_nu) %*% (diag(1,T*N,T*N) - kronecker(diag(1,T,T),matrix(ksi,N,N,byrow = TRUE))) %*% y_nu
    return(g_r1)
  }
  
  res.ad[[1]] <- lasso_results
  
  delta <- lasso_results[(N*N +1) : (N*N + M)]
  beta <- lasso_results[(N*N+ M + 1):( N*N + M + K)]
  
  g_full <- g_function(lasso_results[1:(N*N)])
  
  lasso_x <- t(B_full) %*% Z_full
  lasso_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*% V_full %*% delta - t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))
  fit <- glmnet(lasso_x,lasso_y, alpha = 1, lambda = lambda_T, intercept = FALSE, penalty.factor = v_full)
  
  res.ad[[1]][1:(N * N)] <- as.vector(coef(fit)[-1, ])
  eta_results.ad[[1]] <- as.vector(coef(fit)[-1, ])
  
  beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
    ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res.ad[[1]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) res.ad[[1]][(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
  
  lasso_delta_x <- t(B_full) %*% Z_full %*% V_full - H_full
  lasso_delta_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*%  res.ad[[1]][1 : (N*N)] - g_full
  fit <- glmnet(lasso_delta_x,lasso_delta_y, alpha = 1, lambda = lambda_T1, intercept = FALSE, penalty.factor = u_full)
  delta <- as.vector(coef(fit)[-1,])
  
 
  
  res.ad[[1]][(N*N+ M + 1):( N*N + M + K)] <- beta
  res.ad[[1]][(N*N +1) : (N*N + M)] <- delta
  
  v <- 2 
  
  for (v in 2:max_eval) {
    res.ad[[v]] <- res.ad[[v-1]]
    
    delta <- res.ad[[v]][(N*N +1) : (N*N + M)]
    beta <- res.ad[[v]][(N*N+ M + 1):( N*N + M + K)]
    
    g_full <- g_function(res.ad[[v]][1:(N*N)])
    
    lasso_x <- t(B_full) %*% Z_full
    lasso_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*% V_full %*% delta - t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))
    fit <- glmnet(lasso_x,lasso_y, alpha = 1, lambda = lambda_T, intercept = FALSE, penalty.factor = v_full)
    
    res.ad[[v]][1:(N * N)] <- as.vector(coef(fit)[-1, ])
    eta_results.ad[[v]] <- as.vector(coef(fit)[-1, ])
    
    beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
      ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res.ad[[v]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) res.ad[[v]][(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
    
    lasso_delta_x <- t(B_full) %*% Z_full %*% V_full - H_full
    lasso_delta_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*% res.ad[[v]][1 : (N*N)] - g_full
    fit <- glmnet(lasso_delta_x,lasso_delta_y, alpha = 1, lambda = lambda_T1, intercept = FALSE, penalty.factor = u_full)
    delta <- as.vector(coef(fit)[-1,])
    
    
    res.ad[[v]][(N*N+ M + 1):( N*N + M + K)] <- beta
    res.ad[[v]][(N*N +1) : (N*N + M)] <- delta
    
    if (norm(as.matrix(eta_results.ad[[v]] - eta_results.ad[[v-1]])) < pred_rel  ) {
      print(norm(as.matrix(eta_results.ad[[v]] - eta_results.ad[[v-1]])))
      print("It converges")
      break
    } else {
      print(norm(as.matrix(eta_results.ad[[v]] - eta_results.ad[[v-1]])))
    }
  }
  
  ksi_hat <- matrix(res.ad[[length(res.ad)]][1:(N*N)])
  A_hat <- matrix(res.ad[[length(res.ad)]][1:(N*N)] ,ncol = N , nrow = N, byrow = TRUE)
  beta_hat <- matrix(res.ad[[length(res.ad)]][(N*N+ M + 1):( N*N + M + K)])
  delta_hat <- matrix(res.ad[[length(res.ad)]][(N*N +1) : (N*N + M)])
  lasso_results.ad <- res.ad[[length(res.ad)]]
  
  result_list <- list(
    ksi_hat = ksi_hat,
    A_hat = A_hat,
    beta_hat = beta_hat,
    delta_hat = delta_hat,
    lasso_results_ad = lasso_results.ad
  )
  
  return(result_list)
}

result.2 <- adaptivelasso_stage(lambda_T = 0.01,lambda_T1 = 0.01,result.1 = result.1)



