library(nloptr)
library(numDeriv)
library(progress)
library(glmnet)

lasso_stage <- function(lambda_T) {
  r <- 1
  j <- 1
  max_eval <- 10000
  pred_rel <- 1e-2

  res <- list()
  eta_results <- list()
  
  
  initial_x0 <- rep(0,(N*N + M + K))
  delta <- initial_x0[(N*N +1) : (N*N + M)]
  beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
  lasso_x <- t(B_full) %*% Z_full
  lasso_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*% V_full %*% delta - t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))
  
  res[[1]] <- initial_x0
  fit <- glmnet(lasso_x,lasso_y, alpha = 1, lambda = lambda_T, intercept = FALSE)
  res[[1]][1:(N * N)] <- as.vector(coef(fit)[-1, ])
  eta_results[[1]] <- as.vector(coef(fit)[-1, ])
  
  beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
    ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res[[1]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) res[[1]][(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
  
  delta <- solve(t(H_full - t(B_full) %*% Z_full %*% V_full) %*% (H_full - t(B_full) %*% Z_full %*% V_full)) %*%
    t(H_full - t(B_full) %*% Z_full %*% V_full) %*% ( t(B_full) %*% Z_full %*% res[[1]][1: (N*N)] - t(B_full) %*% y_full +
                                                        K_full %*% (diag(1, T*N, T*N) - kronecker(diag(1,T,T), matrix(res[[1]][1: (N*N)],nrow = N, ncol = N, byrow = TRUE))) %*% y_nu )
  
  res[[1]][(N*N+ M + 1):( N*N + M + K)] <- beta
  res[[1]][(N*N +1) : (N*N + M)] <- delta
  
  v <- 2 
  
  for (v in 2:max_eval) {
    
    res[[v]] <- res[[v-1]]
    beta <- res[[v]][(N*N+ M + 1):( N*N + M + K)]
    delta <- res[[v]][(N*N +1) : (N*N + M)]
    
    lasso_x <- t(B_full) %*% Z_full
    lasso_y <- t(B_full) %*% y_full - t(B_full) %*% Z_full %*% V_full %*% delta - t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))
    fit <- glmnet(lasso_x,lasso_y, alpha = 1, lambda = lambda_T, intercept = FALSE)
  
    res[[v]][1:(N * N)] <- as.vector(coef(fit)[-1, ])
    eta_results[[v]] <- as.vector(coef(fit)[-1, ])
    
    beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
      ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res[[v]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) res[[v]][(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
    
    delta <- solve(t(H_full - t(B_full) %*% Z_full %*% V_full) %*% (H_full - t(B_full) %*% Z_full %*% V_full)) %*%
      t(H_full - t(B_full) %*% Z_full %*% V_full) %*% ( t(B_full) %*% Z_full %*% res[[v]][1: (N*N)] - t(B_full) %*% y_full +
                                                          K_full %*% (diag(1, T*N, T*N) - kronecker(diag(1,T,T), matrix(res[[v]][1: (N*N)],nrow = N, ncol = N, byrow = TRUE))) %*% y_nu )
    
    res[[v]][(N*N+ M + 1):( N*N + M + K)] <- beta
    res[[v]][(N*N +1) : (N*N + M)] <- delta
    
    
    if (norm(as.matrix(eta_results[[v]] - eta_results[[v-1]])) < pred_rel  ) {
      print(norm(as.matrix(eta_results[[v]] - eta_results[[v-1]])))
      print("It converges")
      break
    } else {
      print(norm(as.matrix(eta_results[[v]] - eta_results[[v-1]])))
    }
    
  }
  
  ksi_tilde <- matrix(res[[length(res)]][1:(N*N)])
  A_tilde <- matrix(res[[length(res)]][1:(N*N)] ,ncol = N , nrow = N, byrow = TRUE)
  beta_tilde <- matrix(res[[length(res)]][(N*N+ M + 1):( N*N + M + K)])
  delta_tilde <- matrix(res[[length(res)]][(N*N +1) : (N*N + M)])
  lasso_results <- res[[length(res)]]
  
  result_list <- list(
    ksi_tilde = ksi_tilde,
    A_tilde = A_tilde,
    beta_tilde = beta_tilde,
    delta_tilde = delta_tilde,
    lasso_results = lasso_results
  )
  
  return(result_list)
}

result.1 <- lasso_stage(lambda_T = 0.2)


