#============================== The Block Coordinate Descent Algotihrm (LASSO Stage) ========================
library(nloptr)
library(numDeriv)
library(progress)
r <- 1
j <- 1
max_eval <- 1000
pred_rel <- 1e-4

opts <- list(
  algorithm = "NLOPT_LD_SLSQP", 
  maxeval = 2000,               
  xtol_rel = 1e-1
)


res <- list()
optimization_results <- list()
eta_results <- list()

# =========================第一个循环，优化A第一行元素==========================
initial_x0 <- rep(0,(N*N + M + K))

constraint_functions <- function(x) {
  eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0[(N*N +1) : (N*N + M)]
  beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
  eta[1:N] <- x
  temp <- matrix(eta,N,N, byrow = TRUE)
  for (r in 1:M) {
    temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
  }
  constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
  return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
}

objective_functions <- function(x) {
  eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0[(N*N +1) : (N*N + M)]
  beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
  eta[1:N] <- x
  ob <- (2*T) ^ -1 * norm(
    t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
      t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * sum(abs(eta)) 
  return(ob)
}

objective_functions_grad <- function(x) {
  return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
}

constraint_functions_grad <- function(x) {
  return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
}

x0 <- initial_x0[1:N] # 优化起点

optimization_results[[1]] <- nloptr(
  x0 = x0,                                     # 初始值
  eval_f = objective_functions,           # 目标函数
  eval_grad_f = objective_functions_grad, # 目标函数梯度
  eval_g_ineq = constraint_functions,     # 不等式约束
  eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
  opts = opts                                  # 选项
)

# =========================第一个循环，优化A剩余行的元素========================
for (i in 2:N) {
  constraint_functions <- function(x) {
    eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0[(N*N +1) : (N*N + M)]
    beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
    eta[((i-1) * N + 1):(i * N)] <- x
    temp <- matrix(eta,N,N, byrow = TRUE)
    for (r in 1:M) {
      temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
    }
    constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
    return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
  }
  
  objective_functions <- function(x) {
    eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0[(N*N +1) : (N*N + M)]
    beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
    eta[((i-1) * N + 1):(i * N)] <- x
    ob <- norm(
      t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
        t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * sum(abs(eta))
    return(ob)
  }
  
  objective_functions_grad <- function(x) {
    return(grad(objective_functions,x))
  }
  
  constraint_functions_grad <- function(x) {
    return(jacobian(constraint_functions, x))
  }
  
  x0 <- initial_x0[((i-1) * N + 1):(i * N)]
  
  optimization_results[[i]] <- nloptr(
    x0 = x0,                                     # 初始值
    eval_f = objective_functions,           # 目标函数
    eval_grad_f = objective_functions_grad, # 目标函数梯度
    eval_g_ineq = constraint_functions,     # 不等式约束
    eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
    opts = opts                                  # 选项
  )
  
}
# 将第一次循环的优化结果储存到res[[1]]中
# 这边发现第三个error,beta计算里的temp没修正
# 现在怀疑原文又有一定的typo
res[[1]] <- initial_x0
for (i in 1:N) {
  res[[1]][((i-1) * N + 1):(i * N)] <- optimization_results[[i]]$solution
}

beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
  ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res[[1]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) initial_x0[(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu

delta <- solve(t(H_full - t(B_full) %*% Z_full %*% V_full) %*% (H_full - t(B_full) %*% Z_full %*% V_full)) %*%
  t(H_full - t(B_full) %*% Z_full %*% V_full) %*% ( t(B_full) %*% Z_full %*% res[[1]][1: (N*N)] - t(B_full) %*% y_full +
                                                      K_full %*% (diag(1, T*N, T*N) - kronecker(diag(1,T,T), matrix(res[[1]][1: (N*N)],nrow = N, ncol = N, byrow = TRUE))) %*% y_nu )

res[[1]][(N * N + M + 1):(N * N + K + M)] <- beta
res[[1]][(N * N + 1):(N * N + M)] <- delta
eta_results[[1]] <- res[[1]][1:(N*N)]


# ==========================完整的迭代收敛过程==================================
v <- 2

pb <- progress_bar$new(
  format = "  处理 [:bar] :percent 完成，耗时:elapsed",
  total = max_eval -1 , clear = FALSE, width = 60
)

for (v in 2:max_eval) {
  pb$tick()    #进度条显示，可删除
  initial_x0 <- res[[v-1]] #继承上次循环得到的结果
  
  constraint_functions <- function(x) {
    eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0[(N*N +1) : (N*N + M)]
    beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
    eta[1:N] <- x
    temp <- matrix(eta,N,N, byrow = TRUE)
    for (r in 1:M) {
      temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
    }
    constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
    return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
  }
  
  objective_functions <- function(x) {
    eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0[(N*N +1) : (N*N + M)]
    beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
    eta[1:N] <- x
    ob <- norm(
      t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
        t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * sum(abs(eta)) 
    return(ob)
  }
  
  objective_functions_grad <- function(x) {
    return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
  }
  
  constraint_functions_grad <- function(x) {
    return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
  }
  
  x0 <- initial_x0[1:N] # 优化起点
  
  optimization_results[[1]] <- nloptr(
    x0 = x0,                                     # 初始值
    eval_f = objective_functions,           # 目标函数
    eval_grad_f = objective_functions_grad, # 目标函数梯度
    eval_g_ineq = constraint_functions,     # 不等式约束
    eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
    opts = opts                                  # 选项
  )
  
  for (i in 2:N) {
    
    constraint_functions <- function(x) {
      eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
      delta <- initial_x0[(N*N +1) : (N*N + M)]
      beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
      eta[((i-1) * N + 1):(i * N)] <- x
      temp <- matrix(eta,N,N, byrow = TRUE)
      for (r in 1:M) {
        temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
      }
      constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
      return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
    }
    
    objective_functions <- function(x) {
      eta <- initial_x0[(1): (N*N)]  #储存了A的各元素，按行排列
      delta <- initial_x0[(N*N +1) : (N*N + M)]
      beta <- initial_x0[(N*N+ M + 1):( N*N + M + K)]
      eta[((i-1) * N + 1):(i * N)] <- x
      ob <- norm(
        t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
          t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * sum(abs(eta))
      return(ob)
    }
    
    objective_functions_grad <- function(x) {
      return(grad(objective_functions,x))
    }
    
    constraint_functions_grad <- function(x) {
      return(jacobian(constraint_functions, x))
    }
    
    x0 <- initial_x0[((i-1) * N + 1):(i * N)]
    
    optimization_results[[i]] <- nloptr(
      x0 = x0,                                     # 初始值
      eval_f = objective_functions,           # 目标函数
      eval_grad_f = objective_functions_grad, # 目标函数梯度
      eval_g_ineq = constraint_functions,     # 不等式约束
      eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
      opts = opts                                  # 选项
    )
    
  }
  
  res[[v]] <- initial_x0
  for (i in 1:N) {
    res[[v]][((i-1) * N + 1):(i * N)] <- optimization_results[[i]]$solution
  }
  
  beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
    ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res[[v]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) initial_x0[(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
  
  delta <- solve(t(H_full - t(B_full) %*% Z_full %*% V_full) %*% (H_full - t(B_full) %*% Z_full %*% V_full)) %*%
    t(H_full - t(B_full) %*% Z_full %*% V_full) %*% ( t(B_full) %*% Z_full %*% res[[v]][1: (N*N)] - t(B_full) %*% y_full +
                                                        K_full %*% (diag(1, T*N, T*N) - kronecker(diag(1,T,T), matrix(res[[v]][1: (N*N)],nrow = N, ncol = N, byrow = TRUE))) %*% y_nu )
  
  res[[v]][(N * N + M + 1):(N * N + K + M)] <- beta
  res[[v]][(N * N + 1):(N * N + M)] <- delta
  eta_results[[v]] <- res[[v]][1:(N*N)]
  
  
  if (norm(as.matrix(eta_results[[v]] - eta_results[[v-1]])) < pred_rel  ) {
    print("It converges")
    break
  } else {
    print(norm(as.matrix(eta_results[[v]] - eta_results[[v-1]])))
  }
}

mm <- rep(0,500)
for (v in 2:500) {
  mm[v]<- norm(as.matrix(eta_results[[v]] - eta_results[[v-1]]))
}

ksi_tilde <- matrix(res[[length(res)]][1:(N*N)])
A_tilde <- matrix(res[[length(res)]][1:(N*N)] ,ncol = N , nrow = N, byrow = TRUE)
beta_tilde <- matrix(res[[length(res)]][(N*N+ M + 1):( N*N + N + M)])
delta_tilde <- matrix(res[[length(res)]][(N*N +1) : (N*N + M)])
lasso_results <- res[[length(res)]]

#一个新idea，考虑N + 1个内置的constraint来进行优化
#上边的idea可行，最终可以收敛





