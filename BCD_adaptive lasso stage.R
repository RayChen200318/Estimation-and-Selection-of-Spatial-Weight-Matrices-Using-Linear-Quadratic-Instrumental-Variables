#========================== The Block Coordinate Descent Algorihrm (Adaptive LASSO Stage) ========================
r <- 1
j <- 1
max_eval <- 1000
pred_rel <- 1e-4

opts <- list(
  algorithm = "NLOPT_LD_SLSQP", 
  maxeval = 1000,               
  xtol_rel = 1e-1
)


# adaptive lasso stage notation
res.ad <- list()
eta_results.ad <- list()
optimization_results.delta <- list()
optimization_results.eta <- list()

u_full <- as.matrix(sapply(sapply(delta_tilde,abs),function(element) element^-1),1,M)
v_full <- as.matrix(sapply(sapply(as.vector(t(A_tilde)),abs),function(element) element^-1),1,N*N)


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

# =========================优化3'步，A的第一行元素==============================
initial_x0.ad <- lasso_results

constraint_functions <- function(x) {
  eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
  beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
  eta[1:N] <- x
  temp <- matrix(eta,N,N, byrow = TRUE)
  for (r in 1:M) {
    temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
  }
  constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
  return(constraint)
}

objective_functions <- function(x) {
  eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
  beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
  eta[1:N] <- x
  ob <- norm(
    t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
      t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * t(v_full) %*% abs(eta) 
  return(ob)
}

objective_functions_grad <- function(x) {
  return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
}

constraint_functions_grad <- function(x) {
  return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
}

x0 <- initial_x0.ad[1:N] # 优化起点

optimization_results.eta[[1]] <- nloptr(
  x0 = x0,                                     # 初始值
  eval_f = objective_functions,           # 目标函数
  eval_grad_f = objective_functions_grad, # 目标函数梯度
  eval_g_ineq = constraint_functions,     # 不等式约束
  eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
  opts = opts                                  # 选项
)

# =========================优化3'步，A的剩余行元素==============================
for (i in 2:N) {
  constraint_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    eta[((i-1) * N + 1):(i * N)] <- x
    temp <- matrix(eta,N,N, byrow = TRUE)
    for (r in 1:M) {
      temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
    }
    constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
    return(constraint)
  }
  
  objective_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    eta[((i-1) * N + 1):(i * N)] <- x
    ob <- norm(
      t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
        t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * t(v_full) %*% abs(eta) 
    return(ob)
  }
  
  objective_functions_grad <- function(x) {
    return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
  }
  
  constraint_functions_grad <- function(x) {
    return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
  }
  
  x0 <- initial_x0.ad[((i-1) * N + 1):(i * N)] # 优化起点
  
  optimization_results.eta[[i]] <- nloptr(
    x0 = x0,                                     # 初始值
    eval_f = objective_functions,           # 目标函数
    eval_grad_f = objective_functions_grad, # 目标函数梯度
    eval_g_ineq = constraint_functions,     # 不等式约束
    eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
    opts = opts                                  # 选项
  )
  
}
# 将第一次循环的优化结果储存到res[[1]]中
res.ad[[1]] <- initial_x0.ad
for (i in 1:N) {
  res.ad[[1]][((i-1) * N + 1):(i * N)] <- optimization_results.eta[[i]]$solution
}
eta_results.ad[[1]] <- res[[1]][1:(N * N)]

# ============================优化2',delta的元素================================
constraint_functions <- function(x) {
  eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
  beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
  delta <- x
  temp <- matrix(eta,N,N, byrow = TRUE)
  for (r in 1:M) {
    temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
    #print(temp)
  }
  constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
  return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
}

objective_functions <- function(x) {
  eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
  delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
  beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
  delta <- x
  g_full <- g_function(eta)
  ob <- norm(
    t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - g_full - (t(B_full) %*% Z_full %*% V_full - H_full) %*% delta) ^ 2 + 
    lambda_T1 * t(u_full) %*% abs(delta)
  return(ob)
}

objective_functions_grad <- function(x) {
  return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
}

constraint_functions_grad <- function(x) {
  return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
}

x0 <- initial_x0.ad[(N * N + 1):(N * N + M)] # 优化起点

optimization_results.delta <- nloptr(
  x0 = x0,                                     # 初始值
  eval_f = objective_functions,           # 目标函数
  eval_grad_f = objective_functions_grad, # 目标函数梯度
  eval_g_ineq = constraint_functions,     # 不等式约束
  eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
  opts = opts                                  # 选项
)

beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
  ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res.ad[[1]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) initial_x0[(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu

res.ad[[1]][(N * N + M + 1):(N * N + K + M)] <- beta
res.ad[[1]][(N * N + 1):(N * N + M)] <- optimization_results.delta$solution

# ==============================完整迭代过程====================================
v <- 2

pb <- progress_bar$new(
  format = "  处理 [:bar] :percent 完成，耗时:elapsed",
  total = max_eval -1 , clear = FALSE, width = 60
)

for (v in 2:max_eval) {
  initial_x0.ad <- res[[v-1]]
  
  pb$tick()
  # =========================优化3'步，A的第一行元素==============================
  constraint_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    eta[1:N] <- x
    temp <- matrix(eta,N,N, byrow = TRUE)
    for (r in 1:M) {
      temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
    }
    constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
    return(constraint)
  }
  
  objective_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    eta[1:N] <- x
    ob <- norm(
      t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
        t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * t(v_full) %*% abs(eta) 
    return(ob)
  }
  
  objective_functions_grad <- function(x) {
    return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
  }
  
  constraint_functions_grad <- function(x) {
    return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
  }
  
  x0 <- initial_x0.ad[1:N] # 优化起点
  
  optimization_results.eta[[1]] <- nloptr(
    x0 = x0,                                     # 初始值
    eval_f = objective_functions,           # 目标函数
    eval_grad_f = objective_functions_grad, # 目标函数梯度
    eval_g_ineq = constraint_functions,     # 不等式约束
    eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
    opts = opts                                  # 选项
  )
  
  # =========================优化3'步，A的剩余行元素==============================
  for (i in 2:N) {
    constraint_functions <- function(x) {
      eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
      delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
      beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
      eta[((i-1) * N + 1):(i * N)] <- x
      temp <- matrix(eta,N,N, byrow = TRUE)
      for (r in 1:M) {
        temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
      }
      constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
      return(constraint)
    }
    
    objective_functions <- function(x) {
      eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
      delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
      beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
      eta[((i-1) * N + 1):(i * N)] <- x
      ob <- norm(
        t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - t(B_full) %*% Z_full %*% V_full %*% delta - 
          t(B_full) %*% X_beta_func(beta) %*% as.vector(diag(1,N,N))) ^ 2 + lambda_T * t(v_full) %*% abs(eta) 
      return(ob)
    }
    
    objective_functions_grad <- function(x) {
      return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
    }
    
    constraint_functions_grad <- function(x) {
      return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
    }
    
    x0 <- initial_x0.ad[((i-1) * N + 1):(i * N)] # 优化起点
    
    optimization_results.eta[[i]] <- nloptr(
      x0 = x0,                                     # 初始值
      eval_f = objective_functions,           # 目标函数
      eval_grad_f = objective_functions_grad, # 目标函数梯度
      eval_g_ineq = constraint_functions,     # 不等式约束
      eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
      opts = opts                                  # 选项
    )
    
  }
  # 将第一次循环的优化结果储存到res[[1]]中
  res.ad[[v]] <- initial_x0.ad
  for (i in 1:N) {
    res.ad[[v]][((i-1) * N + 1):(i * N)] <- optimization_results.eta[[i]]$solution
  }
  eta_results.ad[[v]] <- res[[v]][1:(N * N)]
  
  # ============================优化2',delta的元素================================
  constraint_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    delta <- x
    temp <- matrix(eta,N,N, byrow = TRUE)
    for (r in 1:M) {
      temp <- temp + delta[r] * prespecified_spm_sequence[[r]]
      #print(temp)
    }
    constraint <- abs(temp) %*% matrix(1,N,1) - matrix(1,N,1)
    return(rbind(constraint, abs(t(delta) %*% matrix(1,M,1)) - 1))
  }
  
  objective_functions <- function(x) {
    eta <- initial_x0.ad[(1): (N*N)]  #储存了A的各元素，按行排列
    delta <- initial_x0.ad[(N*N +1) : (N*N + M)]
    beta <- initial_x0.ad[(N*N+ M + 1):( N*N + K + M)]
    delta <- x
    g_full <- g_function(eta)
    ob <- norm(
      t(B_full) %*% y_full - t(B_full) %*% Z_full %*% eta - g_full - (t(B_full) %*% Z_full %*% V_full - H_full) %*% delta) ^ 2 + 
      lambda_T1 * t(u_full) %*% abs(delta)
    return(ob)
  }
  
  objective_functions_grad <- function(x) {
    return(grad(objective_functions,x,method.args = list(eps = 1e-4)))
  }
  
  constraint_functions_grad <- function(x) {
    return(jacobian(constraint_functions, x,method.args = list(eps = 1e-4)))
  }
  
  x0 <- initial_x0.ad[(N * N + 1):(N * N + M)] # 优化起点
  
  optimization_results.delta <- nloptr(
    x0 = x0,                                     # 初始值
    eval_f = objective_functions,           # 目标函数
    eval_grad_f = objective_functions_grad, # 目标函数梯度
    eval_g_ineq = constraint_functions,     # 不等式约束
    eval_jac_g_ineq = constraint_functions_grad, # 约束的梯度
    opts = opts                                  # 选项
  )
  
  beta <- solve(t(X_full) %*% B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu) %*%
    ( diag(1,T*N,T*N) - kronecker(diag(1,T,T), matrix(res.ad[[v]][1 : (N*N)],nrow = N, ncol = N, byrow = TRUE)) - Reduce("+", lapply(1:M, function(r) initial_x0[(N * N + 1):(N * N + M)][r] * W_kronecker[[r]]))  ) %*% y_nu
  
  res.ad[[v]][(N * N + M + 1):(N * N + K + M)] <- beta
  res.ad[[v]][(N * N + 1):(N * N + M)] <- optimization_results.delta$solution

  if (norm(as.matrix(eta_results.ad[[v]] - eta_results.ad[[v-1]])) < pred_rel  ) {
    print("It converges")
    break
  } else {
    print(norm(as.matrix(eta_results.ad[[v]] - eta_results.ad[[v-1]])))
  }

}

ksi_hat <- matrix(res.ad[[length(res.ad)]][1:(N*N)])
A_tilde <- matrix(res[[length(res)]][1:(N*N)] ,ncol = N , nrow = N, byrow = TRUE)
beta_tilde <- matrix(res[[length(res)]][(N*N+ M + 1):( N*N + N + M)])
delta_tilde <- matrix(res[[length(res)]][(N*N +1) : (N*N + M)])
lasso_results <- res[[length(res)]]

  