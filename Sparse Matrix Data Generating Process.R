# ======================Sparse Matrix Data Generating Process===================
library(mvtnorm)
library(Matrix)

rm(list=ls())

set.seed(123) 
T <- 20
N <- 10
K <- 10
M <- 10
lambda_T <- 0.25 # tuning parameter
lambda_T1 <- 0.25

# PS:发现一个很纯的问题，之前的协变量，因变量，误差全是随机生成的。。。。
A_star <- matrix(0, nrow = N, ncol = N)
beta_star <- matrix(1,K,1)
num_nonzero <- round(0.05 * N * (N - 1))  # 5%的非对角元素
# 生成所有非对角元素的索引
non_diag_indices <- which(row(A_star) != col(A_star))

selected_indices <- sample(non_diag_indices, num_nonzero)

A_star[selected_indices] <- 0.5

for (i in 1:N) {
  row_sum <- sum(A_star[i, ])
  if (row_sum > 1) {
    A_star[i, ] <- A_star[i, ] / row_sum
  }
}

prespecified_spm_sequence <- list()
for (m in 1:M) {
  W0i <- matrix(0, nrow = N, ncol = N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      binary_relation <- sample(c(0, 1), 1) 
      W0i[i, j] <- binary_relation
      W0i[j, i] <- binary_relation 
    }
  }
  diag(W0i) <- 0
  prespecified_spm_sequence[[m]] <- W0i
}

delta_star <- matrix(0,M,1)
delta_star[c(1,2)] <- 0.2

W_star <- 
  Reduce("+", lapply(1:M, function(r) delta_star[r] * prespecified_spm_sequence[[r]])) + A_star

error_sequence <- list()
# generating variance-covariance matrix of the disturbance term for a single period
generate_cov_matrix <- function(N) {
  cov_matrix <- diag(1, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      if (runif(1) < 0.10) {
        cov_matrix[i, j] <- 0.25
        cov_matrix[j, i] <- 0.25  
      }
    }
  }
  return(cov_matrix)
}
for (t in 1:T) {
  Sigma <- generate_cov_matrix(N)
  e_t <- t(rmvnorm(n = 1, mean = rep(0, N), sigma = Sigma))
  error_sequence[[t]] <- e_t
}


covariate_sequence <- list()
for(t in 1:T){
  X_t <- matrix(rnorm(N*K),nrow = N, ncol = K)
  for (k in 1:K) {
    X_t[,k] <- X_t[,k] + error_sequence[[t]] / 2
  }
  covariate_sequence[[t]] <- X_t
}




dependence_sequence <- list()
for(t in 1:T){
  dependence_sequence[[t]] <- solve(diag(1,N,N) - W_star) %*% covariate_sequence[[t]] %*% beta_star + error_sequence[[t]] 
}

IV_sequence <- list()  #B_t
for (t in 1:T) {
  IV_t <- covariate_sequence[[t]]
  for (m in 1:M) {
    IV_t <- cbind(IV_t,prespecified_spm_sequence[[m]] %*% covariate_sequence[[t]])
  }
  IV_sequence[[t]] <- IV_t
}

IV_average <- matrix(0, nrow = N, ncol = (M+1)*K) #\bar B_t

for (t in 1:T) {
  IV_average <- IV_average + IV_sequence[[t]]
}

IV_average <- IV_average/T #\bar B

IV_demeaned <- list()   # IV_demeaned = B_t - \bar B
for (t in 1:T) {
  IV_demeaned[[t]] <- IV_sequence[[t]] - IV_average
}

Gamma <- matrix(1/((M+1)*K), nrow = (M+1)*K, ncol = 1)

dependence_tilde_sequence <- list()  #\tilde y_i
for (i in 1:N) { # 第四个error
  y_i <- rep(0,N)
  for (t in 1:T) {
    #y_i <- rep(0,N)
    #y_i <- y_i + as.numeric(IV_sequence[[t]][i,] %*% Gamma) * dependence_sequence[[t]]
    y_i <- y_i + as.numeric(IV_demeaned[[t]][i,] %*% Gamma) * dependence_sequence[[t]]
  }
  dependence_tilde_sequence[[i]] <- y_i
}

covariate_tilde_sequence <- list()  #\tilde x_i
for (i in 1:N) { # 第五个error
  x_i <- rep(0,N)
  for (t in 1:T) {
    #x_i <- rep(0,N)
    #x_i <- x_i + as.numeric(IV_sequence[[t]][i,] %*% Gamma) * covariate_sequence[[t]]
    x_i <- x_i + as.numeric(IV_demeaned[[t]][i,] %*% Gamma) * covariate_sequence[[t]]
  }
  covariate_tilde_sequence[[i]] <- x_i
}

# =============================Full Matrix Notation=============================
temp <- IV_demeaned[[1]]
for (t in 2:T) {
  temp <- cbind(temp, IV_demeaned[[t]])
}

B_full <- T ^ (-1/2) * N ^ (-1/2) * kronecker(diag(1,N,N),    #怀疑原文3.1公式有typo
                                              kronecker(diag(1,T,T), t(Gamma)) %*% t(temp))   
B_full <- as(B_full,"sparseMatrix")

y_full <- as.vector(t(do.call(cbind,dependence_sequence)))

epsilon_full <- as.vector(t(do.call(cbind,error_sequence)))

#Z_full <- kronecker(diag(1,N,N), t(do.call(cbind,error_sequence))) 第一个error...
Z_full <- kronecker(diag(1,N,N), t(do.call(cbind,dependence_sequence)))

Z_full <- as(Z_full,"sparseMatrix")

X_beta_func <- function(beta) {
  kronecker(diag(1,N,N),
            kronecker(diag(1,T,T), t(beta)) %*% t(do.call(cbind,covariate_sequence))) 
  
}

X_beta_star <- X_beta_func(beta_star)

dependence_sequence_prim <- lapply(dependence_sequence, t)

y_nu <- t(do.call(cbind, dependence_sequence_prim))

error_sequence_prim <- lapply(error_sequence, t)

epsilon_nu <- t(do.call(cbind, error_sequence_prim))

covariate_sequence_prim <- lapply(covariate_sequence, t)

X_full <- t(do.call(cbind, covariate_sequence_prim))

IV_demeaned_prim <- lapply(IV_demeaned, t)

B_nu <- t(do.call(cbind, IV_demeaned_prim))

temp <- list()
for (t in 1:T) {
  temp[[t]] <- kronecker(covariate_sequence[[t]], IV_demeaned[[t]] %*% Gamma)
}

K_full <- T ^ (-1/2) * N ^ (-1/2) * Reduce('+', temp) %*% solve(t(X_full) %*%
                                                                  B_nu %*% t(B_nu) %*% X_full) %*% t(X_full) %*% B_nu %*% t(B_nu)

W_kronecker <- list()
for (i in 1:M) {
  W_kronecker[[i]] <- as(kronecker(diag(1,T,T), prespecified_spm_sequence[[i]]),"sparseMatrix")
}

H_full <- K_full %*% do.call(cbind,W_kronecker) %*% kronecker(diag(1,M,M), y_nu)

temp <- list()
for (i in 1:M) {
  temp[[i]] <- as.vector(t(prespecified_spm_sequence[[i]]))
}

V_full <- do.call(cbind,temp)

