# =======================Choice of tuning parameter=============================
S_hat <- 5 #The indices for the non-zero values of ksi_hat(the solution for the adaptive lasso)
BIC <- log(T ^ -2 * norm(
  t(B_full) %*% y_full - t(B_full) %*% Z_full %*% ksi_hat - t(B_full) %*% Z_full %*% V_full %*% delta -
    t(B_full) %*% X_beta_func(beta_hat) %*% as.vector(diag(1,N,N))) ^ 2) +
  abs(S_hat) * log(T)/T * log(log(2*N - 2))


