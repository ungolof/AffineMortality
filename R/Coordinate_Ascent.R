#======================================================================================================
# - Set of negative log-likelihood functions, and functions to implement Coordinate Ascent
# - by groups of parameter vectors. All functions are based on the univ. KD log-likelihood function
#======================================================================================================

#================== - BS independent Coordinate ascent algorithm - =========================
nLL_BSi_uKD_CA <- function(x0, delta, kappa, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_BSi <- function(mu_bar, x0=c(6.960591e-03, 9.017154e-03, 5.091784e-03), delta=c(-1.305830e-06, -5.220474e-02, -1.013210e-01), kappa=c(1.162624e-02, 6.787268e-02, 5.061539e-03), sigma=exp(c(-6.806310, -6.790270, -7.559145)), r=exp(c( -3.327060e+01, -6.086479e-01, -1.553156e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  n_factors <- length(kappa)
  # - Matrix for parameter estimates storage
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"), "log_lik.")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_BSi_KD <- optim(x0_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_BSi_KD$par
    print(paste(sprintf("X(0)_%d", c(1:n_factors)), round(x0_par,3)))

    delta_opt_BSi_KD <- optim(delta_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSi_KD$par
    print(paste(sprintf("delta_%d", c(1:n_factors)), round(delta_par,3)))

    kappa_opt_BSi_KD <- optim(kappa_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSi_KD$par
    print(paste(sprintf("kappa_%d", c(1:n_factors)), round(kappa_par,3)))

    l_sigma_opt_BSi_KD <- optim(l_sigma_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_BSi_KD$par
    print(paste(sprintf("sigma_%d", c(1:n_factors)), round(exp(l_sigma_par),3)))

    l_r_opt_BSi_KD <- optim(l_r_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSi_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="BSi_Table.RData")
    }

    if (abs(nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter)){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], sigma=CA_par[iter_count,((n_factors*3 + 1):(n_factors*4))], r1=CA_par[iter_count,(n_factors*4 + 1)], r2=CA_par[iter_count,(n_factors*4 + 2)], rc=CA_par[iter_count,(n_factors*4 + 3)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], sigma=CA_par[iter_count,((n_factors*3 + 1):(n_factors*4))], r1=CA_par[iter_count,(n_factors*4 + 1)], r2=CA_par[iter_count,(n_factors*4 + 2)], rc=CA_par[iter_count,(n_factors*4 + 3)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}


#================== - BS dependent Coordinate ascent algorithm - =========================

nLL_BSd_2F_uKD_CA <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta)

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- t(B_tT[i,]) %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

nLL_BSd_3F_uKD_CA <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta)

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors)
  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}


# - 2-factors
co_asc_BSd_2F <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02), kappa=c(3.455255e-02, 1.075876e-02), sigma_dg=c(6.789215e-04, 2.036748e-03), Sigma_cov=0, r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), "delta_11", 'delta_21', 'delta_22', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_BSd_KD <- optim(x0_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_BSd_KD$par
    print(paste(sprintf("X(0)_%d", c(1:n_factors)), round(x0_par,3)))

    delta_opt_BSd_KD <- optim(delta_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSd_KD$par
    print(paste(c("delta_11", 'delta_21', 'delta_22'), round(delta_par,3)))

    kappa_opt_BSd_KD <- optim(kappa_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSd_KD$par
    print(paste(sprintf("kappa_%d", c(1:n_factors)), round(kappa_par,3)))

    dg_l_Sigma_chol_opt_BSd_KD <- optim(dg_l_Sigma_chol_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_BSd_KD$par

    odg_Sigma_chol_opt_BSd_KD <- optim(odg_Sigma_chol_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_BSd_KD$par
    print(paste(c('sigma_11', 'sigma_21', 'sigma_22'), round(parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par),3)))

    l_r_opt_BSd_KD <- optim(l_r_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSd_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="BSd2_Table.RData")
    }

    if (abs(nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

  #return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:5], kappa=CA_par[iter_count,c(6:7)], Sigma=list(sigma_11 = CA_par[iter_count,8], sigma_21 = CA_par[iter_count,9], sigma_22 = CA_par[iter_count,10]), r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:5], kappa=CA_par[iter_count,c(6:7)], sigma_dg=CA_par[iter_count,c(8,10)], Sigma_cov=CA_par[iter_count,9], r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

# - Directly provide the parameters

co_asc_BSd_3F <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03, 2.711990e-02), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02, -2.299199e-01, 1.383445e-02, -6.310253e-02), kappa=c(3.455255e-02, 1.075876e-02, 1.000030e-02), sigma_dg=c(6.789215e-04, 2.036748e-03, 1.875928e-03), Sigma_cov=c(-1.260778e-06, 1.194974e-06, -3.718267e-06), r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), "delta_11", 'delta_21', 'delta_22', 'delta_31', 'delta_32', 'delta_33', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', 'sigma_31', 'sigma_32', 'sigma_33', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_BSd_3F_KD <- optim(x0_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_BSd_3F_KD$par
    print(paste(sprintf("X(0)_%d", c(1:n_factors)), round(x0_par,3)))

    delta_opt_BSd_3F_KD <- optim(delta_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSd_3F_KD$par
    print(paste(c("delta_11", 'delta_21', 'delta_22', 'delta_31', 'delta_32', 'delta_33'), round(delta_par,3)))

    kappa_opt_BSd_3F_KD <- optim(kappa_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSd_3F_KD$par
    print(paste(sprintf("kappa_%d", c(1:n_factors)), round(kappa_par,3)))

    dg_l_Sigma_chol_opt_BSd_3F_KD <- optim(dg_l_Sigma_chol_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_BSd_3F_KD$par

    odg_Sigma_chol_opt_BSd_3F_KD <- optim(odg_Sigma_chol_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_BSd_3F_KD$par
    print(paste(c('sigma_11', 'sigma_21', 'sigma_22', 'sigma_31', 'sigma_32', 'sigma_33'), round(parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par),3)))

    l_r_opt_BSd_3F_KD <- optim(l_r_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSd_3F_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="BSd3_Table.RData")
    }

    if (abs(nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

  #return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4:9], kappa=CA_par[iter_count,c(10:12)], Sigma=list(sigma_11 = CA_par[iter_count,13], sigma_21 = CA_par[iter_count,14], sigma_22 = CA_par[iter_count,15], sigma_31 = CA_par[iter_count,16], sigma_32 = CA_par[iter_count,17], sigma_33 = CA_par[iter_count,18]), r1=CA_par[iter_count,19], r2=CA_par[iter_count,20], rc=CA_par[iter_count,21]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4:9], kappa=CA_par[iter_count,c(10:12)], sigma_dg = CA_par[iter_count,c(13,15,18)], Sigma_cov=CA_par[iter_count,c(14,16,17)], r1=CA_par[iter_count,19], r2=CA_par[iter_count,20], rc=CA_par[iter_count,21]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

#============= - Gaussian - Makeham (2 factor) - ==============

## - Independent factors
nLL_GMki_uKD_CA <- function(x0, delta, kappa, l_gamma, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, exp(l_sigma), delta, exp(l_gamma))####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_GMki(age, delta, exp(l_gamma))  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_GMki <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03), delta=c(-5.175933e-02, -5.175921e-02), kappa=c(3.455255e-02, 1.075876e-02), gamma=0.11, sigma=c(6.789215e-04, 2.036748e-03), r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, gamma, sigma, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), 'gamma', sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_gamma_par <- log(gamma)
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_KD <- optim(x0_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_KD$par
    print(paste(c("X(0)_Mak", "X(0)_Gom"), round(x0_par,3)))

    delta_opt_KD <- optim(delta_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_KD$par
    print(paste(sprintf("delta_%d", c(1:n_factors)), round(delta_par,3)))

    kappa_opt_KD <- optim(kappa_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, l_gamma=l_gamma_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_KD$par
    print(paste(sprintf("kappa_%d", c(1:n_factors)), round(kappa_par,3)))

    l_gamma_opt_KD <- optim(l_gamma_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_gamma_par <- l_gamma_opt_KD$par
    print(paste("gamma", round(exp(l_gamma_par),3)))

    l_sigma_opt_KD <- optim(l_sigma_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_KD$par
    print(paste(sprintf("sigma_%d", c(1:n_factors)), round(exp(l_sigma_par),3)))

    l_r_opt_KD <- optim(l_r_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, gamma, sigma, r))] <- c(x0_par, delta_par, kappa_par, exp(l_gamma_par), exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma, r))+1] <- -0.5 * nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="GMki_Table.RData")
    }

    if (abs(nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, gamma, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:4], kappa=CA_par[iter_count,c(5:6)], gamma=CA_par[iter_count,7], sigma=CA_par[iter_count,c(8,9)], r1=CA_par[iter_count,10], r2=CA_par[iter_count,11], rc=CA_par[iter_count,12]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:4], kappa=CA_par[iter_count,c(5:6)], gamma=CA_par[iter_count,7], sigma=CA_par[iter_count,c(8,9)], r1=CA_par[iter_count,10], r2=CA_par[iter_count,11], rc=CA_par[iter_count,12]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}


## - Dependent factors
nLL_GMkd_uKD_CA <- function(x0, delta, kappa, l_gamma, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta)

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors)
  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMkd(age, Low_chol, delta_matrix, exp(l_gamma))####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_GMkd(age, delta_matrix, exp(l_gamma))  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_GMkd <- function(mu_bar, x0=c(2.191140e-03, -8.855686e-03), delta=c(-5.175933e-02, 4.578315e-01, -5.175921e-02), kappa=c(3.455255e-02, 1.075876e-02), gamma=0.11, sigma_dg=c(6.789215e-04, 2.036748e-03), Sigma_cov=-1.260778e-06, r=exp(c(-3.345631e+01, -6.015438e-01, -1.557244e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), "delta_11", 'delta_21', 'delta_22', sprintf("kappa_%d", c(1:n_factors)), 'gamma', 'sigma_11', 'sigma_21', 'sigma_22', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_gamma_par <- log(gamma)
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_KD <- optim(x0_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_KD$par

    delta_opt_KD <- optim(delta_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_gamma=l_gamma_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_KD$par

    kappa_opt_KD <- optim(kappa_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, l_gamma=l_gamma_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_KD$par

    l_gamma_opt_KD <- optim(l_gamma_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_gamma_par <- l_gamma_opt_KD$par

    dg_l_Sigma_chol_opt_KD <- optim(dg_l_Sigma_chol_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_KD$par

    odg_Sigma_chol_opt_KD <- optim(odg_Sigma_chol_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_KD$par

    l_r_opt_KD <- optim(l_r_par, nLL_GMkd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, exp(l_gamma_par), parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r))+1] <- -0.5 * nLL_GMkd_uKD_CA(x0_par, delta_par, kappa_par, gamma_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="GMkd_Table.RData")
    }

    if (abs(nLL_GMkd_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_GMkd_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:5], kappa=CA_par[iter_count,c(6:7)], gamma=CA_par[iter_count,8], Sigma=list(sigma_11 = CA_par[iter_count,9], sigma_21 = CA_par[iter_count,10], sigma_22 = CA_par[iter_count,11]), r1=CA_par[iter_count,12], r2=CA_par[iter_count,13], rc=CA_par[iter_count,14]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:5], kappa=CA_par[iter_count,c(6:7)], gamma=CA_par[iter_count,8], sigma_dg=CA_par[iter_count,c(9,11)], Sigma_cov = CA_par[iter_count,10], r1=CA_par[iter_count,12], r2=CA_par[iter_count,13], rc=CA_par[iter_count,14]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}



#================================== - AFNS independent Coordinate ascent algorithm -=====================================

nLL_AFNSi_uKD_CA <- function(x0, delta, kappa, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFNSi <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, -5.990785e-04), delta=-8.304334e-02, kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03), sigma=exp(c(-7.318991, -7.535594, -8.456025)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S', 'x0_C', "delta", 'kappa_L', 'kappa_S', 'kappa_C', 'sigma_L', 'sigma_S', 'sigma_C', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1

  repeat{

    x0_opt_AFNSi_KD <- optim(x0_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFNSi_KD$par
    print(paste(c('x0_L', 'x0_S', 'x0_C'), round(x0_par,3)))
    suppressWarnings(
    delta_opt_AFNSi_KD <- optim(delta_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000)))
    delta_par <- delta_opt_AFNSi_KD$par
    print(paste("delta", round(delta_par,3)))

    kappa_opt_AFNSi_KD <- optim(kappa_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFNSi_KD$par
    print(paste(c('kappa_L', 'kappa_S', 'kappa_C'), round(kappa_par,3)))

    l_sigma_opt_AFNSi_KD <- optim(l_sigma_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFNSi_KD$par
    print(paste(c('sigma_L', 'sigma_S', 'sigma_C'), round(exp(l_sigma_par),3)))

    l_r_opt_AFNSi_KD <- optim(l_r_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFNSi_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFNSi_Table.RData")
    }

    if (abs(nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,1:3], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,5:7], sigma=CA_par[iter_count,8:10], r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,1:3], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,5:7], sigma=CA_par[iter_count,8:10], r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}

#================== - AFNS dependent Coordinate ascent algorithm - =========================

nLL_AFNSd_uKD_CA <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFNSd <- function(mu_bar, x0=c(9.582516e-03, 1.094110e-02, -1.503155e-03), delta=-7.487697e-02, kappa=c(1.389363e-02, 3.525542e-03, 3.004883e-03), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03), Sigma_cov=c(-8.328978e-06, -3.685028e-06, 3.036376e-06), r=exp(c(-3.335725e+01, -6.066149e-01, -1.552061e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S', 'x0_C', "delta", 'kappa_L', 'kappa_S', 'kappa_C', 'sigma_L', 'sigma_LS', 'sigma_S', 'sigma_LC', 'sigma_SC', 'sigma_C', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_AFNSd_KD <- optim(x0_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFNSd_KD$par
    print(paste(c('x0_L', 'x0_S', 'x0_C'), round(x0_par,3)))

    suppressWarnings(
    delta_opt_AFNSd_KD <- optim(delta_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "BFGS", hessian = TRUE, control=list(maxit = 10000)))
    delta_par <- delta_opt_AFNSd_KD$par
    print(paste("delta", round(delta_par,3)))

    kappa_opt_AFNSd_KD <- optim(kappa_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFNSd_KD$par
    print(paste(c('kappa_L', 'kappa_S', 'kappa_C'), round(kappa_par,3)))

    dg_l_Sigma_chol_opt_AFNSd_KD <- optim(dg_l_Sigma_chol_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_AFNSd_KD$par

    odg_Sigma_chol_opt_AFNSd_KD <- optim(odg_Sigma_chol_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_AFNSd_KD$par
    print(paste(c('sigma_L', 'sigma_LS', 'sigma_S', 'sigma_LC', 'sigma_SC', 'sigma_C'), round(parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par),3)))

    l_r_opt_AFNSd_KD <- optim(l_r_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFNSd_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- - 0.5 * nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFNSd_Table.RData")
    }

    if(abs(nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,c(5:7)], Sigma=list(sigma_L = CA_par[iter_count,8], sigma_LS = CA_par[iter_count,9], sigma_S = CA_par[iter_count,10], sigma_LC = CA_par[iter_count,11], sigma_SC = CA_par[iter_count,12], sigma_C = CA_par[iter_count,13]), r1=CA_par[iter_count,14], r2=CA_par[iter_count,15], rc=CA_par[iter_count,16]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,c(5:7)], sigma_dg=CA_par[iter_count,c(8,10,13)], Sigma_cov=CA_par[iter_count,c(9,11,12)], r1=CA_par[iter_count,14], r2=CA_par[iter_count,15], rc=CA_par[iter_count,16]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

#================================== - AFGNS independent Coordinate ascent algorithm -=====================================

## - Version with imposition of X_L(0) + X_S1(0) + X_S2(0) > 0, and delta1 < delta2 optimized simultaneously
nLL_AFGNSi_uKD_CA <- function(x0, delta, kappa, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- c((exp(x0[1]) - x0[2] - x0[3]), x0[2:5]) #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, exp(l_sigma), delta[1], (delta[1] + exp(delta[2])))
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age, (delta[1] + exp(delta[2])))[2], B_AFNS(age,delta[1])[3], B_AFNS(age,(delta[1] + exp(delta[2])))[3])
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFGNSi <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, 5.990785e-04, 0, 0), delta=c(-8.304334e-02, -0.05), kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03, 0.003, 0.001), sigma=exp(c(-7.318991, -7.535594, -8.456025, -7, -7)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S1', 'x0_S2', 'x0_C1', 'x0_C2', "delta1", 'delta2', 'kappa_L', 'kappa_S1', "kappa_S2", 'kappa_C1', 'kappa_C2', 'sigma_L', 'sigma_S1', 'sigma_S2', 'sigma_C1', 'sigma_C2', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- c(log(sum(x0[1:3])), x0[2:5])
  delta_par <- c(delta[1], log(delta[2] - delta[1]))
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1

  repeat{

    x0_opt_AFGNSi_KD <- optim(x0_par, nLL_AFGNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFGNSi_KD$par
    print(paste(c('x0_L', 'x0_S1', 'x0_S2', 'x0_C1', 'x0_C2'), round(c((exp(x0_par[1]) - x0_par[2] - x0_par[3]), x0_par[2:5]))))

    delta_opt_AFGNSi_KD <- optim(delta_par, nLL_AFGNSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFGNSi_KD$par
    print(paste(c("delta1", 'delta2'), round(c(delta_par[1], (delta_par[1] + exp(delta_par[2]))),3)))

    kappa_opt_AFGNSi_KD <- optim(kappa_par, nLL_AFGNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFGNSi_KD$par
    print(paste(c('kappa_L', 'kappa_S1', "kappa_S2", 'kappa_C1', 'kappa_C2'), round(kappa_par,3)))

    l_sigma_opt_AFGNSi_KD <- optim(l_sigma_par, nLL_AFGNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFGNSi_KD$par
    print(paste(c('sigma_L', 'sigma_S1', 'sigma_S2', 'sigma_C1', 'sigma_C2'), round(exp(l_sigma_par),3)))

    l_r_opt_AFGNSi_KD <- optim(l_r_par, nLL_AFGNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFGNSi_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c((exp(x0_par[1]) - sum(x0_par[c(2:3)])), x0_par[2:5], delta_par[1], (delta_par[1] + exp(delta_par[2])), kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_AFGNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFGNSi_Table.RData")
    }

    if (abs(nLL_AFGNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFGNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,1:5], delta=CA_par[iter_count,6:7], kappa=CA_par[iter_count,8:12], sigma=CA_par[iter_count,13:17], r1=CA_par[iter_count,18], r2=CA_par[iter_count,19], rc=CA_par[iter_count,20]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,1:5], delta=CA_par[iter_count,6:7], kappa=CA_par[iter_count,8:12], sigma=CA_par[iter_count,13:17], r1=CA_par[iter_count,18], r2=CA_par[iter_count,19], rc=CA_par[iter_count,20]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}


#================== - AFGNS dependent Coordinate ascent algorithm - =========================

## - Version with imposition of X_L(0) + X_S1(0) + X_S2(0) > 0, and delta1 < delta2 optimized simultaneously

nLL_AFGNSd_uKD_CA <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol1, odg_Sigma_chol2, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- c((exp(x0[1]) - x0[2] - x0[3]), x0[2], x0[3], x0[4], x0[5]) #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  # - Build diffusion process
  ## - Build lower cholesky factor
  Low_chol <- low_trg_fill_0diag(c(odg_Sigma_chol1, odg_Sigma_chol2))
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], (delta[1] + exp(delta[2])))
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,(delta[1] + exp(delta[2])))[2], B_AFNS(age,delta[1])[3], B_AFNS(age,(delta[1] + exp(delta[2])))[3])
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFGNSd <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, 5.990785e-04, 0, 0), delta=c(-8.304334e-02, -0.05), kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03, 0.003, 0.001), sigma_dg=c(3.215422e-03, 2.625474e-03, 1.164715e-03, 0.0003, 0.0001), Sigma_cov=rep(0, 10), r=exp(c(-3.335725e+01, -6.066149e-01, -1.552061e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S1', 'x0_S2', 'x0_C1', 'x0_C2', "delta1", 'delta2', 'kappa_L', 'kappa_S1', "kappa_S2", 'kappa_C1', 'kappa_C2', 'sigma_L', 'sigma_LS1', 'sigma_S1', 'sigma_LS2', 'sigma_S1S2', 'sigma_S2', 'sigma_LC1', 'sigma_S1C1', 'sigma_S2C1', 'sigma_C1', 'sigma_LC2', 'sigma_S1C2', 'sigma_S2C2', 'sigma_C1C2', 'sigma_C2', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- c(log(sum(x0[1:3])), x0[2:5])
  delta_par <- c(delta[1], log(delta[2] - delta[1]))
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol1_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol[1:5]
  odg_Sigma_chol2_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol[6:10]
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    x0_opt_AFGNSd_KD <- optim(x0_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFGNSd_KD$par
    print(paste(c('x0_L', 'x0_S1', 'x0_S2', 'x0_C1', 'x0_C2'), round(c((exp(x0_par[1]) - x0_par[2] - x0_par[3]), x0_par[2:5]))))

    delta_opt_AFGNSd_KD <- optim(delta_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFGNSd_KD$par
    print(paste(c("delta1", 'delta2'), round(c(delta_par[1], (delta_par[1] + exp(delta_par[2]))),3)))

    kappa_opt_AFGNSd_KD <- optim(kappa_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFGNSd_KD$par
    print(paste(c('kappa_L', 'kappa_S1', "kappa_S2", 'kappa_C1', 'kappa_C2'), round(kappa_par,3)))

    dg_l_Sigma_chol_opt_AFGNSd_KD <- optim(dg_l_Sigma_chol_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_AFGNSd_KD$par

    odg_Sigma_chol1_opt_AFGNSd_KD <- optim(odg_Sigma_chol1_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol1_par <- odg_Sigma_chol1_opt_AFGNSd_KD$par

    odg_Sigma_chol2_opt_AFGNSd_KD <- optim(odg_Sigma_chol2_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, odg_Sigma_chol1=odg_Sigma_chol1_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol2_par <- odg_Sigma_chol2_opt_AFGNSd_KD$par
    print(paste(c('sigma_L', 'sigma_LS1', 'sigma_S1', 'sigma_LS2', 'sigma_S1S2', 'sigma_S2', 'sigma_LC1', 'sigma_S1C1', 'sigma_S2C1', 'sigma_C1', 'sigma_LC2', 'sigma_S1C2', 'sigma_S2C2', 'sigma_C1C2', 'sigma_C2'), round(parest2cov(dg_l_Sigma_chol_par, c(odg_Sigma_chol1_par, odg_Sigma_chol2_par)),3)))

    l_r_opt_AFGNSd_KD <- optim(l_r_par, nLL_AFGNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFGNSd_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c((exp(x0_par[1]) - x0_par[2] - x0_par[3]), x0_par[2:5], delta_par[1], (delta_par[1] + exp(delta_par[2])), kappa_par, parest2cov(dg_l_Sigma_chol_par, c(odg_Sigma_chol1_par, odg_Sigma_chol2_par)), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- - 0.5 * nLL_AFGNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol1_par, odg_Sigma_chol2_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFGNSd_Table.RData")
    }

    if(abs(nLL_AFGNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol1_par, odg_Sigma_chol2_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFGNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol1_par, odg_Sigma_chol2_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

#  return(list(par_est = list(x0=CA_par[iter_count,c(1:5)], delta=CA_par[iter_count,c(6:7)], kappa=CA_par[iter_count,c(8:12)], Sigma=list(sigma_L = CA_par[iter_count,13], sigma_LS1 = CA_par[iter_count,14], sigma_S1 = CA_par[iter_count,15], sigma_LS2 = CA_par[iter_count,16], sigma_S1S2 = CA_par[iter_count,17], sigma_S2 = CA_par[iter_count,18], sigma_LC1 = CA_par[iter_count,19], sigma_S1C1 = CA_par[iter_count,20], sigma_S2C1 = CA_par[iter_count,21], sigma_C1 = CA_par[iter_count,22], sigma_LC2 = CA_par[iter_count,23], sigma_S1C2 = CA_par[iter_count,24], sigma_S2C2 = CA_par[iter_count,25], sigma_C1C2 = CA_par[iter_count,26], sigma_C2 = CA_par[iter_count,27]), r1=CA_par[iter_count,28], r2=CA_par[iter_count,29], rc=CA_par[iter_count,30]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
  return(list(par_est = list(x0=CA_par[iter_count,c(1:5)], delta=CA_par[iter_count,c(6:7)], kappa=CA_par[iter_count,c(8:12)], sigma_dg=CA_par[iter_count,c(13,15,18,22,27)],Sigma_cov=CA_par[iter_count,c(14,16,17,19,20,21,23,24,25,26)], r1=CA_par[iter_count,28], r2=CA_par[iter_count,29], rc=CA_par[iter_count,30]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}


#================================== - AF Unrestricted NS independent Coordinate ascent algorithm -=====================================

nLL_AFUNSi_uKD_CA <- function(x0, delta, kappa, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFUNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFUNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFUNSi <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02, -5.990785e-04), delta=c(-1.305830e-06, -5.220474e-02, -1.013210e-01), kappa=c(9.154603e-03, 1.067658e-02, 7.439715e-03), sigma=exp(c(-7.318991, -7.535594, -8.456025)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S', 'x0_C', "delta_1", "delta_2", "delta_3", 'kappa_L', 'kappa_S', 'kappa_C', 'sigma_L', 'sigma_S', 'sigma_C', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1

  repeat{

    x0_opt_AFUNSi_KD <- optim(x0_par, nLL_AFUNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFUNSi_KD$par
    print(paste(c('x0_L', 'x0_S', 'x0_C'), round(x0_par,3)))

    delta_opt_AFUNSi_KD <- optim(delta_par, nLL_AFUNSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFUNSi_KD$par
    print(paste(c('delta_1', 'delta_2', 'delta_3'), round(delta_par,3)))

    kappa_opt_AFUNSi_KD <- optim(kappa_par, nLL_AFUNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFUNSi_KD$par
    print(paste(c('kappa_L', 'kappa_S', 'kappa_C'), round(kappa_par,3)))

    l_sigma_opt_AFUNSi_KD <- optim(l_sigma_par, nLL_AFUNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFUNSi_KD$par
    print(paste(c('sigma_L', 'sigma_S', 'sigma_C'), round(exp(l_sigma_par),3)))

    l_r_opt_AFUNSi_KD <- optim(l_r_par, nLL_AFUNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFUNSi_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_AFUNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFUNSi_Table.RData")
    }

    if (abs(nLL_AFUNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFUNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:3], delta=CA_par[iter_count,4:6], kappa=CA_par[iter_count,7:9], sigma=CA_par[iter_count,10:12], r1=CA_par[iter_count,13], r2=CA_par[iter_count,14], rc=CA_par[iter_count,15]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}


#================================== - AF Reduced NS independent Coordinate ascent algorithm -=====================================

nLL_AFRNSi_uKD_CA <- function(x0, delta, kappa, l_sigma, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFRNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFRNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_AFRNSi <- function(mu_bar, x0=c(1.091714e-02, 1.002960e-02), delta=-8.304334e-02, kappa=c(9.154603e-03, 1.067658e-02), sigma=exp(c(-7.318991, -7.535594)), r=exp(c(-3.371775e+01, -5.887962e-01, -1.548729e+01)), max_iter=200, tol_lik=0.1, workdir=0){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, r))+1)
  colnames(CA_par) <- c('x0_L', 'x0_S', "delta", 'kappa_L', 'kappa_S', 'sigma_L', 'sigma_S', c("r1", "r2", "rc"), "log_lik")

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  x0_par <- x0
  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 1

  repeat{

    x0_opt_AFRNSi_KD <- optim(x0_par, nLL_AFRNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    x0_par <- x0_opt_AFRNSi_KD$par
    print(paste(c('x0_L', 'x0_S'), round(x0_par,3)))

    delta_opt_AFRNSi_KD <- optim(delta_par, nLL_AFRNSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFRNSi_KD$par
    print(paste("delta", round(delta_par,3)))

    kappa_opt_AFRNSi_KD <- optim(kappa_par, nLL_AFRNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFRNSi_KD$par
    print(paste(c('kappa_L', 'kappa_S'), round(kappa_par,3)))

    l_sigma_opt_AFRNSi_KD <- optim(l_sigma_par, nLL_AFRNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFRNSi_KD$par
    print(paste(c('sigma_L', 'sigma_S'), round(exp(l_sigma_par),3)))

    l_r_opt_AFRNSi_KD <- optim(l_r_par, nLL_AFRNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFRNSi_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_AFRNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="AFRNSi_Table.RData")
    }

    if (abs(nLL_AFRNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFRNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:2], delta=CA_par[iter_count,3], kappa=CA_par[iter_count,4:5], sigma=CA_par[iter_count,6:7], r=CA_par[iter_count,8:10]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}


#======================= - CIR Coordinate ascent algorithm with bounded X - =====================

nLL_CIR_uKD_CA_bd <- function(l_x0, delta, l_kappa, l_sigma, l_theta_P, l_r, mu_bar){

  r_1 <- l_r[1]
  r_2 <- l_r[2]
  r_c <- l_r[3]

  n_factors <- length(l_kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- exp(l_x0)
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, (exp(l_theta_P + l_kappa) / delta), exp(l_sigma), delta)
    B_tT[age,] <- B_CIR(age, exp(l_sigma), delta)
  }

  Phi <- diag(exp(-exp(l_kappa)), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  for(t in 1:n_years){

    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag(exp(2 * l_sigma) * ((1 - exp(-exp(l_kappa))) / exp(l_kappa)) * (0.5 * exp(l_theta_P) * (1 - exp(-exp(l_kappa))) + exp(-exp(l_kappa)) * x_ti[1:n_factors]),n_factors)

    # - First observation
    x_ti <- Phi %*% x_ti + exp(l_theta_P) * (1 - exp(-exp(l_kappa)))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  log_lik <- sum(log(F_ti) + (v_ti^2) / F_ti)
  return(log_lik)
}

co_asc_CIR <- function(mu_bar, x0=c(1.611524e-03, 5.763081e-03, 1.208483e-02), delta=c(-0.12379389, -0.06208546, -0.08131285), kappa=c(4.619791e-02, 3.477558e-01, 4.619791e-02), sigma=c(4.143351e-03, 6.242207e-02, 1.797287e-02), theta_P = c(9.322613e-03, 8.457568e-03, 4.661882e-03), r=c(2.952881e-15, 5.445661e-01, 1.493218e-07), max_iter=200, tol_lik=0.1, workdir=0){

  n_factors <- length(kappa)
  # - Matrix for parameter estimates storage
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, sigma, theta_P, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), sprintf("theta_P_%d", c(1:n_factors)), c("r1", "r2", "rc"), "log_lik")
  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  ## - Factor covariance matrix at t=0

  l_x0_par <- log(x0)
  delta_par <- delta
  l_kappa_par <- log(kappa) # - take logs to ensure positivity
  l_sigma_par <- log(sigma)
  l_theta_P_par <- log(theta_P)
  l_r_par <- log(r)

  iter_count <- 1
  repeat{

    l_x0_opt_CIR_KD <- optim(l_x0_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_x0_par <- l_x0_opt_CIR_KD$par
    print(paste(sprintf("X(0)_%d", c(1:n_factors)), round(exp(l_x0_par),3)))

    delta_opt_CIR_KD <- optim(delta_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, l_x0=l_x0_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    delta_par <- delta_opt_CIR_KD$par
    print(paste(sprintf("delta_%d", c(1:n_factors)), round(delta_par,3)))

    l_kappa_opt_CIR_KD <- optim(l_kappa_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_x0=l_x0_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_kappa_par <- l_kappa_opt_CIR_KD$par
    print(paste(sprintf("kappa_%d", c(1:n_factors)), round(exp(l_kappa_par),3)))

    l_sigma_opt_CIR_KD <- optim(l_sigma_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_x0=l_x0_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_CIR_KD$par
    print(paste(sprintf("sigma_%d", c(1:n_factors)), round(exp(l_sigma_par),3)))

    l_theta_P_opt_CIR_KD <- optim(l_theta_P_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_x0=l_x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_theta_P_par <- l_theta_P_opt_CIR_KD$par
    print(paste(sprintf("theta_P_%d", c(1:n_factors)), round(exp(l_theta_P_par), 3)))

    l_r_opt_CIR_KD <- optim(l_r_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_x0=l_x0_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_CIR_KD$par
    print(paste(c("r1", "r2", "rc"), round(exp(l_r_par),3)))

    # - Store par_est
    CA_par[iter_count, 1:length(c(x0, delta, kappa, sigma, theta_P, r))] <- c(exp(l_x0_par), delta_par, exp(l_kappa_par), exp(l_sigma_par), exp(l_theta_P_par), exp(l_r_par))
    CA_par[iter_count, length(c(x0, delta, kappa, sigma, theta_P, r))+1] <- -0.5 * nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if(workdir != 0){
      setwd(dir = workdir)
      part_l <- CA_par[1:iter_count,]
      save(part_l, file="CIR_Table.RData")
    }

    if (abs(nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar)

      print(paste("log_lik", round(CA_par[iter_count, length(c(x0, delta, kappa, sigma, theta_P, r))+1],2)))
      print(paste("Iteration ", iter_count))
      print(paste("---------------------------------"))

      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], sigma=CA_par[iter_count,((n_factors*3 + 1):(n_factors*4))], theta_P=CA_par[iter_count,((n_factors*4 + 1):(n_factors*5))], r1=CA_par[iter_count,(n_factors*5 + 1)], r2=CA_par[iter_count,(n_factors*5 + 2)], rc=CA_par[iter_count,(n_factors*5 + 3)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, theta_P, r))+1], CA_table = CA_par[1:iter_count,]))
}
