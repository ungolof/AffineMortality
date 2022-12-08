###############################################################################################
# - Code to obtain the standard errors of the parameter estimates.
# - This implementation simplifies the ordinary one, since it obtains x0 from the
# - Smoothing distribution. In this way, we account for the variability in parameter
# - estimates due to unknown x0, but we avoid to optimize also over x0, this way simplifying
# - the optimization procedure
###############################################################################################

#library(mvtnorm)
#library(TruncatedNormal)
#library(numDeriv)
#library(MASS)


# - Functions to get the standardized residuals

#======================= - Blackburn-Sherris independent factor model - ===================

std_res_BSi <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1,1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i,1] - B_tT[i,] %*% x_ti

    }
  }

  e_hat <- v_ti / sqrt(F_ti)

  return(e_hat)
}

y_star_BSi <- function(x0, delta, kappa, sigma, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(e_star)   # - Number of ages
  n_years <- ncol(e_star)  # - Number of years

  y_star <- matrix(0, nrow(e_star), ncol(e_star))
  F_ti <- matrix(1, nrow(e_star), ncol(e_star))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}

    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

# - The coordinate ascent function is a slight modification of the original one used for estimation
co_asc_BSi_BS <- function(mu_bar, x0, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

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

    #    x0_opt_BSi_KD <- optim(x0_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    #    x0_par <- x0_opt_BSi_KD$par

    delta_opt_BSi_KD <- optim(delta_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSi_KD$par

    kappa_opt_BSi_KD <- optim(kappa_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSi_KD$par

    l_sigma_opt_BSi_KD <- optim(l_sigma_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_BSi_KD$par

    l_r_opt_BSi_KD <- optim(l_r_par, nLL_BSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSi_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if (abs(nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter)){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], sigma=CA_par[iter_count,((n_factors*3 + 1):(n_factors*4))], r1=CA_par[iter_count,(n_factors*4 + 1)], r2=CA_par[iter_count,(n_factors*4 + 2)], rc=CA_par[iter_count,(n_factors*4 + 3)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}

## - Boostrap estimation of the standard errors (Covariance estimation)

CovEst_BS_BSi <- function(x0, delta, kappa, sigma, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma, r)))
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  colnames(par_table) <- c(sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"))

  res_table <- matrix(NA, n_ages, n_years)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_BSi(x0, delta, kappa, sigma, r, mu_bar)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_BSi(x0_s, delta, kappa, sigma, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_BSi_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1:n_factors] <- par_est_BS_i$par_est$delta
    par_table[i,(n_factors+1):(n_factors*2)] <- par_est_BS_i$par_est$kappa
    par_table[i,(n_factors*2+1):(n_factors*3)] <- par_est_BS_i$par_est$sigma
    par_table[i,(n_factors*3+1):(n_factors*4)] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    #par_table_BSi[i,] <- optim(c(delta, kappa, log(sigma), log(r)), nLL_BSi_uKD_BS, mu_bar_BS=mu_bar_BS, gr = NULL, method = "BFGS", hessian = FALSE, control=list(maxit = 10000))$par # or coordinate ascent
    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "BSiBS_Table.RData")
    }

  }

  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}


#============================== - Blackburn-Sherris dependent factor model - ===================

std_res_BSd_2F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- t(B_tT[i,]) %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  e_hat <- v_ti / sqrt(F_ti)

  return(e_hat)
}

y_star_BSd_2F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta)

  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

co_asc_BSd_2F_BS <- function(mu_bar, x0, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

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

    delta_opt_BSd_KD <- optim(delta_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSd_KD$par

    kappa_opt_BSd_KD <- optim(kappa_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSd_KD$par

    dg_l_Sigma_chol_opt_BSd_KD <- optim(dg_l_Sigma_chol_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_BSd_KD$par

    odg_Sigma_chol_opt_BSd_KD <- optim(odg_Sigma_chol_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_BSd_KD$par

    l_r_opt_BSd_KD <- optim(l_r_par, nLL_BSd_2F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSd_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if (abs(nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSd_2F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,c(1:2)], delta=CA_par[iter_count,3:5], kappa=CA_par[iter_count,c(6:7)], Sigma=list(sigma_11 = CA_par[iter_count,8], sigma_21 = CA_par[iter_count,9], sigma_22 = CA_par[iter_count,10]), r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

CovEst_BS_BSd_2F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  colnames(par_table) <- c("delta_11", 'delta_21', 'delta_22', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', c("r1", "r2", "rc"))
  res_table <- matrix(NA, n_ages, n_years)
  n_factors <- length(kappa)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar=mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_BSd_2F(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - for starting values
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    # - 4) Parameter estimation

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star (small correction wrt the initial code)
    mu_bar_BS <- y_star_BSd_2F(x0_s, delta, kappa, sigma_dg, Sigma_cov, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_BSd_2F_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1:3] <- par_est_BS_i$par_est$delta
    par_table[i,4:5] <- par_est_BS_i$par_est$kappa
    par_table[i,6:8] <- c(par_est_BS_i$par_est$Sigma$sigma_11, par_est_BS_i$par_est$Sigma$sigma_21, par_est_BS_i$par_est$Sigma$sigma_22)
    par_table[i,9:11] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "BSd2BS_Table.RData")
    }

    }

  # 5) Get st. err. parameter estimates
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}


std_res_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- t(B_tT[i,]) %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  e_hat <- v_ti / sqrt(F_ti)

  return(e_hat)
}

y_star_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta)

  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

co_asc_BSd_3F_BS <- function(mu_bar, x0, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

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

    delta_opt_BSd_3F_KD <- optim(delta_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSd_3F_KD$par

    kappa_opt_BSd_3F_KD <- optim(kappa_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSd_3F_KD$par

    dg_l_Sigma_chol_opt_BSd_3F_KD <- optim(dg_l_Sigma_chol_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_BSd_3F_KD$par

    odg_Sigma_chol_opt_BSd_3F_KD <- optim(odg_Sigma_chol_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_BSd_3F_KD$par

    l_r_opt_BSd_3F_KD <- optim(l_r_par, nLL_BSd_3F_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSd_3F_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- -0.5 * nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if (abs(nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_BSd_3F_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4:9], kappa=CA_par[iter_count,c(10:12)], Sigma=list(sigma_11 = CA_par[iter_count,13], sigma_21 = CA_par[iter_count,14], sigma_22 = CA_par[iter_count,15], sigma_31 = CA_par[iter_count,16], sigma_32 = CA_par[iter_count,17], sigma_33 = CA_par[iter_count,18]), r1=CA_par[iter_count,19], r2=CA_par[iter_count,20], rc=CA_par[iter_count,21]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

CovEst_BS_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  colnames(par_table) <- c("delta_11", 'delta_21', 'delta_22', 'delta_31', 'delta_32', 'delta_33', sprintf("kappa_%d", c(1:n_factors)), 'sigma_11', 'sigma_21', 'sigma_22', 'sigma_31', 'sigma_32', 'sigma_33', c("r1", "r2", "rc"))
  res_table <- matrix(NA, n_ages, n_years)
  n_factors <- length(kappa)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar=mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_BSd_3F(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - for starting values
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    # - 4) Parameter estimation

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star (small correction wrt the initial code)
    mu_bar_BS <- y_star_BSd_3F(x0_s, delta, kappa, sigma_dg, Sigma_cov, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_BSd_3F_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1:6] <- par_est_BS_i$par_est$delta
    par_table[i,7:9] <- par_est_BS_i$par_est$kappa
    par_table[i,10:15] <- c(par_est_BS_i$par_est$Sigma$sigma_11, par_est_BS_i$par_est$Sigma$sigma_21, par_est_BS_i$par_est$Sigma$sigma_22, par_est_BS_i$par_est$Sigma$sigma_31, par_est_BS_i$par_est$Sigma$sigma_32, par_est_BS_i$par_est$Sigma$sigma_33)
    par_table[i,16:18] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "BSd3BS_Table.RData")
    }
    }

  # 5) Get st. err. parameter estimates
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}

#================================= - Gompertz - Makeham - ===================

std_res_GMki <- function(x0, delta, kappa, gamma, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
    A_tT[age,1] <- A_GMki(age, sigma, delta, gamma)
    B_tT[age,] <- B_GMki(age,delta, gamma)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1,1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i,1] - B_tT[i,] %*% x_ti

    }
  }

  e_hat <- v_ti / sqrt(F_ti)

  return(e_hat)
}

y_star_GMki <- function(x0, delta, kappa, gamma, sigma, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(e_star)   # - Number of ages
  n_years <- ncol(e_star)  # - Number of years

  y_star <- matrix(0, nrow(e_star), ncol(e_star))
  F_ti <- matrix(1, nrow(e_star), ncol(e_star))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, sigma, delta, gamma)
    B_tT[age,] <- B_GMki(age,delta, gamma)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}

    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

# - The coordinate ascent function is a slight modification of the original one used for estimation
co_asc_GMki_BS <- function(mu_bar, x0, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

  n_factors <- length(kappa)
  # - Matrix for parameter estimates storage
  CA_par <- matrix(NA, nrow=max_iter, ncol=length(c(x0, delta, kappa, gamma, sigma, r))+1)
  colnames(CA_par) <- c(sprintf("x0_%d", c(1:n_factors)), sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), "gamma", sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"), "log_lik.")

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

    #    x0_opt_GMki_KD <- optim(x0_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    #    x0_par <- x0_opt_GMki_KD$par

    delta_opt_GMki_KD <- optim(delta_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_GMki_KD$par

    kappa_opt_GMki_KD <- optim(kappa_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_GMki_KD$par

    l_gamma_opt_GMki_KD <- optim(l_gamma_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, l_sigma=l_sigma_par, kappa=kappa_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_gamma_par <- l_gamma_opt_GMki_KD$par

    l_sigma_opt_GMki_KD <- optim(l_sigma_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_GMki_KD$par

    l_r_opt_GMki_KD <- optim(l_r_par, nLL_GMki_uKD_CA, mu_bar=mu_bar, x0=x0_par, delta=delta_par, kappa=kappa_par, l_gamma=l_gamma_par, l_sigma=l_sigma_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_GMki_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_gamma_par), exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if (abs(nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter)){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_GMki_uKD_CA(x0_par, delta_par, kappa_par, l_gamma_par, l_sigma_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], gamma=CA_par[iter_count,((n_factors*3) + 1)], sigma=CA_par[iter_count,((n_factors*3 + 2):(n_factors*4+1))], r1=CA_par[iter_count,(n_factors*4 + 2)], r2=CA_par[iter_count,(n_factors*4 + 3)], rc=CA_par[iter_count,(n_factors*4 + 4)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}

## - Boostrap estimation of the standard errors (Covariance estimation)

CovEst_BS_GMki <- function(x0, delta, kappa, gamma, sigma, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma, r)))
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  colnames(par_table) <- c(sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), "gamma", sprintf("sigma_%d", c(1:n_factors)), c("r1", "r2", "rc"))

  res_table <- matrix(NA, n_ages, n_years)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_GMki_uKD(x0, delta, kappa, gamma, sigma, r, mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_GMki(x0, delta, kappa, gamma, sigma, r, mu_bar)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_GMki(x0_s, delta, kappa, gamma, sigma, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_GMki_BS(mu_bar_BS, x0=x0_s, delta, kappa, gamma, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1:n_factors] <- par_est_BS_i$par_est$delta
    par_table[i,(n_factors+1):(n_factors*2)] <- par_est_BS_i$par_est$kappa
    par_table[i,(n_factors*2)+1] <- par_est_BS_i$par_est$gamma
    par_table[i,(n_factors*2+2):(n_factors*3+1)] <- par_est_BS_i$par_est$sigma
    par_table[i,(n_factors*3+2):(n_factors*4+1)] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "GMkiBS_Table.RData")
    }

  }

  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}

#============================ - AFNS model with independent factors - ===================

std_res_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

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

  e_hat <- v_ti / sqrt(F_ti)

  return(e_hat)
}

y_star_AFNSi <- function(x0, delta, kappa, sigma, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

co_asc_AFNSi_BS <- function(mu_bar, x0, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

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

    delta_opt_AFNSi_KD <- optim(delta_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFNSi_KD$par

    kappa_opt_AFNSi_KD <- optim(kappa_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, l_sigma=l_sigma_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFNSi_KD$par

    l_sigma_opt_AFNSi_KD <- optim(l_sigma_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFNSi_KD$par

    l_r_opt_AFNSi_KD <- optim(l_r_par, nLL_AFNSi_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFNSi_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))] <- c(x0_par, delta_par, kappa_par, exp(l_sigma_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par))+1] <- - 0.5 * nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)

    if (abs(nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFNSi_uKD_CA(x0_par, delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:3], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,5:7], sigma=CA_par[iter_count,8:10], r1=CA_par[iter_count,11], r2=CA_par[iter_count,12], rc=CA_par[iter_count,13]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, r))+1], CA_table = CA_par[1:iter_count,]))
}

CovEst_BS_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma, r)))
  colnames(par_table) <- c("delta", sprintf("kappa_%s", c("L", 'S', 'C')), sprintf("sigma_%s", c("L", 'S', 'C')), c("r1", "r2", "rc"))
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  res_table <- matrix(NA, n_ages, n_years)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar=mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_AFNSi(x0, delta, kappa, sigma, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_AFNSi(x0_s, delta, kappa, sigma, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_AFNSi_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1] <- par_est_BS_i$par_est$delta
    par_table[i,2:4] <- par_est_BS_i$par_est$kappa
    par_table[i,5:7] <- par_est_BS_i$par_est$sigma
    par_table[i,8:10] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "AFNSiBS_Table.RData")
    }

  }

  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}


#========================= - AFNS model with dependent factors - ===================

std_res_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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

  e_hat <- v_ti / sqrt(F_ti)
  return(e_hat)
}

y_star_AFNSd <- function(x0, delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  y_star <- matrix(0, nrow(mu_bar), ncol(mu_bar))
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
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

co_asc_AFNSd_BS <- function(mu_bar, x0, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

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

    delta_opt_AFNSd_KD <- optim(delta_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, x0=x0_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "BFGS", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFNSd_KD$par

    kappa_opt_AFNSd_KD <- optim(kappa_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, x0=x0_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFNSd_KD$par

    dg_l_Sigma_chol_opt_AFNSd_KD <- optim(dg_l_Sigma_chol_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, x0=x0_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_AFNSd_KD$par

    odg_Sigma_chol_opt_AFNSd_KD <- optim(odg_Sigma_chol_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, x0=x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_AFNSd_KD$par

    l_r_opt_AFNSd_KD <- optim(l_r_par, nLL_AFNSd_uKD_CA, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, x0=x0_par, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFNSd_KD$par

    # - Store par_est
    CA_par[iter_count,1:length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))] <- c(x0_par, delta_par, kappa_par, parest2cov(dg_l_Sigma_chol_par, odg_Sigma_chol_par), exp(l_r_par))
    CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1] <- - 0.5 * nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)


    if(abs(nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_AFNSd_uKD_CA(x0_par, delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,c(1:3)], delta=CA_par[iter_count,4], kappa=CA_par[iter_count,c(5:7)], Sigma=list(sigma_L = CA_par[iter_count,8], sigma_LS = CA_par[iter_count,9], sigma_S = CA_par[iter_count,10], sigma_LC = CA_par[iter_count,11], sigma_SC = CA_par[iter_count,12], sigma_C = CA_par[iter_count,13]), r1=CA_par[iter_count,14], r2=CA_par[iter_count,15], rc=CA_par[iter_count,16]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma_dg, Sigma_cov, r))+1], CA_table = CA_par[1:iter_count,]))
}

CovEst_BS_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  colnames(par_table) <- c("delta", sprintf("kappa_%d", c(1:n_factors)), 'sigma_L', 'sigma_LS', 'sigma_S', 'sigma_LC', 'sigma_SC', 'sigma_C', c("r1", "r2", "rc"))
  res_table <- matrix(NA, n_ages, n_years)
  n_factors <- length(kappa)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar=mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:3]

  # - 1) Get standardized residuals
  std_r <- std_res_AFNSd(x0, delta, kappa, sigma_dg, Sigma_cov, r)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - for starting values
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    # - 4) Parameter estimation

    ## - 4.3) Sample x0 from the smoothing distribution
    x0_s <- as.numeric(mvtnorm::rmvnorm(n = 1, mean=X_t_sm, sigma = S_t_sm))

    # - 3) Get mu_bar_star (small correction wrt the initial code)
    mu_bar_BS <- y_star_AFNSd(x0_s, delta, kappa, sigma_dg, Sigma_cov, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_AFNSd_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1] <- par_est_BS_i$par_est$delta
    par_table[i,2:4] <- par_est_BS_i$par_est$kappa
    par_table[i,5:10] <- c(par_est_BS_i$par_est$Sigma$sigma_11, par_est_BS_i$par_est$Sigma$sigma_21, par_est_BS_i$par_est$Sigma$sigma_22, par_est_BS_i$par_est$Sigma$sigma_31, par_est_BS_i$par_est$Sigma$sigma_32, par_est_BS_i$par_est$Sigma$sigma_33)
    par_table[i,11:13] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "AFNSdBS_Table.RData")
    }

    }

  # 5) Get st. err. parameter estimates
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}


#============================= - CIR - ===================

std_res_CIR <- function(x0, delta, kappa, sigma, theta_P, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- matrix(0, nrow(mu_bar), ncol(mu_bar))
  F_ti <- matrix(1, nrow(mu_bar), ncol(mu_bar))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance

  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- x0
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, (theta_P * kappa / delta), sigma, delta)
    B_tT[age,] <- B_CIR(age, sigma, delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  for(t in 1:n_years){

    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * x_ti[1:n_factors]),n_factors)

    # - First observation
    x_ti <- Phi %*% x_ti + theta_P * (1 - exp(-kappa))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
  }

  e_hat <- v_ti / sqrt(F_ti)
  return(e_hat)
}

y_star_CIR <- function(x0, delta, kappa, sigma, theta_P, r, e_star){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(e_star)   # - Number of ages
  n_years <- ncol(e_star)  # - Number of years

  y_star <- matrix(0, nrow(e_star), ncol(e_star))
  F_ti <- matrix(1, nrow(e_star), ncol(e_star))

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance

  # - Initialize X and Sigma
  #x_ti <- x_0 #init_X
  x_ti <- x0
  P_ti <- 1e-10 * diag(1, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, (theta_P * kappa / delta), sigma, delta)
    B_tT[age,] <- B_CIR(age, sigma, delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  for(t in 1:n_years){

    R[,(t * n_factors + 1):((t+1) * n_factors)] <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * x_ti[1:n_factors]),n_factors)

    # - First observation
    x_ti <- Phi %*% x_ti + theta_P * (1 - exp(-kappa))
    x_ti <- l_bound(x_ti)
    P_ti <- Phi %*% P_ti %*% t(Phi) + R[,(t * n_factors + 1):((t+1) * n_factors)]     # - P_{1,t}
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))
    y_star[1,t] <- A_tT[1,1] + B_tT[1,] %*% x_ti + sqrt(F_ti[1,t]) * e_star[1,t]

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / sqrt(F_ti[i-1,t])) %*% e_star[i-1,t]
      x_ti <- l_bound(x_ti)

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      y_star[i,t] <- A_tT[i,1] + B_tT[i,] %*% x_ti + sqrt(F_ti[i,t]) * e_star[i,t]

    }
  }

  return(y_star)
}

co_asc_CIR_BS <- function(mu_bar, x0, delta, kappa, sigma, theta_P, r, max_iter=200, tol_lik=0.1){

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

    delta_opt_CIR_KD <- optim(delta_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, l_x0=l_x0_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    delta_par <- delta_opt_CIR_KD$par

    l_kappa_opt_CIR_KD <- optim(l_kappa_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_x0=l_x0_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_kappa_par <- l_kappa_opt_CIR_KD$par

    l_sigma_opt_CIR_KD <- optim(l_sigma_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_x0=l_x0_par, l_theta_P=l_theta_P_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_CIR_KD$par

    l_theta_P_opt_CIR_KD <- optim(l_theta_P_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_x0=l_x0_par, l_r=l_r_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_theta_P_par <- l_theta_P_opt_CIR_KD$par

    l_r_opt_CIR_KD <- optim(l_r_par, nLL_CIR_uKD_CA_bd, mu_bar=mu_bar, delta=delta_par, l_kappa=l_kappa_par, l_sigma=l_sigma_par, l_theta_P=l_theta_P_par, l_x0=l_x0_par, gr = NULL, method = "Nelder-Mead", hessian = FALSE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_CIR_KD$par

    # - Store par_est
    CA_par[iter_count, 1:length(c(x0, delta, kappa, sigma, theta_P, r))] <- c(exp(l_x0_par), delta_par, exp(l_kappa_par), exp(l_sigma_par), exp(l_theta_P_par), exp(l_r_par))
    CA_par[iter_count, length(c(x0, delta, kappa, sigma, theta_P, r))+1] <- -0.5 * nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar) - 0.5 * nrow(mu_bar) * ncol(mu_bar)


    if (abs(nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar) - neg_loglikelihood) < tol_lik | (iter_count==max_iter) ){
      break
    }else{
      # - Update log-likelihood
      neg_loglikelihood <- nLL_CIR_uKD_CA_bd(l_x0_par, delta_par, l_kappa_par, l_sigma_par, l_theta_P_par, l_r_par, mu_bar)
      iter_count <- iter_count + 1
    }
  }

  return(list(par_est = list(x0=CA_par[iter_count,1:n_factors], delta=CA_par[iter_count,((n_factors + 1):(n_factors*2))], kappa=CA_par[iter_count,((n_factors*2 + 1):(n_factors*3))], sigma=CA_par[iter_count,((n_factors*3 + 1):(n_factors*4))], theta_P=CA_par[iter_count,((n_factors*4 + 1):(n_factors*5))], r1=CA_par[iter_count,(n_factors*5 + 1)], r2=CA_par[iter_count,(n_factors*5 + 2)], rc=CA_par[iter_count,(n_factors*5 + 3)]), log_lik = CA_par[iter_count,length(c(x0, delta, kappa, sigma, theta_P, r))+1], CA_table = CA_par[1:iter_count,]))
}

CovEst_BS_CIR <- function(x0, delta, kappa, sigma, theta_P, r, mu_bar, n_BS=500, t_ex = 4, max_it=200, tolerance_lev=0.1, workdir=0){
  par_table <- matrix(NA, n_BS, length(c(delta, kappa, sigma, theta_P, r)))
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  colnames(par_table) <- c(sprintf("delta_%d", c(1:n_factors)), sprintf("kappa_%d", c(1:n_factors)), sprintf("sigma_%d", c(1:n_factors)), sprintf("theta_P_%d", c(1:n_factors)), c("r1", "r2", "rc"))

  res_table <- matrix(NA, n_ages, n_years)

  # - 4) Parameter estimation
  ## - 4.1) Get filtered estimates
  Filtering <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_P, r, mu_bar)
  X_t_fil <- Filtering$X_t
  X_t_c_fil <- Filtering$X_t_c
  S_t_fil <- Filtering$S_t
  S_t_c_fil <- Filtering$S_t_c

  ## - 4.2) Get smoothed estimate of x0
  Smoothing <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)
  X_t_sm <- Smoothing$X_t_sm[,1]
  S_t_sm <- Smoothing$S_t_sm[,1:n_factors]

  # - 1) Get standardized residuals
  std_r <- std_res_CIR(x0, delta, kappa, sigma, theta_P, r, mu_bar)
  ## - Fill the first t_ex residuals, which stay the same
  res_table[,1:t_ex] <- std_r[,1:t_ex]

  # - Repeat bootstrapping samples
  for(i in 1:n_BS){

    # - 2) resample standardized residuals
    t_sample <- sample(c((t_ex+1):n_years), (n_years-t_ex), replace=TRUE)  # - starting from 5 to account for KF start up irregularity
    res_table[,(t_ex + 1):n_years] <- std_r[,t_sample]

    ## - 4.3) Sample x0 from the smoothing distribution - adjusted to ensure a lower bound of zero
    x0_s <- as.numeric(TruncatedNormal::rtmvnorm(n=1, mu=X_t_sm, sigma=S_t_sm, lb=rep(1e-10, n_factors)))  #l_bound(rmvnorm(n = 1, mu=X_t_sm, Sigma = S_t_sm)))

    # - 3) Get mu_bar_star
    mu_bar_BS <- y_star_CIR(x0_s, delta, kappa, sigma, theta_P, r, res_table)

    ## - 4.4) Get bootstrapped parameter estimates
    par_est_BS_i <- co_asc_CIR_BS(mu_bar_BS, x0=x0_s, delta, kappa, sigma, theta_P, r, max_iter=max_it, tol_lik=tolerance_lev)
    par_table[i,1:n_factors] <- par_est_BS_i$par_est$delta
    par_table[i,(n_factors+1):(n_factors*2)] <- par_est_BS_i$par_est$kappa
    par_table[i,(n_factors*2+1):(n_factors*3)] <- par_est_BS_i$par_est$sigma
    par_table[i,(n_factors*3+1):(n_factors*4)] <- par_est_BS_i$par_est$theta_P
    par_table[i,(n_factors*4+1):(n_factors*4 + 3)] <- c(par_est_BS_i$par_est$r1, par_est_BS_i$par_est$r2, par_est_BS_i$par_est$rc)

    if(workdir != 0){
      setwd(dir = workdir)
      save(par_table[1:i,], "CIRBS_Table.RData")
    }

    print(paste(i, "% Bootstrap estimate"))
  }

  # 5) Get st. err. parameter estimates (unconstrained parameters)
  cov_pe <- cov(par_table)
  serr_pe <- sqrt(diag(cov_pe))

  return(list(Par_Table=par_table, Cov = cov_pe, St.err = serr_pe))
}




