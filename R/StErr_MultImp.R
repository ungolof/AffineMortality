#======= Estimation of standard errors using multiple imputation

#library(mvtnorm)
#library(TruncatedNormal)
#library(numDeriv)
#library(MASS)

# ============ - Blackburn-Sherris independent model (version without x0)

# - Step 0: Build log-likelihood

nLL_BSi_serr <- function(delta, kappa, l_sigma, l_r, mu_bar, X_rnd_nLL){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_BSi_serr_H <- function(vdParameters, mu_bar, X_rnd_H){
  n_factors <- nrow(X_rnd_H)

  delta = vdParameters[1:n_factors]
  kappa = vdParameters[(n_factors+1):(n_factors*2)]
  l_sigma = vdParameters[(n_factors*2+1):(n_factors*3)]
  l_r = vdParameters[(n_factors*3+1):(n_factors*4)]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_BSi_MI <- function(mu_bar, X_random, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 0

  repeat{
    neg_loglikelihood <- nLL_BSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_BSi_KD <- optim(delta_par, nLL_BSi_serr, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSi_KD$par

    kappa_opt_BSi_KD <- optim(kappa_par, nLL_BSi_serr, delta=delta_par, l_sigma=l_sigma_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSi_KD$par

    l_sigma_opt_BSi_KD <- optim(l_sigma_par, nLL_BSi_serr, delta=delta_par, kappa=kappa_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_BSi_KD$par

    l_r_opt_BSi_KD <- optim(l_r_par, nLL_BSi_serr, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSi_KD$par

    neg_loglikelihood_condition <- nLL_BSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random)

    if((abs(neg_loglikelihood_condition - neg_loglikelihood) < tol_lik) | (iter_count==max_iter)){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par))
}


# - Function for the calculation of the covariance matrix and of the standard errors
# - It takes the parameter values and the number of draws as input (default value set equal to 50)
CovEst_MI_BSi <- function(x0, delta, kappa, sigma, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  X_t_fil <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t
  X_t_c_fil <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t_c
  S_t_fil <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t
  S_t_c_fil <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma, r)), length(c(delta, kappa, sigma, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_BSi_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_BSi_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  nabla_grad <- diag(c(rep(1, n_factors*2), sigma, r))

  Cov_est_orig <- t(nabla_grad) %*% Cov_est %*% nabla_grad

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  return(list(Cov_par_est = Cov_est_orig, St_err=list(delta=serr_est_orig[1:n_factors], kappa=serr_est_orig[(n_factors+1):(n_factors*2)], sigma=serr_est_orig[(2*n_factors+1):(n_factors*3)], r1=serr_est_orig[(3*n_factors+1)], r2=serr_est_orig[(3*n_factors+2)], rc=serr_est_orig[(3*n_factors+3)])))

}

#============= - Blackburn-Sherris dependent 3 factors - ==========================

# - Version with variance
parest2cov_jac <- function(dg_l_Sigma_chol_odg_Sigma_chol){

  n_factors <- (-1 + sqrt(1 + 8 * length(dg_l_Sigma_chol_odg_Sigma_chol))) * 0.5
  dg_l_Sigma_chol <- dg_l_Sigma_chol_odg_Sigma_chol[1:n_factors]
  odg_Sigma_chol <- dg_l_Sigma_chol_odg_Sigma_chol[(n_factors+1):length(dg_l_Sigma_chol_odg_Sigma_chol)]

  Low_Chol <- matrix(0, n_factors, n_factors)
  diag(Low_Chol) <- exp(dg_l_Sigma_chol)
  Low_Chol <- Low_Chol + low_trg_fill_0diag(odg_Sigma_chol)

  Sigma_diffusion <- matrix(NA, n_factors, n_factors)
  Sigma_diffusion <- Low_Chol %*% t(Low_Chol)
  #  diag(Sigma_diffusion) <- sqrt(diag(Sigma_diffusion))

  Sigma_el <- rep(0, (length(dg_l_Sigma_chol_odg_Sigma_chol)))
  Sigma_el[1:n_factors] <- diag(Sigma_diffusion)
  count_vec <- n_factors + 1
  for(row in 2:n_factors){
    for(col in 1:(row-1)){
      Sigma_el[count_vec] <- Sigma_diffusion[row,col]
      count_vec <- count_vec + 1
    }
  }

  return(Sigma_el)
}

nLL_BSd_3F_serr <- function(delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar, X_rnd_nLL){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_BSd_3F_serr_H <- function(vdParameters, mu_bar, X_rnd_H){

  n_factors <- 3
  # - Parameters
  delta <- vdParameters[1:6]
  kappa <- vdParameters[7:9] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[10:12] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[13:15] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[18]
  r_1 <- vdParameters[16]
  r_2 <- vdParameters[17]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_BSd_3F_MI <- function(mu_bar, X_random, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 0
  repeat{
    neg_loglikelihood <- nLL_BSd_3F_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_BSd_3F_KD <- optim(delta_par, nLL_BSd_3F_serr, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSd_3F_KD$par

    kappa_opt_BSd_3F_KD <- optim(kappa_par, nLL_BSd_3F_serr, delta=delta_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSd_3F_KD$par

    dg_l_Sigma_chol_opt_BSd_3F_KD <- optim(dg_l_Sigma_chol_par, nLL_BSd_3F_serr, delta=delta_par, kappa=kappa_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_BSd_3F_KD$par

    odg_Sigma_chol_opt_BSd_3F_KD <- optim(odg_Sigma_chol_par, nLL_BSd_3F_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_BSd_3F_KD$par

    l_r_opt_BSd_3F_KD <- optim(l_r_par, nLL_BSd_3F_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSd_3F_KD$par

    if (abs(nLL_BSd_3F_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar, X_rnd_nLL = X_random) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par))
}

CovEst_MI_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  pe_dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  pe_odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  X_t_fil <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t
  X_t_c_fil <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t_c
  S_t_fil <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t
  S_t_c_fil <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r)), length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_BSd_3F_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_BSd_3F_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  ## - Build sub-jacobian for Sigma (use created function parest2cov in Est_fun_v0)

  subjac <- numDeriv::jacobian(parest2cov_jac, c(pe_dg_l_Sigma_chol, pe_odg_Sigma_chol), method="Richardson", side=NULL)

  Jac <- matrix(0, 18, 18)
  diag(Jac[1:9, 1:9]) <- rep(1, 9) # the untransformed delta
  Jac[10:15, 10:15] <- subjac
  diag(Jac[16:18, 16:18]) <- r

  Cov_est_orig <- t(Jac) %*% Cov_est %*% Jac

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  ## - Second Jacobian step to reconcile sigma_dg^2 to sigma_dg
  nabla_grad <- matrix(0, 18, 18)
  diag(nabla_grad[c(1:9, 13:18), c(1:9, 13:18)]) <- 1
  diag(nabla_grad[c(10:12), c(10:12)]) <- 0.5 / sigma_dg

  Cov_est_orig2 <- t(nabla_grad) %*% Cov_est_orig %*% nabla_grad

  serr_est_orig2 <- sqrt(diag(Cov_est_orig2))

  return(list(Cov_par_est = Cov_est_orig2, St_err=list(delta=serr_est_orig2[1:6], kappa=serr_est_orig2[7:9], sigma_dg=serr_est_orig2[10:12], Sigma_cov=serr_est_orig2[13:15], r1=serr_est_orig2[16], r2=serr_est_orig2[17], rc=serr_est_orig2[18])))

}


# ============ - AFNS independent model - =====================

nLL_AFNSi_serr <- function(delta, kappa, l_sigma, l_r, mu_bar, X_rnd_nLL){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- rep(0, n_ages)
  for(age in 1:n_ages){
    # - Original diag from Blackburn-Sherris (2013)
    H[age] <- exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:age))) / age
  }
  # H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(H)) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / H) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_AFNSi_serr_H <- function(vdParameters, mu_bar, X_rnd_H){

  n_factors <- 3

  delta = vdParameters[1]
  kappa = vdParameters[2:4]
  l_sigma = vdParameters[5:7]
  l_r = vdParameters[8:10]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_AFNSi_MI <- function(mu_bar, X_random, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 0
  repeat{
    neg_loglikelihood <- nLL_AFNSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_BSi_KD <- optim(delta_par, nLL_AFNSi_serr, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_BSi_KD$par

    kappa_opt_BSi_KD <- optim(kappa_par, nLL_AFNSi_serr, delta=delta_par, l_sigma=l_sigma_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_BSi_KD$par

    l_sigma_opt_BSi_KD <- optim(l_sigma_par, nLL_AFNSi_serr, delta=delta_par, kappa=kappa_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_BSi_KD$par

    l_r_opt_BSi_KD <- optim(l_r_par, nLL_AFNSi_serr, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_BSi_KD$par

    if (abs(nLL_AFNSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random) - neg_loglikelihood) < tol_lik | (iter_count==max_iter)){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par))
}

CovEst_MI_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  X_t_fil <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t
  X_t_c_fil <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t_c
  S_t_fil <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t
  S_t_c_fil <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma, r)), length(c(delta, kappa, sigma, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_AFNSi_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_AFNSi_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  nabla_grad <- diag(c(rep(1, 4), sigma, r))

  Cov_est_orig <- t(nabla_grad) %*% Cov_est %*% nabla_grad

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  return(list(Cov_par_est = Cov_est_orig, St_err=list(delta=serr_est_orig[1], kappa=serr_est_orig[2:4], sigma=serr_est_orig[5:7], r1=serr_est_orig[8], r2=serr_est_orig[9], rc=serr_est_orig[10])))

}

# ============ - AFNS dependent model - =====================

nLL_AFNSd_serr <- function(delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol, l_r, mu_bar, X_rnd_nLL){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_AFNSd_serr_H <- function(vdParameters, mu_bar, X_rnd_H){

  n_factors <- 3
  # - Parameters
  delta <- vdParameters[1]
  kappa <- vdParameters[2:4] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[5:7] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol <- vdParameters[8:10] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[13]
  r_1 <- vdParameters[11]
  r_2 <- vdParameters[12]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_AFNSd_MI <- function(mu_bar, X_random, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  l_r_par <- log(r)

  iter_count <- 0
  repeat{
    neg_loglikelihood <- nLL_AFNSd_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_AFNSd_KD <- optim(delta_par, nLL_AFNSd_serr, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFNSd_KD$par

    kappa_opt_AFNSd_KD <- optim(kappa_par, nLL_AFNSd_serr, delta=delta_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFNSd_KD$par

    dg_l_Sigma_chol_opt_AFNSd_KD <- optim(dg_l_Sigma_chol_par, nLL_AFNSd_serr, delta=delta_par, kappa=kappa_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_AFNSd_KD$par

    odg_Sigma_chol_opt_AFNSd_KD <- optim(odg_Sigma_chol_par, nLL_AFNSd_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol_par <- odg_Sigma_chol_opt_AFNSd_KD$par

    l_r_opt_AFNSd_KD <- optim(l_r_par, nLL_AFNSd_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFNSd_KD$par

    if (abs(nLL_AFNSd_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol_par, l_r_par, mu_bar, X_rnd_nLL = X_random) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol=odg_Sigma_chol_par, l_r=l_r_par))
}

CovEst_MI_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  pe_dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  pe_odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  X_t_fil <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t
  X_t_c_fil <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t_c
  S_t_fil <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t
  S_t_c_fil <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r)), length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_AFNSd_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_AFNSd_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  ## - Build sub-jacobian for Sigma (use created function parest2cov in Est_fun_v0)

  subjac <- numDeriv::jacobian(parest2cov_jac, c(pe_dg_l_Sigma_chol, pe_odg_Sigma_chol), method="Richardson", side=NULL)

  Jac <- matrix(0, 13, 13)
  diag(Jac[1:4, 1:4]) <- rep(1, 4) # the untransformed delta
  Jac[5:10, 5:10] <- subjac
  diag(Jac[11:13, 11:13]) <- r

  Cov_est_orig <- t(Jac) %*% Cov_est %*% Jac

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  ## - Second Jacobian step to reconcile sigma_dg^2 to sigma_dg
  nabla_grad <- matrix(0, 13, 13)
  diag(nabla_grad[c(1:4, 8:13), c(1:4, 8:13)]) <- 1
  diag(nabla_grad[c(5:7), c(5:7)]) <- 0.5 / sigma_dg

  Cov_est_orig2 <- t(nabla_grad) %*% Cov_est_orig %*% nabla_grad

  serr_est_orig2 <- sqrt(diag(Cov_est_orig2))

  return(list(Cov_par_est = Cov_est_orig2, St_err=list(delta=serr_est_orig2[1], kappa=serr_est_orig2[2:4], sigma_dg=serr_est_orig2[5:7], Sigma_cov=serr_est_orig2[8:10], r1=serr_est_orig2[11], r2=serr_est_orig2[12], rc=serr_est_orig2[13])))

}

# ============ - AFGNS independent model - =====================

nLL_AFGNSi_serr <- function(delta, kappa, l_sigma, l_r, mu_bar, X_rnd_nLL){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- rep(0, n_ages)
  for(age in 1:n_ages){
    # - Original diag from Blackburn-Sherris (2013)
    H[age] <- exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:age))) / age
  }
  # H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, exp(l_sigma), delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age, delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(H)) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / H) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_AFGNSi_serr_H <- function(vdParameters, mu_bar, X_rnd_H){

  n_factors <- 5

  delta = vdParameters[c(1,2)]
  kappa = vdParameters[3:7]
  l_sigma = vdParameters[8:12]
  l_r = vdParameters[13:15]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, exp(l_sigma), delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age, delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((exp(l_sigma*2)) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_AFGNSi_MI <- function(mu_bar, X_random, delta, kappa, sigma, r, max_iter=200, tol_lik=0.1){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  l_sigma_par <- log(sigma)
  l_r_par <- log(r)

  iter_count <- 0
  repeat{
    neg_loglikelihood <- nLL_AFGNSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_AFGNSi_KD <- optim(delta_par, nLL_AFGNSi_serr, mu_bar=mu_bar, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFGNSi_KD$par

    kappa_opt_AFGNSi_KD <- optim(kappa_par, nLL_AFGNSi_serr, mu_bar=mu_bar, delta=delta_par, l_sigma=l_sigma_par, l_r=l_r_par, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFGNSi_KD$par

    l_sigma_opt_AFGNSi_KD <- optim(l_sigma_par, nLL_AFGNSi_serr, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_r=l_r_par, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_sigma_par <- l_sigma_opt_AFGNSi_KD$par

    l_r_opt_AFGNSi_KD <- optim(l_r_par, nLL_AFGNSi_serr, mu_bar=mu_bar, delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFGNSi_KD$par

    if (abs(nLL_AFGNSi_serr(delta_par, kappa_par, l_sigma_par, l_r_par, mu_bar, X_rnd_nLL = X_random) - neg_loglikelihood) < tol_lik | (iter_count==max_iter)){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, l_sigma=l_sigma_par, l_r=l_r_par))
}

CovEst_MI_AFGNSi <- function(x0, delta, kappa, sigma, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  X_t_fil <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t
  X_t_c_fil <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t_c
  S_t_fil <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t
  S_t_c_fil <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma, r)), length(c(delta, kappa, sigma, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_AFGNSi_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_AFGNSi_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$l_sigma, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  nabla_grad <- diag(c(rep(1, 7), sigma, r))

  Cov_est_orig <- t(nabla_grad) %*% Cov_est %*% nabla_grad

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  return(list(Cov_par_est = Cov_est_orig, St_err=list(delta=serr_est_orig[1:2], kappa=serr_est_orig[3:7], sigma=serr_est_orig[8:12], r1=serr_est_orig[13], r2=serr_est_orig[14], rc=serr_est_orig[15])))

}

# ============ - AFGNS dependent model - =====================

nLL_AFGNSd_serr <- function(delta, kappa, dg_l_Sigma_chol, odg_Sigma_chol1, odg_Sigma_chol2, l_r, mu_bar, X_rnd_nLL){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Parameters
  r_c <- l_r[3]
  r_1 <- l_r[1]
  r_2 <- l_r[2]

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_nLL

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age, delta[2])[2], B_AFNS(age, delta[1])[3], B_AFNS(age, delta[2])[3])
  }

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

nLL_AFGNSd_serr_H <- function(vdParameters, mu_bar, X_rnd_H){

  n_factors <- 5
  # - Parameters
  delta <- vdParameters[1:2]
  kappa <- vdParameters[3:7] # - take logs to ensure positivity
  dg_l_Sigma_chol <- vdParameters[8:12] # - log diagonal elements of the lower cholesky dec.
  odg_Sigma_chol1 <- vdParameters[13:17] # - off diagonal elements of the lower cholesky dec.
  odg_Sigma_chol2 <- vdParameters[18:22] # - off diagonal elements of the lower cholesky dec.
  r_c <- vdParameters[23]
  r_1 <- vdParameters[24]
  r_2 <- vdParameters[25]

  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  delta_matrix <- low_trg_fill(delta)

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t <- X_rnd_H

  H <- matrix(0, n_ages, n_ages) # - Age covariance
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  # - Likelihood tools
  log_lik <- rep(0, n_years)

  # - Setting H (covariance of measurement error - indep. among ages)
  H <- meas_err_BS(r_1, r_2, r_c, mu_bar)

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

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
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  for(t in 1:n_years){
    log_lik[t] <- sum(log(diag(H))) + sum(((mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1])^2) / diag(H)) +  #t(mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) %*% diag(1 / diag(H)) %*% (mu_bar[,t] - A_tT - B_tT %*% X_t[,t+1]) +
      sum(log(diag(R))) + t(X_t[,t+1] - Phi %*% X_t[,t]) %*% diag(1 / diag(R)) %*% (X_t[,t+1] - Phi %*% X_t[,t])
  }

  return(sum(log_lik))
}

co_asc_AFGNSd_MI <- function(mu_bar, X_random, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=200, tol_lik=0.1){

  # - Matrix for parameter estimates storage
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  # - Initialize log-likelihood
  neg_loglikelihood <- 0

  delta_par <- delta
  kappa_par <- kappa
  dg_l_Sigma_chol_par <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol1_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol[1:5]
  odg_Sigma_chol2_par <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol[6:10]
  l_r_par <- log(r)

  iter_count <- 0
  repeat{
    neg_loglikelihood <- nLL_AFGNSd_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol1_par, odg_Sigma_chol2_par, l_r_par, mu_bar, X_rnd_nLL = X_random)
    iter_count <- iter_count + 1

    delta_opt_AFGNSd_KD <- optim(delta_par, nLL_AFGNSd_serr, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    delta_par <- delta_opt_AFGNSd_KD$par

    kappa_opt_AFGNSd_KD <- optim(kappa_par, nLL_AFGNSd_serr, delta=delta_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    kappa_par <- kappa_opt_AFGNSd_KD$par

    dg_l_Sigma_chol_opt_AFGNSd_KD <- optim(dg_l_Sigma_chol_par, nLL_AFGNSd_serr, delta=delta_par, kappa=kappa_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    dg_l_Sigma_chol_par <- dg_l_Sigma_chol_opt_AFGNSd_KD$par

    odg_Sigma_chol1_opt_AFGNSd_KD <- optim(odg_Sigma_chol1_par, nLL_AFGNSd_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol1_par <- odg_Sigma_chol1_opt_AFGNSd_KD$par

    odg_Sigma_chol2_opt_AFGNSd_KD <- optim(odg_Sigma_chol2_par, nLL_AFGNSd_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, l_r=l_r_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    odg_Sigma_chol2_par <- odg_Sigma_chol2_opt_AFGNSd_KD$par

    l_r_opt_AFGNSd_KD <- optim(l_r_par, nLL_AFGNSd_serr, delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, mu_bar, X_rnd_nLL = X_random, gr = NULL, method = "Nelder-Mead", hessian = TRUE, control=list(maxit = 10000))
    l_r_par <- l_r_opt_AFGNSd_KD$par

    if (abs(nLL_AFGNSd_serr(delta_par, kappa_par, dg_l_Sigma_chol_par, odg_Sigma_chol1_par, odg_Sigma_chol2_par, l_r_par, mu_bar, X_rnd_nLL = X_random) - neg_loglikelihood) < tol_lik  | (iter_count==max_iter) ){
      break
    }
  }

  return(list(delta=delta_par, kappa=kappa_par, dg_l_Sigma_chol=dg_l_Sigma_chol_par, odg_Sigma_chol1=odg_Sigma_chol1_par, odg_Sigma_chol2=odg_Sigma_chol2_par, l_r=l_r_par))
}

CovEst_MI_AFGNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, D_se=50, max_it=200, tolerance_lev=0.1){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  pe_dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  pe_odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol

  X_t_fil <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t
  X_t_c_fil <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t_c
  S_t_fil <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t
  S_t_c_fil <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$S_t_c

  ## - Step 1.2 Get smoothed estimate which are the mean and variance of the normal distribution for sampling the states
  X_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$X_t_sm
  S_t_sm <- RTS_sm_bas(X_t_fil, X_t_c_fil, S_t_fil, S_t_c_fil, kappa, n_years)$S_t_sm

  # - Step 2 and 3: sample states (in iteration) and estimate
  ## - Step 2.1: Set D_se and define vector of sampled states
  #  D_se <- 50
  X_rnd <- matrix(0, n_factors, n_years+1)

  ## - Step 2.2: Initialise V_bar,
  V_bar <- matrix(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r)), length(c(delta, kappa, sigma_dg, Sigma_cov, r)))
  theta_bar <- rep(0, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - vector for the mean over the parameter estimates
  B_table <- matrix(NA, D_se, length(c(delta, kappa, sigma_dg, Sigma_cov, r))) # - Table for parameter storage

  for(d in 1:D_se){

    ## - Step 2.3: Sample states
    for(t in 1:(n_years+1)){
      X_rnd[,t] <- mvtnorm::rmvnorm(n = 1, mean=X_t_sm[,t], sigma = S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)])
    }

    # - Step 3: optimization and parameter storage (the covariance is not stored) - Call Coordinate Ascent
    p_opt_se <- co_asc_AFGNSd_MI(mu_bar=mu_bar, X_random=X_rnd, delta, kappa, sigma_dg, Sigma_cov, r, max_iter=max_it, tol_lik=tolerance_lev)
    hessian_comp <- numDeriv::hessian(nLL_AFGNSd_serr_H, c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol1, p_opt_se$odg_Sigma_chol2, p_opt_se$l_r), method="Richardson", mu_bar=mu_bar, X_rnd_H=X_rnd)

    V_bar <- V_bar + MASS::ginv(hessian_comp) / D_se
    B_table[d,] <- c(p_opt_se$delta, p_opt_se$kappa, p_opt_se$dg_l_Sigma_chol, p_opt_se$odg_Sigma_chol1, p_opt_se$odg_Sigma_chol2, p_opt_se$l_r)
    print(paste(d, "% imputation"))
  }

  # - Step 4: Estimation of the covariance matrix of the parameters

  B_cov <- cov(B_table) # (check whether it takes population covariance or sample covariance)

  Cov_est <- V_bar + B_cov * (1 + 1/D_se)

  serr_est <- sqrt(diag(Cov_est))

  # - Step 5: Back-transformation for the original parameters

  ## - Build sub-jacobian for Sigma (use created function parest2cov in Est_fun_v0)

  subjac <- numDeriv::jacobian(parest2cov_jac, c(pe_dg_l_Sigma_chol, pe_odg_Sigma_chol), method="Richardson", side=NULL)

  Jac <- matrix(0, 25, 25)
  diag(Jac[1:7, 1:7]) <- rep(1, 7) # the untransformed delta
  Jac[8:22, 8:22] <- subjac
  diag(Jac[23:25, 23:25]) <- r

  Cov_est_orig <- t(Jac) %*% Cov_est %*% Jac

  serr_est_orig <- sqrt(diag(Cov_est_orig))

  ## - Second Jacobian step to reconcile sigma_dg^2 to sigma_dg
  nabla_grad <- matrix(0, 25, 25)
  diag(nabla_grad[c(1:7, 13:25), c(1:7, 13:25)]) <- 1
  diag(nabla_grad[c(8:12), c(8:12)]) <- 0.5 / sigma_dg

  Cov_est_orig2 <- t(nabla_grad) %*% Cov_est_orig %*% nabla_grad

  serr_est_orig2 <- sqrt(diag(Cov_est_orig2))

  return(list(Cov_par_est = Cov_est_orig2, St_err=list(delta=serr_est_orig2[1:2], kappa=serr_est_orig2[3:7], sigma_dg=serr_est_orig2[8:12], Sigma_cov=serr_est_orig2[13:22], r1=serr_est_orig2[23], r2=serr_est_orig2[24], rc=serr_est_orig2[25])))

}
