
#====================== - Projected survival curves - ================

# - Blackburn - Sherris independent factor model
S_t_BSi_proj <- function(x0, delta, kappa, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age, delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - Blackburn - Sherris dependent factor model
S_t_BSd_2F_proj <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  delta_matrix <- low_trg_fill(delta)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

S_t_BSd_3F_proj <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  delta_matrix <- low_trg_fill(delta)

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - Blackburn - Sherris independent factor model
S_t_GMki_proj <- function(x0, delta, kappa, gamma, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_GMki_uKD(x0, delta, kappa, gamma, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, sigma, delta, gamma)
    B_tT[age,] <- B_GMki(age, delta, gamma)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - Blackburn - Sherris dependent factor model
S_t_GMkd_proj <- function(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  delta_matrix <- low_trg_fill(delta)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_GMkd_uKD(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMkd(age, Low_chol, delta_matrix, gamma)
    B_tT[age,] <- B_GMkd(age, delta_matrix, gamma)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - AFNS independent factor model
S_t_AFNSi_proj <- function(x0, delta, kappa, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)
    B_tT[age,] <- B_AFNS(age,delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - AFNS dependent factor model
S_t_AFNSd_proj <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)
    B_tT[age,] <- B_AFNS(age,delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - AFGNS independent factor model
S_t_AFGNSi_proj <- function(x0, delta, kappa, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, sigma, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

S_t_AFUNSi_proj <- function(x0, delta, kappa, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFUNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFUNSi(age, sigma, delta)
    B_tT[age,] <- B_AFUNS(age,delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

S_t_AFRNSi_proj <- function(x0, delta, kappa, sigma, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFRNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFRNSi(age, sigma, delta)
    B_tT[age,] <- B_AFRNS(age,delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - AFGNS dependent factor model
S_t_AFGNSd_proj <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}

# - CIR model
S_t_CIR_proj <- function(x0, delta, kappa, sigma, theta_P, r, mu_bar, proj_years){

  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)
  n_years <- ncol(mu_bar)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)
  X_t_last <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_P, r, mu_bar)$X_t[,n_years+1]

  # - Partial sum from j=0 to (proj years - 1)
  matrix_sum <- matrix(1, n_factors, proj_years)
  for(year in 1:proj_years){
    matrix_sum[,year] <- exp(-(year-1) * kappa)
  }
  summatory <- rowSums(matrix_sum)

  E_X_t1 <- exp(-kappa * proj_years) * X_t_last + theta_P * (1 - exp(-kappa)) * summatory #exp(-kappa) * X_t_last + theta_P * (1 - exp(-kappa)) # - FIX THE RECURSION TO GENERALIZE

  S_projected <- matrix(NA, n_ages, 1)
  mu_hat_prj <- matrix(NA, n_ages, 1)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, (theta_P * kappa / delta), sigma, delta)
    B_tT[age,] <- B_CIR(age, sigma, delta)
    S_projected[age, 1] <- exp(-age * (A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1))
    mu_hat_prj[age,1] <- A_tT[age,1] + t(B_tT[age,]) %*% E_X_t1
  }
  return(S_projected)
}



