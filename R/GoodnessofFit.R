#source('FilterSmoother.R')

#==========================================
# - Fitted rates (mu_bar_hat)
#==========================================

# - Blackburn-Sherris
## - Independent
mu_bar_hat_BSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_BSd_2F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- 2
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  KF_outcome <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

mu_bar_hat_BSd_3F <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- 3
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  KF_outcome <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - Gompertz-Makeham
## - Independent
mu_bar_hat_GMki <- function(x0, delta, kappa, gamma, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_GMki_uKD(x0, delta, kappa, gamma, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, exp(l_sigma), delta, exp(l_gamma))####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_GMki(age, delta, exp(l_gamma))  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_GMkd <- function(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  KF_outcome <- KF_GMkd_uKD(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMkd(age, Low_chol, delta_matrix, exp(l_gamma))####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_GMkd(age, delta_matrix, exp(l_gamma))  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - AFNS
## - Independent
mu_bar_hat_AFNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)
    B_tT[age,] <- B_AFNS(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_AFNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - AFGNS
## - Independent
mu_bar_hat_AFGNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, sigma, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - Dependent
mu_bar_hat_AFGNSd <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - AFUNS
## - Independent
mu_bar_hat_AFUNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFUNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFUNSi(age, sigma, delta)
    B_tT[age,] <- B_AFUNS(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

# - AFRNS
## - Independent
mu_bar_hat_AFRNSi <- function(x0, delta, kappa, sigma, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_AFRNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFRNSi(age, sigma, delta)
    B_tT[age,] <- B_AFRNS(age,delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}

## - CIR
mu_bar_hat_CIR <- function(x0, delta, kappa, sigma, theta_P, r, mu_bar){
  n_factors <- length(kappa)
  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  KF_outcome <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_P, r, mu_bar)
  X_t <- KF_outcome$X_t

  mu_bar_hat <- matrix(NA, n_ages, n_years)
  rownames(mu_bar_hat) <- AgeRange
  colnames(mu_bar_hat) <- CohortRange

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_CIR(age, (theta_P * kappa / delta), sigma, delta)
    B_tT[age,] <- B_CIR(age, sigma, delta)
  }

  for(t in 1:n_years){
    mu_bar_hat[,t] <- A_tT + B_tT %*% X_t[,t+1]
  }
  return(mu_bar_hat)
}


#==========================================================
# - Model Selection (AIC and BIC) - Build cross-validation
#==========================================================

AIC_BIC <- function(log_likelihood, n_par, n_obs){

  AIC <- - 2 * log_likelihood + 2 * n_par
  BIC <- - 2 * log_likelihood + n_par * log(n_obs)

  return(list(AIC=AIC, BIC=BIC))
}

#=======================================================
# - Fitted rates (useful for residuals plot and RMSE)
#=======================================================

# - mu_bar hat

# - Residuals
residuals_f <- function(observed, estimated){
  residuals <- observed - estimated
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

# - Standardized residuals
residuals_std_BSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_BSi(age,delta)
  }

  residuals <- observed - mu_bar_hat_BSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)
  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)
  #Chol_dec <- chol(Var_y)
  premultiplier <- solve(expm::sqrtm(Var_y))#solve(t(Chol_dec))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_BSd_2F <- function(observed, x0, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)
  }

  residuals <- observed - mu_bar_hat_BSd_3F(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)
  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_BSd_3F <- function(observed, x0, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)
  }

  residuals <- observed - mu_bar_hat_BSd_3F(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)
  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_GMki <- function(observed, x0, delta, kappa, gamma, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_GMki(age,delta, gamma)
  }

  residuals <- observed - mu_bar_hat_GMki(x0, delta, kappa, gamma, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)
  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)
  #Chol_dec <- chol(Var_y)
  premultiplier <- solve(expm::sqrtm(Var_y))#solve(t(Chol_dec))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_GMkd <- function(observed, x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_GMkd(age, delta_matrix, gamma)
  }

  residuals <- observed - mu_bar_hat_GMkd(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)
  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFNSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFNS(age,delta)
  }

  residuals <- observed - mu_bar_hat_AFNSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFNSd <- function(observed, x0, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFNS(age,delta)
  }

  residuals <- observed - mu_bar_hat_AFNSd(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)
  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFGNSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  residuals <- observed - mu_bar_hat_AFGNSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFGNSd <- function(observed, x0, delta, kappa, sigma_dg, Sigma_cov, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  # - Build diffusion process
  ## - Build lower cholesky factor
  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)
  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  R_matrix <- matrix(0, n_factors, n_factors)
  # - Get Sigma (covariance matrix of the diffusion process)
  Sigma_diff <- Low_chol %*% t(Low_chol) # matrix(0, n_factors, n_factors)

  # - Get R (covariance of the state variable)
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R_matrix[row,col] <- Sigma_diff[row,col] * (1 - exp(- kappa[row] - kappa[col])) / (kappa[row] + kappa[col])
    }
  }

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  residuals <- observed - mu_bar_hat_AFGNSd(x0, delta, kappa, sigma_dg, Sigma_cov, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)
  Var_y <- H + B_tT %*% R_matrix %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFRNSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFRNS(age,delta)
  }

  residuals <- observed - mu_bar_hat_AFRNSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_AFUNSi <- function(observed, x0, delta, kappa, sigma, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_AFUNS(age,delta)
  }

  residuals <- observed - mu_bar_hat_AFUNSi(x0, delta, kappa, sigma, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  Var_y <- H + B_tT %*% diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors) %*% t(B_tT)

  premultiplier <- solve(expm::sqrtm(Var_y))

  for(t in 1:n_years){
    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

residuals_std_CIR <- function(observed, x0, delta, kappa, sigma, theta_P, r){
  n_factors <- length(kappa)
  n_ages <- nrow(observed)
  n_years <- ncol(observed)

  B_tT <- matrix(NA, n_ages, n_factors)

  l_r1 <- log(r[1])
  l_r2 <- log(r[2])
  l_rc <- log(r[3])

  B_tT <- matrix(NA, n_ages, n_factors)
  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    B_tT[age,] <- B_CIR(age, sigma, delta)
  }

  residuals <- observed - mu_bar_hat_CIR(x0, delta, kappa, sigma, theta_P, r, observed)
  H <- meas_err_BS(l_r1, l_r2, l_rc, observed)

  X_t <- KF_CIR_uKD(x0, delta, kappa, sigma, theta_P, r, observed)$X_t
  P_ti <- 1e-10 * diag(1, n_factors)

  Var_y <- matrix(0, n_years, n_years)


  for(t in 1:n_years){
    R <- diag((sigma^2) * ((1 - exp(-kappa)) / kappa) * (0.5 * theta_P * (1 - exp(-kappa)) + exp(-kappa) * X_t[,t+1]), n_factors)
    Var_y <- H + (B_tT %*% R %*% t(B_tT))
    premultiplier <- solve(expm::sqrtm(Var_y))

    residuals[,t] <- premultiplier %*% residuals[,t]
  }

  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}


# - Poisson standardized residuals

poi_r_std <- function(mu_bar_hat, deaths, exposures){
  mu_hat <- mu_bar_hat
  for(i in 2:nrow(mu_hat)){
    mu_hat[i,] <- i * mu_bar_hat[i,] - (i-1) * mu_bar_hat[i-1,]
  }
  D_hat <- mu_hat * exposures

  poi_res_std <- (deaths - D_hat) / sqrt(D_hat)

  return(poi_res_std)
}

#=======================================================
# - Probability of negative rates
#=======================================================

pr_neg_rates_f_BSi <- function(n_sim=100000, x0, delta, kappa, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_BSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean_BSi <- exp(- kappa * yrs_proj) * X_t_last
  var_BSi <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_BSi, var_BSi)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age, delta)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_BSi <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_BSi)
}

pr_neg_rates_f_BSd_2F <- function(n_sim=100000, x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, yrs_proj=1){

  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  Sigma <- low_trg_fill_0diag(Sigma_cov) + t(low_trg_fill_0diag(Sigma_cov)) + diag(sigma_dg^2)

  X_t_last <- KF_BSd_2F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  R <- Sigma
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma[row,col] * (1 - exp(- (kappa[row] + kappa[col]) * yrs_proj)) / (kappa[row] + kappa[col])
    }
  }

  mean_BSd <- exp(- kappa * yrs_proj) * X_t_last

  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_BSd, R)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]# + rnorm(1, 0, sd=sqrt(exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age))
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_BSd <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_BSd)
}

pr_neg_rates_f_BSd_3F <- function(n_sim=100000, x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, yrs_proj=1){

  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  Sigma <- low_trg_fill_0diag(Sigma_cov) + t(low_trg_fill_0diag(Sigma_cov)) + diag(sigma_dg^2)

  X_t_last <- KF_BSd_3F_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  R <- Sigma
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma[row,col] * (1 - exp(- (kappa[row] + kappa[col]) * yrs_proj)) / (kappa[row] + kappa[col])
    }
  }

  mean_BSd <- exp(- kappa * yrs_proj) * X_t_last

  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_BSd, R)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_3F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_3F(age, delta_matrix)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]# + rnorm(1, 0, sd=sqrt(exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age))
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_BSd <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_BSd)
}

pr_neg_rates_f_GMki <- function(n_sim=100000, x0, delta, kappa, gamma, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_GMki_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean_BSi <- exp(- kappa * yrs_proj) * X_t_last
  var_BSi <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_BSi, var_BSi)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, sigma, delta, gamma)
    B_tT[age,] <- B_GMki(age, delta, gamma)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x)
}

pr_neg_rates_f_GMkd <- function(n_sim=100000, x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar, yrs_proj=1){

  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  Sigma <- low_trg_fill_0diag(Sigma_cov) + t(low_trg_fill_0diag(Sigma_cov)) + diag(sigma_dg^2)

  X_t_last <- KF_GMkd_uKD(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  R <- Sigma
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma[row,col] * (1 - exp(- (kappa[row] + kappa[col]) * yrs_proj)) / (kappa[row] + kappa[col])
    }
  }

  mean_BSd <- exp(- kappa * yrs_proj) * X_t_last

  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_BSd, R)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSd_2F(age, Low_chol, delta_matrix)
    B_tT[age,] <- B_BSd_2F(age, delta_matrix)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]# + rnorm(1, 0, sd=sqrt(exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age))
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_BSd <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_BSd)
}

pr_neg_rates_f_AFNSi <- function(n_sim=100000, x0, delta, kappa, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_AFNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean_AFNSi <- exp(- kappa * yrs_proj) * X_t_last
  var_AFNSi <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_AFNSi, var_AFNSi)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSi(age, sigma, delta)
    B_tT[age,] <- B_AFNS(age, delta)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_AFNSi <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_AFNSi)
}

pr_neg_rates_f_AFNSd <- function(n_sim=100000, x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, yrs_proj=1){

  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  Sigma <- low_trg_fill_0diag(Sigma_cov) + t(low_trg_fill_0diag(Sigma_cov)) + diag(sigma_dg^2)

  X_t_last <- KF_AFNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  R <- Sigma
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma[row,col] * (1 - exp(- (kappa[row] + kappa[col]) * yrs_proj)) / (kappa[row] + kappa[col])
    }
  }

  mean_AFNSd <- exp(- kappa * yrs_proj) * X_t_last

  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean_AFNSd, R)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFNSg(age, Low_chol, delta) ####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFNS(age, delta)  ###### B_ind(age,delta)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]# + rnorm(1, 0, sd=sqrt(exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age))
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x_AFNSd <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x_AFNSd)
}

pr_neg_rates_f_AFGNSi <- function(n_sim=100000, x0, delta, kappa, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_AFGNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean <- exp(- kappa * yrs_proj) * X_t_last
  var <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean, var)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, sigma, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x)
}

pr_neg_rates_f_AFGNSd <- function(n_sim=100000, x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar, yrs_proj=1){

  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  Sigma <- low_trg_fill_0diag(Sigma_cov) + t(low_trg_fill_0diag(Sigma_cov)) + diag(sigma_dg^2)

  X_t_last <- KF_AFGNSd_uKD(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar)$X_t[,n_years+1]

  R <- Sigma
  for(row in 1:n_factors){
    for(col in 1:n_factors){
      R[row,col] <- Sigma[row,col] * (1 - exp(- (kappa[row] + kappa[col]) * yrs_proj)) / (kappa[row] + kappa[col])
    }
  }

  mean <- exp(- kappa * yrs_proj) * X_t_last

  dg_l_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$dg_l_Sigma_chol
  odg_Sigma_chol <- cov2par(c(sigma_dg^2, Sigma_cov))$odg_Sigma_chol
  Low_chol <- low_trg_fill_0diag(odg_Sigma_chol)
  diag(Low_chol) <- exp(dg_l_Sigma_chol)

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean, R)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]# + rnorm(1, 0, sd=sqrt(exp(rc) + exp(r1) * sum(exp(exp(r2) * c(1:age))) / age))
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x)
}

pr_neg_rates_f_AFRNSi <- function(n_sim=100000, x0, delta, kappa, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_AFRNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean <- exp(- kappa * yrs_proj) * X_t_last
  var <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean, var)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFRNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFRNS(age,delta)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x)
}

pr_neg_rates_f_AFUNSi <- function(n_sim=100000, x0, delta, kappa, sigma, r, mu_bar, yrs_proj=1){
  n_years <- ncol(mu_bar)
  n_ages <- nrow(mu_bar)
  n_factors <- length(kappa)

  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  X_t_last <- KF_AFUNSi_uKD(x0, delta, kappa, sigma, r, mu_bar)$X_t[,n_years+1]

  mean <- exp(- kappa * yrs_proj) * X_t_last
  var <- diag((sigma^2) / (2 * kappa) * (1 - exp(-2 * kappa * yrs_proj)))

  X_rnd <- matrix(NA, n_factors, n_sim)

  X_rnd <- rmvnorm(n = n_sim, mean, var)

  mu_rates_tb <- matrix(NA, n_ages, n_sim)

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFUNSi(age, exp(l_sigma), delta)
    B_tT[age,] <- B_AFUNS(age,delta)
    for(sim in 1:n_sim){
      mu_rates_tb[age,sim] <- A_tT[age,1] + t(B_tT[age,]) %*% X_rnd[sim,]
    }
  }

  mu_rates_tb_ind <- ifelse(mu_rates_tb<0,1,0)

  mu_pcg_x <- 100 * rowSums(mu_rates_tb_ind) / n_sim

  return(mu_pcg_x)
}







