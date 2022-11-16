#################################################################################

# - Filtering and smoothing

#################################################################################


#============================== - Filtering functions and smoothing functions for EM - ==============================


## - Kalman filtering function
### - Independent BS
KF_BSi_uKD <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_BSi(age, sigma, delta)
    B_tT[age,] <- B_BSi(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Dependent BS

#### - Univariate Koopman-Durbin
KF_BSd_2F_uKD <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- diag(1, n_factors) * 1e-10

  S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
  S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      #      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + max(1e-8, r_c + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i) #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
  }

  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

KF_BSd_3F_uKD <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- diag(1, n_factors) * 1e-10

  S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
  S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
  }

  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

# - Gompertz-Makeham
### - Independent BS
KF_GMki_uKD <- function(x0, delta, kappa, gamma, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_GMki(age, sigma, delta, gamma)
    B_tT[age,] <- B_GMki(age,delta, gamma)
  }

  Phi <- diag(exp(-kappa), n_factors)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

#### - Univariate Koopman-Durbin
KF_GMkd_uKD <- function(x0, delta, kappa, gamma, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  delta_matrix <- low_trg_fill(delta) # - create lower diagonal matrix

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- diag(1, n_factors) * 1e-10

  S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
  S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    A_tT[age,1] <- A_GMkd(age, Low_chol, delta_matrix, gamma)####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- B_GMkd(age, delta_matrix, gamma)  ###### B_ind(age,delta)
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      #      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + max(1e-8, r_c + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i) #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
  }

  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}


### - Independent AFNS
KF_AFNSi_uKD <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Dependent AFNS
KF_AFNSd_uKD <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- diag(1, n_factors) * 1e-10

  S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
  S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
  }

  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Independent AFNS
KF_AFGNSi_uKD <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFGNSi(age, sigma, delta[1], delta[2])
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Dependent AFNS
KF_AFGNSd_uKD <- function(x0, delta, kappa, sigma_dg, Sigma_cov, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- diag(1, n_factors) * 1e-10

  S_t_c[,(1:n_factors)] <- diag(1, n_factors) * 1e-10
  S_t[,(1:n_factors)] <- diag(1, n_factors) * 1e-10

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

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
    A_tT[age,1] <- A_AFGNSg(age, Low_chol, delta[1], delta[2])####### A_ind(age, exp(l_sigma), delta)
    B_tT[age,] <- c(B_AFNS(age,delta[1])[c(1,2)], B_AFNS(age,delta[2])[2], B_AFNS(age,delta[1])[3], B_AFNS(age,delta[2])[3])
  }

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))    #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)
    #    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + max(1e-8, r_c + exp(r_1) * exp(exp(r_2))) #max(0.000001, r_c + r_1 * exp(r_2))#exp(r_1 + r_2)  In case of the error specification of Xu et al. (2019)

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti
  }

  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Independent AFUNS
KF_AFUNSi_uKD <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFUNSi(age, sigma, delta)
    B_tT[age,] <- B_AFUNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - Independent AFRNS
KF_AFRNSi_uKD <- function(x0, delta, kappa, sigma, r, mu_bar){

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_factors <- length(kappa)

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  v_ti <- mu_bar
  F_ti <- mu_bar

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <-1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  # - Initialize X and Sigma
  x_ti <- x0 #init_X
  P_ti <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors) # - Factor covariance

  for(age in 1:n_ages){    # - scroll over the ages
    A_tT[age,1] <- A_AFRNSi(age, sigma, delta)
    B_tT[age,] <- B_AFRNS(age,delta)
  }

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)
  R <- diag(((sigma^2) / (2 * kappa)) * (1 - exp(-2 * kappa)), n_factors)

  for(t in 1:n_years){

    # - First observation
    x_ti <- Phi %*% x_ti    # - x_{1,t}
    P_ti <- Phi %*% P_ti %*% t(Phi) + R     # - P_{1,t}
    v_ti[1,t] <- mu_bar[1,t] - A_tT[1] - B_tT[1,] %*% x_ti
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t] # * (v_ti[i-1,t] / F_ti[i-1,t])   ### - Same change from + to minus following Xu, Sherris and Ziveyi 2019SAJ               # - a_{t,i+1}

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti) #exp(r_1 + r_2 * (i-1)) %*% t(K_ti) ####(max(0.001, r_c + r_1 * sum(exp(r_2 * c(1:(i-1)))) / (i-1))) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }

    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

### - CIR
KF_CIR_uKD <- function(x0, delta, kappa, sigma, theta_P, r, mu_bar){
  n_factors <- length(kappa)

  r_1 <- log(r[1])
  r_2 <- log(r[2])
  r_c <- log(r[3])

  n_ages <- nrow(mu_bar)   # - Number of ages
  n_years <- ncol(mu_bar)  # - Number of years

  v_ti <- mu_bar
  F_ti <- mu_bar

  # - Def. variables
  ## - State variables
  X_t <- matrix(NA, n_factors, (n_years+1)) # - state (including state at t=0)
  X_t_c <- matrix(NA, n_factors, (n_years+1)) # - conditional state (including state at t=0; conditional to t-1)

  S_t <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - covariance matrix of factors
  S_t_c <- matrix(NA, n_factors, n_factors * (n_years + 1)) # - cond'l covariance matrix of factors

  ## - Factor loading matrices
  A_tT <- matrix(0, n_ages, 1)
  B_tT <- matrix(NA, n_ages, n_factors)

  # - Initial values of states and of covariance
  X_t_c[,1] <- x0
  X_t[,1] <- x0

  S_t_c[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)
  S_t[,(1:n_factors)] <- 1e-10 * diag(1, n_factors)

  R <- matrix(0, n_factors, n_factors * (n_years + 1)) # - Factor covariance

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
    F_ti[1,t] <- B_tT[1,] %*% P_ti %*% B_tT[1,] + exp(r_c) + exp(r_1 + exp(r_2))

    X_t_c[,t+1] <- x_ti
    S_t_c[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

    for(i in 2:n_ages){
      x_ti <- x_ti + P_ti %*% B_tT[i-1,] %*% (1 / F_ti[i-1,t]) %*% v_ti[i-1,t]
      x_ti <- l_bound(x_ti)

      # - Joseph formula univariate
      K_ti <- P_ti %*% B_tT[i-1,] / F_ti[i-1,t]
      P_ti <- (diag(1, n_factors) - K_ti %*% B_tT[i-1,]) %*% P_ti %*% t(diag(1, n_factors) - K_ti %*% B_tT[i-1,]) + K_ti %*% (exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:(i-1)))) / (i-1)) %*% t(K_ti)

      # - log-likelihood values from Koopman and Durbin
      F_ti[i,t] <- B_tT[i,] %*% P_ti %*% B_tT[i,] + exp(r_c) + exp(r_1) * sum(exp(exp(r_2) * c(1:i))) / i #exp(r_1 + r_2 * i)
      v_ti[i,t] <- mu_bar[i,t] - A_tT[i] - B_tT[i,] %*% x_ti

    }
    X_t[,t+1] <- x_ti
    S_t[,(t * n_factors + 1):((t+1) * n_factors)] <- P_ti

  }
  return(list(X_t=X_t, X_t_c=X_t_c, S_t=S_t, S_t_c=S_t_c))
}

# - Rauch-Tung-Streibel smoother (Valid for all models)
RTS_sm_bas <- function(X_t, X_t_c, S_t, S_t_c, kappa, n_years){

  n_factors <- length(kappa)
  X_t_sm <- matrix(NA, n_factors, (n_years+1))
  S_t_sm <- matrix(NA, n_factors, n_factors * (n_years + 1))
  G_t <- matrix(NA, n_factors, n_factors * n_years)   # G_t goes is for t = 0,...,T-1

  Phi <- diag(exp(-kappa), n_factors) # K_p <- diag(kappa, 2) ## exp(-K_p)

  # - Set last value
  X_t_sm[,n_years+1] <- X_t[,n_years+1] # because it goes from t=0 to t=T, hence this line reads X_T_sm = X_T
  S_t_sm[,(n_factors * n_years + 1):(n_factors * (n_years + 1))] <- S_t[,(n_factors * n_years + 1):(n_factors * (n_years + 1))]

  # - Smoothing recursions
  for(t in n_years:1){

    G_t[,((t-1) * n_factors + 1):(t * n_factors)] <- S_t[,((t-1) * n_factors + 1):(t * n_factors)] %*% t(Phi) %*% solve(S_t_c[,(t * n_factors + 1):((t+1) * n_factors)])
    X_t_sm[,t] <- X_t[,t] + G_t[,((t-1) * n_factors + 1):(t * n_factors)] %*% (X_t_sm[,t+1] - X_t_c[,t+1]) # - Update step
    S_t_sm[,((t-1) * n_factors + 1):(t * n_factors)] <- S_t[,((t-1) * n_factors + 1):(t * n_factors)] + G_t[,((t-1) * n_factors + 1):(t * n_factors)] %*% (S_t_sm[,(t * n_factors + 1):((t+1) * n_factors)] - S_t_c[,(t * n_factors + 1):((t+1) * n_factors)]) %*% t(G_t[,((t-1) * n_factors + 1):(t * n_factors)])
  }
  return(list(X_t_sm=X_t_sm, S_t_sm=S_t_sm))
}



