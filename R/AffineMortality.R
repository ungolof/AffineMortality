# - Functions for AffineMortality package

#source("Coordinate_Ascent.R")
#source("Est_fun.R")
#source("FilterSmoother.R")
#source("GoodnessofFit.R")
#source("Projection.R")
#source("StErr_Bootstrap.R")
#source("StErr_MultImp.R")


## - Fit
# - Starting values must be provided by the user in the format of a list, eg.: list(x0=c(...))
## - eventually add the option for using subplex optimization
#' @title affine_fit
#'
#' @description Estimation of the parameters of affine mortality models using the gradient-free Nelder-Mead algorithm for groups of parameters.
#'
#' @param model Specific model to be fit. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence (fact_dep=FALSE)
#' @param n_factors Number of factors. For some models, these are set by default (e.g. n_factors=3 for AFNS models)
#' @param data Table with the average mortality rates (ages on the rows and years on the columns)
#' @param st_val Starting value for the parameters to be supplied as list
#' @param max_iter Maximum number of iterations for the Coordinate Ascent estimation algorithm in absence of convergence. Default value set to 200
#' @param tolerance Minimum value of increase of the log-likelihood value before convergence. Default value set to 0.1
#' @param wd Working directory for saving the partial output of the estimation process
#'
#' @return A list with components:
#' \item{model}{ Name of the model}
#' \item{fit}{ List with parameter estimates `par_est`, log-likelihood function value `log_lik`, and table of their value at each iteration `CA_table`}
#' \item{n.parameters}{ Number of parameters of the model}
#' \item{AIC}{ Value of the Akaike Information Criterion of the model}
#' \item{BIC}{ Value of the Bayesian Information Criterion of the model}

#' @examples
#' data(mu_bar) # - load the age-cohort US dataset of males aged 50-99 born between 1883 and 1915, and load the default starting values (`sv_default`)
#' starting_values <- sv_default$BSd # - list of default starting values as provided by the package
#' pe_BSd_3F <- affine_fit(model="BS", fact_dep=TRUE, n_factors=3, data=mu_bar, st_val=starting_values, max_iter=5, tolerance=0.1, wd="working_folder_directory")

#' @export
affine_fit <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, data, st_val, max_iter=200, tolerance=0.1, wd=0){

  if(model=="AFNS"){
    if(fact_dep==TRUE){
      fit <- co_asc_AFNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
      return(list(model=model, fit=fit, n.parameters=16, AIC=AIC_BIC(fit$log_lik, 16, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 16, (nrow(data) * ncol(data)))$BIC))
    }
    else{
      fit <- co_asc_AFNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
      return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
    }
  }else{
    if(model=="BS"){
      if(fact_dep==TRUE){
        if(n_factors==2){
          fit <- co_asc_BSd_2F(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
          return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
        } else{
          fit <- co_asc_BSd_3F(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
          return(list(model=model, fit=fit, n.parameters=21, AIC=AIC_BIC(fit$log_lik, 21, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 21, (nrow(data) * ncol(data)))$BIC))
        }
      }else{ #i.e. factor independence
        fit <- co_asc_BSi(mu_bar = data, x0=st_val$x0[1:n_factors], delta=st_val$delta[1:n_factors], kappa=st_val$kappa[1:n_factors], sigma = st_val$sigma[1:n_factors], r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
        return(list(model=model, fit=fit, n.parameters=(3 + n_factors*4), AIC=AIC_BIC(fit$log_lik, (3 + n_factors*4), (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, (3 + n_factors*4), (nrow(data) * ncol(data)))$BIC))
      }}else{
        if(model=="CIR"){
          fit <- co_asc_CIR(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma=st_val$sigma, theta_P=st_val$theta_P, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
          return(list(model=model, fit=fit, n.parameters=(3 + n_factors*5), AIC=AIC_BIC(fit$log_lik, (3 + n_factors*5), (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, (3 + n_factors*5), (nrow(data) * ncol(data)))$BIC))
        }else{
          if(model=="AFUNS"){
            if(fact_dep==TRUE){
              fit <- co_asc_AFUNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
              return(list(model=model, fit=fit, n.parameters=18, AIC=AIC_BIC(fit$log_lik, 18, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 18, (nrow(data) * ncol(data)))$BIC))
            }else{
              fit <- co_asc_AFUNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
              return(list(model=model, fit=fit, n.parameters=15, AIC=AIC_BIC(fit$log_lik, 15, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 15, (nrow(data) * ncol(data)))$BIC))
            }}else{
              if(model=="AFRNS"){
                if(fact_dep==TRUE){
                  fit <- co_asc_AFRNSd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                  return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
                }else{
                  fit <- co_asc_AFRNSi(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                  return(list(model=model, fit=fit, n.parameters=10, AIC=AIC_BIC(fit$log_lik, 10, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 10, (nrow(data) * ncol(data)))$BIC))
                }
              }else{
                if(model=="GMk"){
                  if(fact_dep==TRUE){
                    fit <- co_asc_GMkd(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, gamma=st_val$gamma, sigma_dg = st_val$sigma_dg, Sigma_cov = st_val$Sigma_cov, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                    return(list(model=model, fit=fit, n.parameters=13, AIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 13, (nrow(data) * ncol(data)))$BIC))
                  }else{
                    fit <- co_asc_GMki(mu_bar = data, x0=st_val$x0, delta=st_val$delta, kappa=st_val$kappa, sigma = st_val$sigma, r=c(st_val$r1, st_val$r2, st_val$rc), max_iter=max_iter, tol_lik=tolerance, workdir=wd)
                    return(list(model=model, fit=fit, n.parameters=11, AIC=AIC_BIC(fit$log_lik, 11, (nrow(data) * ncol(data)))$AIC, BIC=AIC_BIC(fit$log_lik, 11, (nrow(data) * ncol(data)))$BIC))
                  }
                }
              }


            }
        }
      }}}


#======================== - Filtering distribution - ===================================

# - Filtering distribution (xfilter)
#' @title xfilter
#'
#' @description Estimation of mean and variance of the latent states using the univariate Kalman filtering procedure of Koopman and Durbin (2000)
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Value of the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return A list with the following components:
#' \item{X_t}{ Mean of the latent process X(t) for the update step}
#' \item{X_t_c}{ Mean of the latent process X(t) for the prediction step}
#' \item{S_t}{ Covariance matrix of the latent process X(t) for the update step}
#' \item{S_t_c.}{ Covariance matrix of the latent process X(t) for the prediction step}

#' @examples
#' X_filtered <- xfilter(model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar)

#' @export
xfilter <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      filter <- KF_AFNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(filter)
    } else{
      filter <- KF_AFNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(filter)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            filter <- KF_BSd_2F_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          } else{
            filter <- KF_BSd_3F_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          }
      } else{ # - it is the Blackburn-Sherris model with independent factors
        filter <- KF_BSi_uKD(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(filter)

      }
    } else{
      if(model=="CIR"){
        filter <- KF_CIR_uKD(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(filter)
      } else{
        if(model=="AFGNS"){# - if none of the previous model was selected, then it is an AFGNS
          if(fact_dep==TRUE){
            filter <- KF_AFGNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          } else{
            filter <- KF_AFGNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(filter)
          }} else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                filter <- KF_AFUNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(filter)
              } else{
                filter <- KF_AFUNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(filter)
              }
            }else{
              if(model=="AFRNS"){
                if(fact_dep==TRUE){
                  filter <- KF_AFRNSd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(filter)
                } else{
                  filter <- KF_AFRNSi_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(filter)
                }
              } else{
                if(model=="GMk"){
                  if(fact_dep==TRUE){
                    filter <- KF_GMkd_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(filter)
                  } else{
                    filter <- KF_GMki_uKD(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(filter)
                  }
                }
              }
            }
          }
      }
    }
  }
}


#======================== - Smoothing - ===================================
#' @title xsmooth
#'
#' @description Estimation of mean and covariance of the distribution of the smoothed latent process X(t) using the Rauch-Tung-Striebel procedure
#'
#' @param filterobject A `xfilter` object
#' @param kappa parameter kappa from the real-world dynamics of the latent variables
#'
#' @return List with mean `X_t_sm` and covariance `S_t_sm` of the distribution of the latent process X(t)
#' @examples
#' X_smoothed <- xsmooth(filterobject=X_filtered, kappa=pe_BSd_3F$fit$par_est$kappa)

#' @export
xsmooth <- function(filterobject, kappa){
  smooth <- RTS_sm_bas(filterobject$X_t, filterobject$X_t_c, filterobject$S_t, filterobject$S_t_c, kappa, (ncol(filterobject$X_t)-1))
  return(smooth)
}

#======================== - Goodness of fit - ===================================

## - Fitted rates
#' @title mubar_hat
#'
#' @description Function returning the fitted values of the average mortality rates
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return Matrix with the fitted `mu_bar` rates given model and parameters
#' @examples
#' fitted_BSd <- mubar_hat(model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar)

#' @export
mubar_hat <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      fitted <- mu_bar_hat_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(fitted)
    } else{
      fitted <- mu_bar_hat_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
      return(fitted)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            fitted <- mu_bar_hat_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          } else{ # - the factors must be three
            fitted <- mu_bar_hat_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          }
        } else{
        fitted <- mu_bar_hat_BSi(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(fitted)
      }
    } else{
      if(model=="CIR"){
        fitted <- mu_bar_hat_CIR(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
        return(fitted)
      } else{
        if(model=="AFGNS"){    # - if none of the previous model was selected, then it is an AFGNS
          if(fact_dep==TRUE){
            fitted <- mu_bar_hat_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          } else{
            fitted <- mu_bar_hat_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
            return(fitted)
          }}else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                fitted <- mu_bar_hat_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(fitted)
              } else{
                fitted <- mu_bar_hat_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                return(fitted)
              }
            }else{
              if(model=="AFRNS"){
                if(fact_dep==TRUE){
                  fitted <- mu_bar_hat_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(fitted)
                } else{
                  fitted <- mu_bar_hat_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                  return(fitted)
                }
              }else{
                if(model=="GMk"){
                  if(fact_dep==TRUE){
                    fitted <- mu_bar_hat_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(fitted)
                  } else{
                    fitted <- mu_bar_hat_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data)
                    return(fitted)
                  }
                }
              }
            }
          }
      }
    }
  }
}

## - Fitted rates
#' @title std_res
#'
#' @description Function returning the standardized residuals using the method of Ungolo et. al. (2023)
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#'
#' @return Matrix with the standardized residuals given model and parameters.
#' @examples
#' std_resid <- std_res(model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar)

#' @export
std_res <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
    } else{
      std_res <- residuals_std_AFNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
      return(std_res)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            std_res <- residuals_std_BSd_2F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          } else{
            std_res <- residuals_std_BSd_3F(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          }
        }else{ # - it is the Blackburn-Sherris model with independent factors
          std_res <- residuals_std_BSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
          return(std_res)
        }
      }
     else{
      if(model=="CIR"){
        std_res <- residuals_std_CIR(data, x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc))
        return(std_res)
      } else{
        if(model=="AFGNS"){
          if(fact_dep==TRUE){
            std_res <- residuals_std_AFGNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          } else{
            std_res <- residuals_std_AFGNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
            return(std_res)
          }}else{
            if(model=="AFRNS"){
              if(fact_dep==TRUE){
                std_res <- residuals_std_AFRNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              } else{
                std_res <- residuals_std_AFRNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
                return(std_res)
              }
            }else{
              if(model=="AFUNS"){
                if(fact_dep==TRUE){
                  std_res <- residuals_std_AFUNSd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
                  return(std_res)
                } else{
                  std_res <- residuals_std_AFUNSi(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
                  return(std_res)
                }
              }else{
                if(model=="GMk"){
                  if(fact_dep==TRUE){
                    std_res <- residuals_std_GMkd(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc))
                    return(std_res)
                  } else{
                    std_res <- residuals_std_GMki(data, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc))
                    return(std_res)
                  }
                }}
            }
          }
      }
    }
  }
}


## - Probability of negative rates
#' @title prob_neg_mu
#'
#' @description Probability of negative rates when projected over `n` years, based on simulated values of the latent variables
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#' @param years_proj Number of years of ahead-projected rates (default 1 year)
#' @param n_simulations Number of simulations (default value 100000)
#'
#' @return A vector of probabilities of negative rates for each age considered in the dataset.
#' @examples
#' neg_mu_p <- prob_neg_mu(model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar, years_proj = 10, n_simulations=1000)

#' @export
prob_neg_mu <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data, years_proj=1, n_simulations=100000){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      pr_neg <- pr_neg_rates_f_AFNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
      return(pr_neg)
    } else{
        pr_neg <- pr_neg_rates_f_AFNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
        return(pr_neg)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            pr_neg <- pr_neg_rates_f_BSd_2F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          } else{
            pr_neg <- pr_neg_rates_f_BSd_3F(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          }
      } else{ # - it is the Blackburn-Sherris model with independent factors
          pr_neg <- pr_neg_rates_f_BSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
          return(pr_neg)
        }
      }
     else{
      if(model=="CIR"){
          pr_neg <- 0
          return(pr_neg)
      } else{
        if(model=="AFGNS"){
        if(fact_dep==TRUE){
            pr_neg <- pr_neg_rates_f_AFGNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
        } else{
            pr_neg <- pr_neg_rates_f_AFGNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
        }}else{
        if(model=="AFRNS"){
          if(fact_dep==TRUE){
            pr_neg <- pr_neg_rates_f_AFRNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          } else{
            pr_neg <- pr_neg_rates_f_AFRNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
            return(pr_neg)
          }
        }else{
          if(model=="AFUNS"){
            if(fact_dep==TRUE){
              pr_neg <- pr_neg_rates_f_AFUNSd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
              return(pr_neg)
            } else{
              pr_neg <- pr_neg_rates_f_AFUNSi(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
              return(pr_neg)
            }
          }else{
            if(model=="GMk"){
              if(fact_dep==TRUE){
                pr_neg <- pr_neg_rates_f_GMkd(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
                return(pr_neg)
              } else{
                pr_neg <- pr_neg_rates_f_GMki(n_sim=n_simulations, x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, yrs_proj=years_proj)
                return(pr_neg)
              }
            }
          }
        }
      }
    }
  }
}}



#======================== - Projection - ===================================

# Improve output presentation with a smarter use of print
#' @title affine_project
#'
#' @description Projected survival curves and average mortality rates over n years ahead projections
#'
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (`fact_dep=TRUE`) or independence
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters. If not set, then default values will be used
#' @param data Table with the average mortality rates
#' @param years_proj Number of years of ahead-projected rates
#'
#' @return A vector of the projected survival rates for each age
#' @examples
#' BSd_3F_proj <- affine_project(model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar, years_proj = 1)

#' @export
affine_project <- function(model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data, years_proj=1){
  if(model=="AFNS"){
    if(fact_dep==TRUE){
      projection <- S_t_AFNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
      return(projection)
    } else{
      projection <- S_t_AFNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
      return(projection)
    }
  } else{
    if(model=="BS"){
      if(fact_dep==TRUE){
          if(n_factors==2){
            projection <- S_t_BSd_2F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          } else{
            projection <- S_t_BSd_3F_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          }
        }else{
        projection <- S_t_BSi_proj(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma = parameters$sigma[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
        return(projection)
      }} else{
      if(model=="CIR"){
        projection <- S_t_CIR_proj(x0=parameters$x0[1:n_factors], delta=parameters$delta[1:n_factors], kappa=parameters$kappa[1:n_factors], sigma=parameters$sigma[1:n_factors], theta_P=parameters$theta_P[1:n_factors], r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
        return(projection)
      } else{
        if(model=="AFGNS"){
          # - if none of the previous model was selected, then it is an AFGNS
          if(fact_dep==TRUE){
            projection <- S_t_AFGNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          } else{
            projection <- S_t_AFGNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
            return(projection)
          }
        }else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              projection <- S_t_AFRNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
              return(projection)
            } else{
              projection <- S_t_AFRNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
              return(projection)
            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                projection <- S_t_AFUNSd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                return(projection)
              } else{
                projection <- S_t_AFUNSi_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                return(projection)
              }
            }else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  projection <- S_t_GMkd_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                  return(projection)
                } else{
                  projection <- S_t_GMki_proj(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), data, proj_years=years_proj)
                  return(projection)
                }
              }
            }
          }
        }
      }
    }
  }}



#======================= - Graphics - ======================================

## - Heatmap of residuals
# Improve output presentation with a smarter use of print
#' @title heatmap_res
#'
#' @description Returns the heatmaps created by using heatmaply
#'
#' @param residuals Table of residuals obtainable using std.res
#' @param color TRUE if colored (default) or FALSE if black and white
#'
#' @return Plot of the residuals heatmap
#' @examples
#' std_resid <- std_res(model="AFNS", fact_dep=FALSE, parameters=pe_AFNSi$fit$par_est, data=mu_bar) # - Get matrix of standardized residuals
#' heatmap_res(residuals=std_resid, color=FALSE)

#' @export
heatmap_res <- function(residuals, color=TRUE){
  if(color==FALSE){
    heatmaply::heatmaply(residuals[c(nrow(residuals):1),], dendrogram = "none", xlab = "Year", ylab="Age", dynamicTicks=FALSE, hide_colorbar=FALSE, margins= c(0,0,0,0), color = c("white", "black"))
  } else{
    heatmaply::heatmaply(residuals[c(nrow(residuals):1),], dendrogram = "none", xlab = "Year", ylab="Age", dynamicTicks=FALSE, hide_colorbar=FALSE, margins= c(0,0,0,0))
  }
}

#======================= - Parameter uncertainty - ======================================

## - Fix by adding a print function for partial results
#' @title par_cov
#'
#' @description Estimation of the uncertainty of the parameters by Bootstrap or Multiple imputation
#'
#' @param method MI if Multiple Imputation or Bootstrap if the Bootstrap is desired
#' @param model Affine model. E.g. for the Blackburn-Sherris model we have model="BS", and so on
#' @param fact_dep Boolean parameter indicating whether estimate models with factor dependence (fact_dep=TRUE) or independence (fact_dep=FALSE)
#' @param n_factors Number of factors. For some models, these are set by default
#' @param parameters Starting value for the parameters.
#' @param data Table with the average mortality rates
#' @param D_se Number of Iterations if Multiple Imputation is used (default set to 50)
#' @param BS_s Number of Boostrap samples if this method is used
#' @param t_excl Number of time-periods to be excluded for stability of the Bootstrap method
#' @param max_iter Maximum number of iterations for the Coordinate Ascent estimation algorithm in absence of convergence
#' @param tolerance Minimum value of improvement of the log-likelihood value before convergence
#' @param wd Working directory for saving the partial output of the estimation process
#'
#' @return A list with the covariance matrix of the parameter estimates `Cov_par_est` and their standard errors `St_err`
#' @examples
#' par_unc_MI <- par_cov(method="MI", model="BS", fact_dep=TRUE, n_factors=3, parameters=pe_BSd_3F$fit$par_est, data=mu_bar, D_se=5, max_iter=10, tolerance=0.1)

#' @export
par_cov <- function(method=c("MI", 'Bootstrap'), model=c("BS", "AFNS", "AFGNS", "AFUNS", "AFRNS", "CIR", "GMk"), fact_dep=c(FALSE, TRUE), n_factors=3, parameters, data, D_se=50, BS_s=500, t_excl=4, max_iter=200, tolerance=0.1, wd=0){

  if(method=="MI"){
    if(model=="AFNS"){
      if(fact_dep==TRUE){
          covariance <- CovEst_MI_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      } else{
          covariance <- CovEst_MI_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)
      }
    } else{
      if(model=="BS"){
        if(fact_dep==TRUE){
            if(n_factors==2){
              covariance <- CovEst_MI_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            } else{ # - the factors must be three
              covariance <- CovEst_MI_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
           }

        } else{ # - it is the Blackburn-Sherris model with independent factors

            covariance <- CovEst_MI_BSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        }
      } else{
        if(model=="CIR"){
            covariance <- CovEst_MI_CIR(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, theta_P = parameters$theta_P, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        } else{
          if(model=="AFGNS"){   # - if none of the previous model was selected, then it is an AFGNS
          if(fact_dep==TRUE){
              covariance <- CovEst_MI_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          } else{
              covariance <- CovEst_MI_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          }
          }else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              covariance <- CovEst_MI_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            } else{
              covariance <- CovEst_MI_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                covariance <- CovEst_MI_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              } else{
                covariance <- CovEst_MI_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              }
            }else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  covariance <- CovEst_MI_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                } else{
                  covariance <- CovEst_MI_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, D_se=D_se, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                }
              }
            }
          }
        }
      }
    }}
    }else{ #### - Perform Bootstrap
    if(model=="AFNS"){
      if(fact_dep==TRUE){
          covariance <- CovEst_BS_AFNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      } else{
          covariance <- CovEst_BS_AFNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
          return(covariance)

      }
    } else{
      if(model=="BS"){
        if(fact_dep==TRUE){
            if(n_factors==2){
              covariance <- CovEst_BS_BSd_2F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            } else{ # - the factors must be three
              covariance <- CovEst_BS_BSd_3F(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            }

        } else{ # - it is the Blackburn-Sherris model with independent factors
            covariance <- CovEst_BS_BSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        }
      } else{
        if(model=="CIR"){
            covariance <- CovEst_BS_CIR(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, theta_P = parameters$theta_P, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
            return(covariance)

        } else{
          if(model=="AFGNS"){
          if(fact_dep==TRUE){
              covariance <- CovEst_BS_AFGNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

          } else{
              covariance <- CovEst_BS_AFGNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
          }
          }else{
          if(model=="AFRNS"){
            if(fact_dep==TRUE){
              covariance <- CovEst_BS_AFRNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)

            } else{
              covariance <- CovEst_BS_AFRNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
              return(covariance)
            }
          }else{
            if(model=="AFUNS"){
              if(fact_dep==TRUE){
                covariance <- CovEst_BS_AFUNSd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)

              } else{
                covariance <- CovEst_BS_AFUNSi(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                return(covariance)
              }
            }else{
              if(model=="GMk"){
                if(fact_dep==TRUE){
                  covariance <- CovEst_BS_GMkd(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma_dg = parameters$sigma_dg, Sigma_cov = parameters$Sigma_cov, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)

                } else{
                  covariance <- CovEst_BS_GMki(x0=parameters$x0, delta=parameters$delta, kappa=parameters$kappa, gamma=parameters$gamma, sigma = parameters$sigma, r=c(parameters$r1, parameters$r2, parameters$rc), mu_bar=data, n_BS=BS_s, t_ex = t_excl, max_it=max_iter, tolerance_lev=tolerance, workdir=wd)
                  return(covariance)
                }
              }
            }
          }
        }
      }
    }
  }}
}


#========= - Other utility functions - ============

#==================================================
# - Root Mean Squared Error (RMSE)
#==================================================

# - First reproduce mu_bar_hat for all models
#' @title RMSE
#'
#' @description Calculation of the Root Mean Squared Error
#'
#' @param observed the table of observed average mortality rates
#' @param estimated table of the fitted average mortality rates as obtainable through the function mubar_hat
#'
#' @return Returns a scalar denoting the Root Mean Squared Error from the fitted model
#' @examples
#' RMSE(mu_bar, fitted_BSd)

#' @export
RMSE <- function(observed, estimated){
  RMSE <- sqrt(mean((observed - estimated)^2))
  return(RMSE)
}


#=======================================================
# - Mean Absolute Percentage Error (MAPE)
#=======================================================
# - First reproduce mu_bar_hat for all models
#' @title MAPE_age
#'
#' @description Calculation of the Mean Absolute Percentage Error by age (i.e. by each row of the average mortality table)
#'
#' @param observed the table of observed average mortality rates
#' @param estimated table of the fitted average mortality rates as obtainable through the function mubar_hat
#'
#' @return Returns a vector with the Mean Absolute Percentage Error by age
#' @examples
#' MAPE_age(mu_bar, fitted_BSd)

#' @export
MAPE_age <- function(observed, estimated){
  MAPE <- rowMeans(abs((observed - estimated)/observed))
  return(MAPE)
}

## - 0/1 residuals (0 if negative, 1 if positive)
#' @title residuals_01
#'
#' @description Computation of the positive (1) and negative (0) residuals. Useful to understand the presence of the patterns within the residuals
#'
#' @param observed the table of observed average mortality rates
#' @param estimated table of the fitted average mortality rates as obtainable through the function mubar_hat
#'
#' @return Returns a matrix with 0-1 residuals by age and year
#' @examples
#' residuals_01(mu_bar, fitted_BSd)

#' @export
residuals_01 <- function(observed, estimated){
  residuals <- ifelse(observed - estimated > 0, 1, 0)
  colnames(residuals) <- colnames(observed)
  rownames(residuals) <- rownames(observed)
  return(residuals)
}

#' @title avg2rates
#'
#' @description Converts a table with average mortality rates into mortality rates
#'
#' @param mu_bar the table of average mortality rates
#'
#' @return Returns a matrix of the same dimension as mu_bar with mortality rates
#' @examples
#' mu_rates <- avg2rates(mu_bar)

#' @export
avg2rates <- function(mu_bar){
  mu <- mu_bar
  for(i in 2:nrow(mu_bar)){
    mu[i,] <- i * mu_bar[i,] - (i-1) * mu_bar[i-1,]
  }
  return(mu)
}

#' @title rates2avg
#'
#' @description Converts a table with mortality rates into average mortality rates for use in estimation of affine mortality models
#'
#' @param mu the table of average mortality rates
#'
#' @return Returns a matrix of the same dimension as mu with average mortality rates
#' @examples
#' mu_bar <- rates2avg(mu_rates)

#' @export
rates2avg <- function(mu){
  mu_bar <- mu
  for(row in 1:nrows(mu)){
    for(col in 1:ncol(mu)){
      mu_bar[row,col] <- mean(mu[1:row,col])
    }
  }
  return(mu_bar)
}





