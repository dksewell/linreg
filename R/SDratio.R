#' @name BF
#' @aliases SDratio
#'
#' 
#' @title Bayes factors for lm_b and glm_b
#' 
#' @description
#' Bayes factors for Bayesian regression objects using the Savage-Dickey ratio
#' 
#' 
#' @details
#' 
#' Bayes factors are given in terms of favoring the two-tailed alternative hypothesis 
#' vs. the null hypothesis that the regression coefficient equals zero.  Currently only 
#' implemented for lm_b objects.  Use bridgesampling package for aov_b.
#' 
#' Interpretation is taken from Kass and Raftery.
#' 
#' @param object lm_b or glm_b object
#' 
#' @references 
#' 
#' Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association, 90(430), 773–795.
#' 
#' @importFrom ks kde
#' 
#' @export

#' @export
SDratio = function(object){
  UseMethod("SDratio")
}


#' @rdname IC
#' @export
SDratio.lm_b = function(object){
  if(is.na(object$hyperparameters[1])) stop("Cannot compute Bayes factors with an improper prior.")
  
  log_numerators = 
    dlst(0.0,
         df = object$hyperparameters$a,
         mu = object$hyperparameters$mu,
         sigma = 
           sqrt(object$hyperparameters$b / 
                  object$hyperparameters$a * 
                  diag(qr.solve(object$hyperparameters$V))),
         log = TRUE)
  
  log_denominators = 
    dlst(0.0,
         df = object$posterior_parameters$a,
         mu = object$posterior_parameters$mu,
         sigma = 
           sqrt(object$posterior_parameters$b / 
                  object$posterior_parameters$a * 
                  diag(qr.solve(object$posterior_parameters$V))),
         log = TRUE)
  
  results = 
    tibble(Variable = colnames(object$posterior_parameters$V_tilde),
           `BF favoring alternative` = 
             exp(log_numerators - log_denominators)) |>
    mutate(bf_max = pmax(`BF favoring alternative`,
                         1.0 / `BF favoring alternative`)) |> 
    mutate(Interpretation = 
             ifelse(bf_max <= 3.2,
                    "Not worth more than a bare mention",
                    ifelse(bf_max <= 10,
                           "Substantial",
                           ifelse(bf_max <= 100,
                                  "Strong",
                                  "Decisive"))),
           .keep = "unused")
  
  return(results)
}




#' @rdname BF
#' @export
SDratio.glm_b = function(object){
  if(is.na(object$hyperparameters[1])) stop("Cannot compute Bayes factors with an improper prior.")
  
  # Evaluate prior density at zero
  prior_cov = NULL
  try({
    prior_cov = 
      chol2inv(chol(object$hyperparameters$prior_beta_precision))
  },silent = TRUE)
  if(is.null(prior_cov)){
    try({
      prior_cov = 
        qr.solve(object$hyperparameters$prior_beta_precision)
    },silent = TRUE)
  }
  if(is.null(prior_cov)){
    try({
      prior_cov = 
        solve(object$hyperparameters$prior_beta_precision)
    },silent = TRUE)
  }
  
  log_numerators = 
    dnorm(0.0,
          object$hyperparameters$prior_beta_mean,
          sqrt(diag(prior_cov)),
          log = TRUE)
  
  # Evaluate posterior density at zero
  if("posterior_covariance" %in% names(object)){
    log_denominators = 
      dnorm(0.0,
            object$summary$`Post Mean`,
            sqrt(diag(object$posterior_covariance)),
            log = TRUE)
  }else{
    kde_helper = function(i){
      if( (0 < min(object$proposal_draws[,i])) |
          (0 > max(object$proposal_draws[,i]))){
        return(.Machine$double.eps)
      }else{
        kde_fit = 
          ks::kde(object$proposal_draws[,i],
                  w = object$importance_sampling_weights * NROW(object$importance_sampling_weights))
        
        return(kde_fit$estimate[which.min(abs(kde_fit$eval.points - 0.0))])
      }
    }
    
    log_denominators = 
      sapply(1:ncol(object$proposal_draws), kde_helper) |> 
      pmax(.Machine$double.eps) |> 
      log()
  }
  
  
  results = 
    tibble(Variable = object$summary$Variable,
           `BF favoring alternative` = 
             exp(log_numerators - log_denominators)) |>
    mutate(bf_max = pmax(`BF favoring alternative`,
                         1.0 / `BF favoring alternative`)) |> 
    mutate(Interpretation = 
             ifelse(bf_max <= 3.2,
                    "Not worth more than a bare mention",
                    ifelse(bf_max <= 10,
                           "Substantial",
                           ifelse(bf_max <= 100,
                                  "Strong",
                                  "Decisive"))),
           .keep = "unused")
  
  return(results)
}
