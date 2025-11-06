#' Bayes factors 
#' 
#' Bayes factors for Bayesian regression objects using the Savage-Dickey ratio
#' 
#' Bayes factors are given in terms of favoring the two-tailed alternative hypothesis 
#' vs. the null hypothesis that the regression coefficient equals zero.  Currently only 
#' implemented for lm_b objects.  Use bridgesampling package for aov_b.
#' 
#' @param object lm_b object
#' 
#' @export



SDratio = function(object){
  if(is.na(object$hyperparms[1])) stop("Cannot compute Bayes factors with an improper prior.")
  
  log_numerators = 
    dlst(0.0,
         df = object$hyperparms$a,
         mu = object$hyperparms$mu,
         sigma = 
           sqrt(object$hyperparms$b / 
                  object$hyperparms$a * 
                  diag(qr.solve(object$hyperparms$V))),
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
             exp(log_numerators - log_denominators))
  
  return(results)
}

