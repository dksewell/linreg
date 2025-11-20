#' Bayesian Linear Models
#' 
#' lm_b is used to fit linear models.  It can be used to carry out 
#' regression, single stratum analysis of variance and analysis of 
#' covariance (although aov_b may provide a more convenient 
#' interface for ANOVA.)
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.
#' @param prior character.  One of "zellner", "conjugate", or "improper", giving 
#' the type of prior used on the regression coefficients.
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  Ignored unless prior = "zellner".
#' @param prior_beta_mean numeric vector of same length as regression coefficients 
#' (denoted p). Unless otherwise specified, automatically set to rep(0,p).  Ignored 
#' unless prior = "conjugate".
#' @param prior_beta_precision pxp matrix giving a priori precision matrix to be 
#' scaled by the residual precision.
#' @param prior_var_shape numeric. Twice the shape parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param prior_var_rate  numeric. Twice the rate parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param ROPE vector of positive values giving ROPE boundaries for each regression 
#' coefficient.  Optionally, you can not include a ROPE boundary for the intercept. 
#' If missing, defaults go to those suggested by Kruchke (2018).
#' @param CI_level numeric. Credible interval level.
#' 
#' @return lm_b() returns an object of class "lm_b", which behaves as a list with 
#' the following elements:
#' \itemize{
#'  \item summary - tibble giving results for regression coefficients
#'  \item posterior_parameters - list giving the posterior parameters
#'  \item hyperparameters - list giving the user input (or default) hyperparameters used
#'  \item fitted - posterior mean of the individuals' means
#'  \item residuals - posterior mean of the residuals
#'  \item formula, data - input by user
#' }
#' 
#' 
#' \strong{ROPE:}
#' 
#' If missing, the ROPE bounds will be given under the principle of "half of a 
#' small effect size."  Using Cohen's D of 0.2 as a small effect size, the ROPE is 
#'  built under the principle that moving the full range of X (i.e., \eqn{\pm 2 s_x}) 
#'  will not move the mean of y by more than the overall mean of \eqn{y} 
#'  minus \eqn{0.1s_y} to the overall mean of \eqn{y} plus \eqn{0.1s_y}.  
#'  The result is a ROPE equal to \eqn{|\beta_j| < 0.05s_y/s_j}.  If the covariate is 
#'  binary, then this is simply \eqn{|\beta_j| < 0.2s_y}.
#'  
#' 
#' @import dplyr
#' @import coda
#' @import extraDistr
#' @import tibble
#' @export
#' @exportClass lm_b


lm_b = function(formula,
                data,
                prior = c("zellner","conjugate","improper")[1],
                zellner_g,
                prior_beta_mean,
                prior_beta_precision,
                prior_var_shape = 0.001,
                prior_var_rate = 0.001,
                ROPE,
                CI_level = 0.95){
  
  # Extract 
  y = model.response(model.frame(formula, data))
  X = model.matrix(formula,data)
  p = ncol(X)
  alpha = 1 - CI_level
  XtX = crossprod(X)
  XtX_inv = qr.solve(XtX)
  s_y = sd(y)
  s_j = apply(X[,-1,drop = FALSE],2,sd)
  N = nrow(X)
  
  # Get ROPE 
  if(missing(ROPE)){
    
    ROPE = 
      c(NA,
        0.2 * s_y / ifelse(apply(X[,-1],2,
                                 function(z) isTRUE(all.equal(0:1,
                                                              sort(unique(z))))),
                           1.0,
                           4.0 * s_j))
    # From Kruchke (2018) on standardized regression
    # This considers a small change in y to be \pm 0.1s_y (half of Cohen's D small effect).
    # So this is the change due to moving through the range of x (\pm 2s_X).
    # For binary (or one-hot) use 1.
    
    
  }else{
    if( !(length(ROPE) %in% (ncol(X) - 1:0))) stop("Length of ROPE, if supplied, must match the number of regression coefficients (or one less, if no ROPE is for the intercept).")
    if( any(ROPE) <= 0 ) stop("User supplied ROPE values must be positive.  The ROPE is assumed to be +/- the user supplied values.")
    if(length(ROPE) == ncol(X) - 1) ROPE = c(NA,ROPE)
  }
  ROPE_bounds = 
    paste("(",-round(ROPE,3),",",round(ROPE,3),")",sep="")
  
  
  
  
  prior = 
    c("zellner","conjugate","improper")[pmatch(tolower(prior),
                                               c("zellner","conjugate","improper"),duplicates.ok = FALSE)]
  
  if(prior == "zellner"){
    if(missing(zellner_g)){
      message("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
      zellner_g = N
    }
    if(missing(prior_var_shape)){
      message("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      message("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      b = 0.001
    }else{
      b = prior_var_rate
    }
    
    V = 1/zellner_g * XtX
    
    V_tilde = (zellner_g + 1.0) / zellner_g * XtX
    V_tilde_inv = qr.solve(V_tilde)
    mu_tilde = V_tilde_inv %*% crossprod(X,y)
    a_tilde = a + N
    b_tilde = b + crossprod(y) - crossprod(mu_tilde, V_tilde %*% mu_tilde)
    b_tilde = drop(b_tilde)
    
    results = 
      tibble(Variable = rownames(mu_tilde),
             `Post Mean` = mu_tilde[,1],
             Lower = qlst(alpha/2,
                          a_tilde,
                          mu_tilde,
                          sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             Upper = qlst(1 - alpha/2,
                          a_tilde,
                          mu_tilde,
                          sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             `Prob Dir` = 
               sapply(plst(0,
                           a_tilde,
                           mu_tilde,
                           sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                      function(x) max(x,1-x)),
             ROPE = 
               plst(ROPE,
                    a_tilde,
                    mu_tilde,
                    sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
               plst(-ROPE,
                    a_tilde,
                    mu_tilde,
                    sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                             V_tilde = V_tilde,
                             a_tilde = a_tilde,
                             b_tilde = b_tilde),
           hyperparameters = list(mu = rep(0,p),
                             V = V,
                             a = a,
                             b = b))
  }
  
  
  if(prior == "conjugate"){
    
    if(missing(prior_beta_mean)){
      message("\nThe mu hyperparameter in the normal prior is not specified.  It will be set automatically to 0.\n")
      mu = rep(0,p)
    }else{
      mu = prior_beta_mean
    }
    if(missing(prior_var_shape)){
      message("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      message("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      b = 0.001
    }else{
      b = prior_var_rate
    }
    if(missing(prior_beta_precision)){
      message("\nThe V hyperparameter in the normal prior is not specified.  It will be set automatically to 4/25Diag(s^2_{X_j})")
      V = 
        diag(1.0 / c(2.5 * s_y,2.5 * s_y / s_j)^2)
    }else{
      V = prior_beta_precision
    }
    
    V_tilde = V + XtX
    V_tilde_inv = qr.solve(V_tilde)
    mu_tilde = V_tilde_inv %*% (V%*%mu + crossprod(X,y))
    a_tilde = a + N
    b_tilde = b + crossprod(y) - crossprod(mu_tilde, V_tilde %*% mu_tilde)
    b_tilde = drop(b_tilde)
    
    results = 
      tibble(Variable = rownames(mu_tilde),
             `Post Mean` = mu_tilde[,1],
             Lower = qlst(alpha/2,
                          a_tilde,
                          mu_tilde,
                          sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             Upper = qlst(1 - alpha/2,
                          a_tilde,
                          mu_tilde,
                          sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             `Prob Dir` = 
               sapply(plst(0,
                           a_tilde,
                           mu_tilde,
                           sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                      function(x) max(x,1-x)),
             ROPE = 
               plst(ROPE,
                    a_tilde,
                    mu_tilde,
                    sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
               plst(-ROPE,
                    a_tilde,
                    mu_tilde,
                    sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
             `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                             V_tilde = V_tilde,
                             a_tilde = a_tilde,
                             b_tilde = b_tilde),
           hyperparameters = list(mu = mu,
                             V = V,
                             a = a,
                             b = b))
  }
  
  
  if(prior == "improper"){
    
    mod = lm(formula,data)
    mu_tilde = coef(mod)
    Sigma = sigma(mod)^2 * XtX_inv
    results = 
      tibble(Variable = names(mu_tilde),
             `Post Mean` = mu_tilde,
             Lower = qlst(alpha/2,
                          N - p,
                          mu_tilde,
                          sqrt(diag(Sigma))),
             Upper = qlst(1 - alpha/2,
                          N - p,
                          mu_tilde,
                          sqrt(diag(Sigma))),
             `Prob Dir` = 
               sapply(plst(0,
                           N - p,
                           mu_tilde,
                           sqrt(diag(Sigma))),
                      function(x) max(x,1-x)),
             ROPE = 
               plst(ROPE,
                    N - p,
                    mu_tilde,
                    sqrt(diag(Sigma))) -
               plst(-ROPE,
                    N - p,
                    mu_tilde,
                    sqrt(diag(Sigma))),
             `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                             V_tilde = XtX,
                             Sigma = Sigma,
                             a_tilde = N - p,
                             b_tilde = sum(resid(mod)^2)),
           hyperparameters = NA)
  }
  
  
  return_object$fitted = 
    drop(X %*% return_object$summary$`Post Mean`)
  return_object$residuals = 
    drop(y - return_object$fitted)
  
  return_object$formula = formula
  return_object$data = data
  return_object$prior = prior
  return_object$ROPE = ROPE
  
  rownames(return_object$summary) = NULL
  
  class(return_object) = "lm_b"
  return(return_object)
}