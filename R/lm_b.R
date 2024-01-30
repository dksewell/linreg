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
#' @param prior_var_shape numeric. Shape parameter for the inverse gamma prior on
#' the residual variance.
#' @param prior_var_rate  numeric. Rate parameter for the inverse gamma prior on
#' the residual variance.
#' @param CI_level numeric. Credible interval level.
#' 
#' @return lm_b returns an object of class "lm_b", which behaves as a list with 
#' the following elements:
#' \itemize{
#'  \item summary - data.frame giving results for regression coefficients
#'  \item post_parms - list giving the posterior parameters
#'  \item hyperparms - list giving the user input (or default) hyperparameters used
#'  \item fitted - posterior mean of the individuals' means
#'  \item residuals - posterior mean of the residuals
#'  \item formula, data - input by user
#' }
#' 
#' @import magrittr
#' @import dplyr
#' @import coda
#' @import extraDistr
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
  
  prior = 
    c("zellner","conjugate","improper")[pmatch(tolower(prior),
                                               c("zellner","conjugate","improper"),duplicates.ok = FALSE)]
  
  if(prior == "zellner"){
    if(missing(zellner_g)){
      warning("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
      zellner_g = N
    }
    if(missing(prior_var_shape)){
      warning("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      warning("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
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
      data.frame(Variable = rownames(mu_tilde),
                 `Post Mean` = mu_tilde,
                 Lower = qlst(alpha/2,
                              a_tilde,
                              mu_tilde,
                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 Upper = qlst(1 - alpha/2,
                              a_tilde,
                              mu_tilde,
                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 ROPE = 
                   plst(0.05 * s_y / c(1,s_j),
                        a_tilde,
                        mu_tilde,
                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
                   plst(-0.05 * s_y / c(1,s_j),
                        a_tilde,
                        mu_tilde,
                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 `Prob Dir` = 
                   sapply(plst(0,
                               a_tilde,
                               mu_tilde,
                               sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                          function(x) max(x,1-x))
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           post_parms = list(mu_tilde = mu_tilde,
                             V_tilde = V_tilde,
                             a_tilde = a_tilde,
                             b_tilde = b_tilde),
           hyperparms = list(mu = rep(0,p),
                             V = V,
                             a = a,
                             b = b))
  }
  
  
  if(prior == "conjugate"){
    
    if(missing(prior_beta_mean)){
      warning("\nThe mu hyperparameter in the normal prior is not specified.  It will be set automatically to 0.\n")
      mu = rep(0,p)
    }else{
      mu = prior_beta_mean
    }
    if(missing(prior_var_shape)){
      warning("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      warning("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      b = 0.001
    }else{
      b = prior_var_rate
    }
    if(missing(prior_beta_precision)){
      warning("\nThe V hyperparameter in the normal prior is not specified.  It will be set automatically to 4/25Diag(s^2_{X_j})")
      V = 4/25 * diag(c(25/4*0.01,s_j^2))
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
      data.frame(Variable = rownames(mu_tilde),
                 `Post Mean` = mu_tilde,
                 Lower = qlst(alpha/2,
                              a_tilde,
                              mu_tilde,
                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 Upper = qlst(1 - alpha/2,
                              a_tilde,
                              mu_tilde,
                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 ROPE = 
                   plst(0.05 * s_y / c(1,s_j),
                        a_tilde,
                        mu_tilde,
                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
                   plst(-0.05 * s_y / c(1,s_j),
                        a_tilde,
                        mu_tilde,
                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                 `Prob Dir` = 
                   sapply(plst(0,
                               a_tilde,
                               mu_tilde,
                               sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                          function(x) max(x,1-x))
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           post_parms = list(mu_tilde = mu_tilde,
                             V_tilde = V_tilde,
                             a_tilde = a_tilde,
                             b_tilde = b_tilde),
           hyperparms = list(mu = mu,
                             V = V,
                             a = a,
                             b = b))
  }
  
  
  if(prior == "improper"){
    
    mod = lm(formula,data)
    mu_tilde = coef(mod)
    Sigma = sigma(mod)^2 * XtX_inv
    results = 
      data.frame(Variable = names(mu_tilde),
                 `Post Mean` = mu_tilde,
                 Lower = qlst(alpha/2,
                              N - p,
                              mu_tilde,
                              sqrt(diag(Sigma))),
                 Upper = qlst(1 - alpha/2,
                              N - p,
                              mu_tilde,
                              sqrt(diag(Sigma))),
                 ROPE = 
                   plst(0.05 * s_y / c(1,s_j),
                        N - p,
                        mu_tilde,
                        sqrt(diag(Sigma))) -
                   plst(-0.05 * s_y / c(1,s_j),
                        N - p,
                        mu_tilde,
                        sqrt(diag(Sigma))),
                 `Prob Dir` = 
                   sapply(plst(0,
                               N - p,
                               mu_tilde,
                               sqrt(diag(Sigma))),
                          function(x) max(x,1-x))
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           post_parms = list(mu_tilde = mu_tilde,
                             V_tilde = XtX,
                             Sigma = Sigma,
                             a_tilde = N - p,
                             b_tilde = sum(resid(mod)^2)),
           hyperparms = NA)
  }
  
  
  return_object$fitted = 
    X %*% return_object$summary$Post.Mean
  return_object$residuals = 
    y - return_object$fitted
  
  return_object$formula = formula
  return_object$data = data
  
  rownames(return_object$summary) = NULL
  
  class(return_object) = "lm_b"
  return(return_object)
}