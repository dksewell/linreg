#' Bayesian Generalized Linear Models
#' 
#' glm_b is used to fit linear models.  It can be used to carry out 
#' regression, single stratum analysis of variance and analysis of 
#' covariance (although aov_b may provide a more convenient 
#' interface for ANOVA.)
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.
#' @param family A description of the error distribution and link function 
#' to be used in the model. See ?glm for more information.
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param prior character.  One of "zellner", "conjugate", or "improper", giving 
#' the type of prior used on the regression coefficients.
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  Ignored unless prior = "zellner". Default is the number of observations.
#' @param prior_beta_mean numeric vector of same length as regression coefficients 
#' (denoted p). Unless otherwise specified, automatically set to rep(0,p).  Ignored 
#' unless prior = "conjugate".
#' @param prior_beta_precision pxp matrix giving a priori precision matrix to be 
#' scaled by the residual precision.
#' @param ROPE vector of positive values giving ROPE boundaries for each regression 
#' coefficient.  Optionally, you can not include a ROPE boundary for the intercept. 
#' If missing, defaults go to those suggested by Kruchke (2018).
#' @param n_draws Number of draws for importance sampling.  If set to be NA, large sample 
#' normal approximation will be used.
#' @param proposal_df degrees of freedom used in the multivariate t proposal distribution.
#' @param CI_level numeric. Credible interval level.
#' @param seed integer.  Always set your seed!!!
#' 
#' @return glm_b() returns an object of class "glm_b", which behaves as a list with 
#' the following elements:
#' \itemize{
#'  \item summary - tibble giving results for regression coefficients
#'  \item posterior_draws - list giving the posterior parameters
#'  \item hyperparms - list giving the user input (or default) hyperparameters used
#'  \item fitted - posterior mean of the individuals' means
#'  \item residuals - posterior mean of the residuals
#'  \item formula, data - input by user
#' }
#' 
#' \strong{Importance sampling:}
#' 
#' \code{glm_b} will, unless n_draws is set to be NA, perform importance sampling. 
#' The proposal will use a multivariate t distribution, centered at the 
#' posterior mode, with the negative hessian as its precision matrix.  Do NOT
#' treat the proposal_draws as posterior draws.
#' 
#' \strong{Priors:}
#' 
#' If the prior is set to be either "zellner" or "normal", a normal distribution 
#' will be used as the prior of \eqn{\beta}, specifically 
#' \deqn{\beta \sim N(\mu, V)}
#' where \eqn{\mu} is the prior_beta_mean and V is the prior_beta_precision (not covariance) matrix. 
#' \itemize{
#'  \item{\code{zellner}: \code{glm_b} sets \eqn{\mu=0} and \eqn{V = \frac{1}{g} X^{\top} X}.}
#'  \item{\code{normal}: If missing, \code{glm_b} sets \eqn{\mu=0}, and V will be a diagonal matrix.  The first element, 
#'  corresponding to the intercept, will be \eqn{(2.5\times \max{s_y,1})^{-2}}, where 
#'  \eqn{s_y} is the standard deviation of \eqn{y}.  Remaining diagonal elements 
#'  will equal \eqn{(2.5 s_y/s_x)^{-2}}, where \eqn{s_x} is the standard deviations 
#'  of the covariates.  This equates to being 95% certain a priori that a change in 
#'  x by one standard deviation (\eqn{s_x}) would not lead to a change in the linear predictor of 
#'  more than 5 standard deviations (\eqn{5s_y}).  This imposes weak regularization that adapts to the scale 
#'  of the data elements.}
#' }
#' 
#' \strong{ROPE:}
#' 
#' If missing, the ROPE bounds will be given under the principle of "half of a 
#' small effect size."  
#' \itemize{
#'  \item{Gaussian.  Using Cohen's D of 0.2 as a small effect size, the ROPE is 
#'  built under the principle that moving the full range of X (i.e., \eqn{\pm 2 s_x}) 
#'  will not move the mean of y by more than the overall mean of \eqn{y} 
#'  minus \eqn{0.1s_y} to the overall mean of \eqn{y} plus \eqn{0.1s_y}.  
#'  The result is a ROPE equal to \eqn{|\beta_j| < 0.05s_y/s_j}.  If the covariate is 
#'  binary, then this is simply \eqn{|\beta_j| < 0.2s_y}.}
#'  \item{Poisson. FDA guidance suggests a small effect is a rate ratio less 
#'  than 1.25.  We use half this effect: 1.125, and consider ROPE to indicate 
#'  that a moving the full range of X (\eqn{\pm 2s_x} will not change the rate 
#'  ratio by more than this amount.  Thus the ROPE for the regression 
#'  coefficient equals \eqn{|\beta| < \frac{\log(1.125)}{4s_x}}. For binary 
#'  covariates, this is simply \eqn{|\beta| < \log(1.125)}.}
#' }
#' 
#' 
#' 
#' @references 
#' 
#' Kruschke JK. Rejecting or Accepting Parameter Values in Bayesian Estimation. Advances in Methods and Practices in Psychological Science. 2018;1(2):270-280. doi:10.1177/2515245918771304
#' 
#' @import magrittr
#' @import dplyr
#' @import coda
#' @import extraDistr
#' @import tibble
#' @import mvtnorm
#' @export
#' @exportClass glm_b


glm_b = function(formula,
                data,
                family,
                trials,
                prior = c("zellner","normal","improper")[1],
                zellner_g,
                prior_beta_mean,
                prior_beta_precision,
                ROPE,
                n_draws = 5e3,
                proposal_df = 5,
                CI_level = 0.95,
                seed = 1){
  
  set.seed(seed)
  
 # Get correct family
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  # Extract
  mframe = model.frame(formula, data)
  y = model.response(mframe)
  X = model.matrix(formula,data)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  s_j = apply(X[,-1,drop = FALSE],2,sd)
  s_y = sd(y)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  
  # Get ROPE 
  if(missing(ROPE)){
    
    if(family$family == "gaussian"){
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
    }
    if(family$family %in% c("poisson","binomial")){
      ROPE = 
        c(NA,
          log(1.0 + 0.25/2) / ifelse(apply(X[,-1],2,
                                   function(z) isTRUE(all.equal(0:1,
                                                                sort(unique(z))))),
                                   1.0,
                                   4 * s_j))
      # From Kruchke (2018) on rate ratios from FDA <1.25.  
      # Use the same thing for odds ratios.
      # So this is the change due to moving through the range of x (\pm 2s_X).
      # For binary (or one-hot) use 1.
    }
    
  }else{
    if( !(length(ROPE) %in% (ncol(X) - 1:0))) stop("Length of ROPE, if supplied, must match the number of regression coefficients (or one less, if no ROPE is for the intercept).")
    if( any(ROPE) <= 0 ) stop("User supplied ROPE values must be positive.  The ROPE is assumed to be +/- the user supplied values.")
    if(length(ROPE) == ncol(X) - 1) ROPE = c(NA,ROPE)
  }
  ROPE_bounds = 
    paste("(",-round(ROPE,3),",",round(ROPE,3),")",sep="")
  

  # Get prior on $\beta$
  prior =
    c("zellner","normal","improper")[pmatch(tolower(prior),
                                               c("zellner","normal","improper"),duplicates.ok = FALSE)]


  if(prior == "zellner"){
    if(missing(zellner_g)){
      message("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
      zellner_g = N
    }

    prior_beta_mean = rep(0,p)
    prior_beta_precision = 1/zellner_g * crossprod(X)

  }

  if(prior == "normal"){

    if(missing(prior_beta_mean)){
      message("\nThe mu hyperparameter in the normal prior is not specified.  It will be set automatically to 0.\n")
      prior_beta_mean = rep(0,p)
    }

    if(missing(prior_beta_precision)){
      message("\nThe V hyperparameter in the normal prior is not specified.  It will be set automatically to 4/25Diag(s^2_{X_j})")

      prior_beta_precision =
        diag(1.0 / c(2.5 * max(s_y,1),2.5 * s_y / s_j)^2)

    }

  }
  
  if(prior == "improper"){
    prior_beta_mean = prior_beta_precision = NA
  }


  # Check for errors in family type and outcome
  ## Binomial
  if(family$family == "binomial"){
    if(length(unique(y)) == 2){
      if( (class(y) %in% c("numeric","integer") ) &
          ( !isTRUE(all.equal(0:1,sort(unique(y)))) ) ){
        paste0("Treating ",
               min(y),
               " as '0' and ",
               max(y),
               " as '1'") %>%
          message()
        y = ifelse(y == min(y),0,1)
      }
      if(class(y) == "character") y = factor(y)
      if(class(y) == "factor"){
        paste0("Treating ",
               levels(y)[2],
               " as '1' and ",
               levels(y)[1],
               " as '0'") %>%
          message()
        y = ifelse(y == levels(y)[1],0,1)
      }
    }
    if(missing(trials)){
      message("Assuming all observations correspond to Bernoulli, i.e., Binomial with one trial.")
      trials = rep(1.0,N)
    }else{
      if(class(trials) == "character") trials = data[[trials]]
      trials = as.numeric(trials)
    }
  }

  ## Poisson
  if(family$family == "poisson"){
    if( !(class(y) %in% c("numeric","integer")) ) stop("Outcome must be numeric or integer.")
    if(min(y) < 0) stop("Minumum of outcome must not be negative")
    trials = rep(1.0,N)
  }
  

  
  # Get initial start from glm
  if(family$family == "binomial"){
    init = glm(cbind(y,trials) ~ X[,-1] + offset(os),
               family)
  }else{
    init = glm(y ~ X[,-1] + offset(os),
               family)
  }

  # Get actual posterior mode from this
  ## Create posterior function
  log_posterior = function(x){
    eta = drop(X %*% x) + os
    mu = family$linkinv(eta)

    if(family$family == "binomial"){
      lpost =
        sum(
          dbinom(y,trials,mu,log = TRUE)
        )
    }
    if(family$family == "poisson"){
      lpost =
        sum(
          dpois(y,mu,log = TRUE)
        )
    }
    
    if(prior != "improper"){
      lpost =
        lpost -
        0.5 *
        drop(crossprod(
          x - prior_beta_mean,
          prior_beta_precision %*% (x - prior_beta_mean)
        ))
    }

    return(lpost)
  }
  
  ## Create gradient
  nabla_log_posterior = function(x){ 
    eta = drop(X %*% x) + os
    mu = family$linkinv(eta)

    if(family$family == "binomial"){
      grad_lpost =
        (y - trials * mu) %*% X
    }
    if(family$family == "poisson"){
      grad_lpost =
        (y - mu) %*% X
    }
    
    if(prior != "improper"){
      grad_lpost =
        drop( grad_lpost ) -
        drop( prior_beta_precision %*% (x - prior_beta_mean) )
    }

    return(drop(grad_lpost))
  }
  
  # Get posterior mode and hessian
  ## Run optim
  if( ( (family$family == "binomial") & (family$link == "logit") ) | 
      ( (family$family == "poisson") & (family$link == "log") ) ){
    opt = 
      optim(coef(init),
            log_posterior,
            nabla_log_posterior,
            method = "BFGS",
            hessian = TRUE,
            control = list(fnscale = -1))
  }else{
    opt = 
      optim(coef(init),
            log_posterior,
            method = "Nelder-Mead",
            hessian = TRUE,
            control = list(fnscale = -1))
  }
  ## Compute covariance matrix
  covmat = NULL
  try({
    covmat = 
      chol2inv(chol(-opt$hessian))
  }, silent=TRUE)
  if(is.null(covmat)){
    try({
      covmat = 
        qr.solve(-opt$hessian)
    }, silent=TRUE)
  }
  if(is.null(covmat)){
    try({
      covmat = 
        solve(-opt$hessian)
    }, silent=TRUE)
  }
  if(is.null(covmat)) stop("Hessian is not invertible.")
  
  
  if(is.na(n_draws)){
    # Collate results from large sample approx to return
    return_object = list()
    ## Summary
    return_object$summary = 
      tibble(Variable = colnames(X),
             `Post Mean` = opt$par,
             Lower = 
               qnorm(alpha / 2,
                     opt$par,
                     sd = sqrt(diag(covmat))),
             Upper = 
               qnorm(1 - alpha / 2,
                     opt$par,
                     sd = sqrt(diag(covmat))),
             `Prob Dir` = 
               pnorm(0,
                     opt$par,
                     sd = sqrt(diag(covmat))),
             ROPE = 
               pnorm(ROPE,
                     opt$par,
                     sqrt(diag(covmat))) - 
               pnorm(-ROPE,
                     opt$par,
                     sqrt(diag(covmat))),
             `ROPE bounds` = ROPE_bounds)
    
    return_object$summary$`Prob Dir` = 
      ifelse(return_object$summary$`Prob Dir` > 0.5,
             return_object$summary$`Prob Dir`,
             1.0 - return_object$summary$`Prob Dir`)
    
    ## Posterior covariance matrix
    return_object$posterior_covariance = covmat
    
    
  }else{#End: large sample summary
    
    # Perform importance sampling
    ## Get draws from the proposal
    proposal_draws = 
      mvtnorm::rmvt(n_draws,
                    delta = opt$par,
                    sigma = covmat,
                    df = proposal_df)
    # Get IS weights
    is_weights = 
      future_sapply(1:n_draws,
                    function(i){
                      log_posterior(proposal_draws[i,])
                    }
      ) -
      sum(
        mvtnorm::dmvt(proposal_draws,
                      delta = opt$par,
                      sigma = covmat,
                      df = proposal_df,
                      log = TRUE)
      )
    
    ## Stabilize and Normalize
    is_weights = exp(is_weights - max(is_weights))
    is_weights = is_weights / sum(is_weights)
    
    ## Get effective sample size
    ESS = 1.0 / sum( is_weights^2 )
    
    # Create helper function to get CI bounds
    CI_from_weighted_sample = function(x,w){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      LB = max(which(w <= 0.5 * alpha))
      UB = min(which(w >= 1.0 - 0.5 * alpha))
      return(c(lower = x[LB],
               upper = x[UB]))
    }
    ## Go ahead and get bounds
    CI_bounds = 
      apply(proposal_draws,2,
            CI_from_weighted_sample,
            w = is_weights)
    
    # Create helper function to get ROPE
    ROPE_from_weighted_sample = function(x,w,R){
      sum(w[(x > -R) & (x < R)])
    }
    
    # Create helper function to get PDir
    PDir_from_weighted_sample = function(x,w){
      max(sum(w[x > 0]),
          sum(w[x < 0]))
    }
    
    results =
      tibble(Variable = colnames(X), 
             `Post Mean` = apply(proposal_draws,2,
                                 weighted.mean,
                                 w = is_weights),
             Lower = CI_bounds["lower",],
             Upper = CI_bounds["upper",],
             `Prob Dir` = apply(proposal_draws,2,
                                PDir_from_weighted_sample,
                                w = is_weights),
             ROPE =
               sapply(1:ncol(X),
                      function(j){
                        ROPE_from_weighted_sample(proposal_draws[,j],
                                                  is_weights,
                                                  ROPE[j])
                      }),
             `ROPE bounds` = ROPE_bounds
      )
    
    return_object =
      list(summary = results,
           proposal_draws = proposal_draws,
           importance_sampling_weights = is_weights,
           effective_sample_size = ESS)
    
    
  }#End: importance sampling summary
  
  
  # Return hyperparameters and other user-inputs
  if(prior != "improper"){
    return_object$hyperparameters = 
      list(prior_beta_mean = prior_beta_mean,
           prior_beta_precision = prior_beta_precision)
  }else{
    return_object$hyperparameters = NA
  }
  return_object$formula = formula
  return_object$data = data
  return_object$family = family
  return_object$prior = prior
  return_object$ROPE = ROPE
  
  
  # Get fitted objects
  eta = drop(X %*% return_object$summary$`Post Mean`) + os
  mu = family$linkinv(eta)
  return_object$fitted =
    trials * mu
  
  # Get Pearson residuals
  SDs = 
    sqrt(
      trials * family$variance(mu)
    )
  return_object$residuals = 
    (y - trials * mu) / SDs
  
  rownames(return_object$summary) = NULL

  class(return_object) = "glm_b"
  return(return_object)

}