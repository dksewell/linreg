#' Non-parametric linear models
#' 
#' np_lm_b uses general Bayesian inference with loss-likelihood bootstrap. 
#' This is, as implemented here, a Bayesian non-parametric linear models 
#' inferential engine.  
#' 
#' 
#' Applicable data types are continuous (use family = 
#' gaussian()), count (use family = poisson()), or binomial 
#' (use family = binomial()). 
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.
#' @param family A description of the error distribution and link function 
#' to be used in the model. See ?glm for more information.
#' @param loss Either "selfinformation", "gls" (for generalized squared error), 
#' or a function that takes in two arguments, the first of which should 
#' be the vector of outcomes and the second should be the expected value of y; 
#' The outcome of the function should be the loss evaluated for each observation. 
#' By default, the self-information loss is used (i.e., the negative log-likelihood).  
#' Note: I really do mean the expected value of y, even for binomial (i.e., n*p).
#' @param loss_gradient If loss is a user-defined function (as opposed to 
#' "selfinformation" or "gls"), supplying the gradient to the loss will 
#' speed up the algorithm. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param n_draws integer.  Number of posterior draws to obtain.  If left missing, 
#' the large sample approximation will be used.
#' @param CI_level numeric. Credible interval level.
#' @param seed integer.  Always set your seed!!!
#' 
#' @return np_lm_b() returns an object of class "np_lm_b", which behaves as
#' a list with the following elements:
#' \itemize{
#'  \item summary - a tibble giving results for regression coefficients.
#' }
#' 
#' @import future
#' @import future.apply
#' @import Matrix
#' @export 
#' @exportClass np_lm_b

np_lm_b = function(formula,
                   data,
                   family,
                   loss = "selfinformation",
                   loss_gradient,
                   trials,
                   n_draws,
                   CI_level = 0.95,
                   seed = 1){
  
  # Get correct family
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  # Get loss function
  if(class(loss) == "function"){
    loss_fun = loss
    loss = "custom"
  }
  
  if(loss == "selfinformation"){
    if(family$family == "gaussian"){
      loss_fun = function(y,mu) (y - mu)^2
    }
    if(family$family == "binomial"){
      loss_fun = function(y,mu) -dbinom(y,trials,mu/trials,log=T)
    }
    if(family$family == "poisson"){
      loss_fun = function(y,mu) -dpois(y,mu,log=T)
    }
  }
  
  if(loss == "gls"){
    loss_fun = function(y,mu) (y - mu)^2  / (trials * family$variance(mu))
  }
  
  # Extract 
  mframe = model.frame(formula, data)
  y = model.response(mframe)
  X = model.matrix(formula,data)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  
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
  
  ## Gaussian
  if(family$family == "gaussian"){
    if( !(class(y) %in% c("numeric","integer")) ) stop("Outcome must be numeric.")
    trials = rep(1.0,N)
  }
  
  
  
  # Create total loss function
  loss_wrapper = function(x,w){
    eta = drop(X %*% x) + os
    mu = family$linkinv(eta)
    
    weighted.mean(
      loss_fun(y, trials * mu),
      w
    )
  }
  
  # Get total loss gradient
  if(missing(loss_gradient)) loss_gradient = NULL
  ## self-information loss
  if(loss == "selfinformation"){
    if( (family$family == "gaussian") & (family$link == "identity")){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        mu = family$linkinv(eta)
        
        apply(-2.0 * (y - mu) * X,
              2,
              weighted.mean,
              w = w)
      }
    }
    if( (family$family == "poisson") & (family$link == "log") ){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        mu = family$linkinv(eta)
        
        apply( (mu - y) * X,
               2,
               weighted.mean,
               w = w)
      }
    }
    if( (family$family == "binomial") & (family$link == "logit") ){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        probs = family$linkinv(eta)
        
        apply( ( trials * probs - y ) * X,
               2,
               weighted.mean,
               w = w)
      }
    }
  }
  
  ## GLS loss
  if(loss == "gls"){
    if( (family$family == "gaussian") & (family$link == "identity")){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        mu = family$linkinv(eta)
        
        apply(-2.0 * (y - mu) * X,
              2,
              weighted.mean,
              w = w)
      }
    }
    if( (family$family == "poisson") & (family$link == "log") ){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        mu = family$linkinv(eta)
        
        apply( ((mu^2 - y^2) / mu) * X,
               2,
               weighted.mean,
               w = w)
      }
    }
    if( (family$family == "binomial") & (family$link == "logit") ){
      loss_gradient = function(x,w){
        eta = drop(X %*% x) + os
        mu = trials * family$linkinv(eta)
        
        apply( ( (mu - y) * ( trials * (mu + y) - 2.0 * mu * y ) /
                   ( mu * (trials - mu) ) ) * X,
               2,
               weighted.mean,
               w = w)
      }
    }
  }
  
  
  
  # Get empirical minimizer
  if(!is.null(loss_gradient)){
    empir_min = 
      optim(glm(I(y/trials) ~ X[,-1],
                family = family,
                weights = trials) %>% coef(),
            loss_wrapper,
            loss_gradient,
            method = "CG",
            w = rep(1.0,N),
            control = list(maxit = 1e4))$par
  }else{
    empir_min = 
      optim(glm(I(y/trials) ~ X[,-1],
                family = family,
                weights = trials) %>% coef(),
            loss_wrapper,
            method = "Nelder-Mead",
            w = rep(1.0,N),
            control = list(maxit = 1e4))$par
  }
  
  
  # Obtain general posterior parameters
  if(missing(n_draws)){
    
    ## Find I and J matrices
    ### self-information
    if(loss == "selfinformation"){
      if(family$family == "gaussian"){
        eta = drop(X %*% empir_min) + os
        mu = family$linkinv(eta)
        
        I_beta0 = 
          crossprod(-2.0 * (y - mu) * X) / N
        J_beta0 = 
          crossprod(X) / N
      }
      
      if(family$family == "poisson"){
        eta = drop(X %*% empir_min) + os
        mu = family$linkinv(eta)
        W1 = 
          Diagonal(x = mu - y)
        W2 = 
          Diagonal(x = mu)
        
        I_beta0 = 
          crossprod(W1 %*% X) / N
        J_beta0 = 
          crossprod(X,W2 %*% X) / N
      }
      
      if(family$family == "binomial"){
        eta = drop(X %*% empir_min) + os
        probs = family$linkinv(eta)
        mu = trials * probs
        W1 = 
          Diagonal(x = mu - y)
        W2 = 
          Diagonal(x = mu * (1.0 - probs))
        
        I_beta0 = 
          crossprod(W1 %*% X) / N
        J_beta0 = 
          crossprod(X,W2 %*% X) / N
      }
    }
    
    ### GLS
    if(loss == "gls"){
      if(family$family == "gaussian"){
        eta = drop(X %*% empir_min) + os
        mu = family$linkinv(eta)
        
        I_beta0 = 
          crossprod(-2.0 * (y - mu) * X) / N
        J_beta0 = 
          crossprod(X) / N
      }
      
      if(family$family == "poisson"){
        eta = drop(X %*% empir_min) + os
        mu = family$linkinv(eta)
        W = 
          Diagonal(x = 
                     (mu^2 + y^2) / mu
          )
        
        I_beta0 = 
          crossprod(((mu^2 - y^2) / mu) * X) / N
        J_beta0 = 
          crossprod(X,W %*% X) / N
      }
      
      if(family$family == "binomial"){
        eta = drop(X %*% empir_min) + os
        mu = trials * family$linkinv(eta)
        W = 
          Diagonal(x = 
                     ( 
                       trials^2 * (mu^2 + y^2) - 
                         2.0 * trials * mu * y * (mu + y) +
                         2.0 * mu^2 * y^2
                     ) /
                     (
                       trials * mu * (trials - mu)
                     )
          )
        
        I_beta0 = 
          crossprod(( (mu - y) * ( trials * (mu + y) - 2.0 * mu * y ) /
                        ( mu * (trials - mu) ) ) * X) / N
        J_beta0 = 
          crossprod(X,W %*% X) / N
      }
    }
    
    ## Compute covariance matrix
    J_inv = NULL
    try({
      J_inv = 
        chol2inv(chol(J_beta0))
    }, silent=TRUE)
    if(is.null(J_inv)){
      try({
        J_inv = 
          qr.solve(J_beta0)
      }, silent=TRUE)
    }
    if(is.null(J_inv)){
      try({
        J_inv = 
          solve(J_beta0)
      }, silent=TRUE)
    }
    if(is.null(J_inv)) stop("J(beta_0) is not invertible.  Use bootstrapping instead.")
    
    covmat = 
      J_inv %*% I_beta0 %*% J_inv
    
    
    # Collate results from large sample approx to return
    results = list()
    ## Summary
    results$summary = 
      tibble(Variable = colnames(X),
             `Post Mean` = empir_min,
             Lower = 
               qnorm(alpha / 2,
                     empir_min,
                     sd = sqrt(diag(covmat) / N)),
             Upper = 
               qnorm(1 - alpha / 2,
                     empir_min,
                     sd = sqrt(diag(covmat) / N)),
             `Prob Dir` = 
               pnorm(0,
                     empir_min,
                     sd = sqrt(diag(covmat) / N)))
    
    results$summary$`Prob Dir` = 
      ifelse(results$summary$`Prob Dir` > 0.5,
             results$summary$`Prob Dir`,
             1.0 - results$summary$`Prob Dir`)
    
    ## Posterior covariance matrix
    results$posterior_covariance = covmat / N
    
    
  }else{
    
    
    # Perform Bayesian loss-likelihood bootstrap
    if("sequential" %in% class(plan())){
      # Get weights for Bayesian bootstrap
      set.seed(seed)
      dirichlet_draws = 
        rdirichlet(n_draws,
                   rep(1.0,N))
      beta_draws = matrix(0.0,n_draws,p)
      for(i in 1:n_draws){
        if(!is.null(loss_gradient)){
          temp =
            optim(empir_min,
                  loss_wrapper,
                  loss_gradient,
                  method = "CG",
                  w = dirichlet_draws[i,],
                  control = list(maxit = 1e4))
        }else{
          temp =
            optim(empir_min,
                  loss_wrapper,
                  method = "Nelder-Mead",
                  w = dirichlet_draws[i,],
                  control = list(maxit = 1e4))
        }
        if(temp$conv == 0){
          beta_draws[i,] = temp$par
        }else{
          beta_draws[i,] = NA
        }
      }
      
    }else{
      
      helper = function(i){
        if(!is.null(loss_gradient)){
          temp =
            optim(empir_min,
                  loss_wrapper,
                  loss_gradient,
                  method = "CG",
                  w = drop( rdirichlet(1, rep(1.0,N)) ),
                  control = list(maxit = 1e4))
        }else{
          temp =
            optim(empir_min,
                  loss_wrapper,
                  method =  "Nelder-Mead",
                  w = drop( rdirichlet(1, rep(1.0,N)) ),
                  control = list(maxit = 1e4))
        }
        if(temp$conv == 0){
          return(temp$par)
        }else{
          return(rep(NA,p))
        }
      }
      
      beta_draws = 
        future_sapply(1:n_draws,
                      helper,
                      future.seed = seed) %>% 
        t() %>% 
        na.omit()
    }
    
    
    
    # Collate results from loss-lik bootstrapping to return
    results = list()
    ## Summary
    results$summary = 
      tibble(Variable = colnames(X),
             `Post Mean` = colMeans(beta_draws),
             Lower = 
               beta_draws %>% 
               na.omit() %>% 
               apply(2,quantile,prob = alpha / 2),
             Upper = 
               beta_draws %>% 
               apply(2,quantile,prob = 1.0 - alpha / 2),
             `Prob Dir` = 
               beta_draws %>%  
               na.omit() %>% 
               apply(2,function(x) mean(x > 0))
      )
    results$summary$`Prob Dir` = 
      ifelse(results$summary$`Prob Dir` > 0.5,
             results$summary$`Prob Dir`,
             1.0 - results$summary$`Prob Dir`)
    
    ## Posterior draws
    results$posterior_draws = na.omit(beta_draws)
    colnames(results$posterior_draws) = results$summary$Variable
    
  }
  
  
  ## Fitted values
  eta = X %*% results$summary$`Post Mean` + os
  results$fitted =
    trials * family$linkinv(eta)
  
  ## Pearson residuals
  results$residuals = 
    (y - results$fitted) /
    sqrt(trials * family$variance(family$linkinv(eta)))
  
  ## Input values
  results$formula = formula
  results$data = data
  results$family = family
  
  
  class(results) = "np_lm_b"
  
  return(results)
}