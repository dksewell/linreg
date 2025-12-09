#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param seed integer.  Always set your seed!!!
#' @param ... optional arguments.
#' 
#' @return tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @import mvtnorm
#' @exportS3Method predict glm_b


predict.glm_b = function(object,
                         newdata,
                         trials,
                         CI_level = 0.95,
                         PI_level = 0.95,
                         seed = 1,
                         ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  
  if(missing(newdata)){
    newdata = object$data
  }
  
  # Extract 
  mframe = model.frame(object$formula, newdata)
  X = model.matrix(object$formula,newdata)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  # Get trials variables sorted
  if(object$family$family == "binomial"){
    if(missing(trials)){
      message("Assuming all observations correspond to Bernoulli, i.e., Binomial with one trial.")
      trials = rep(1.0,N)
    }else{
      if(class(trials) == "character") trials = newdata[[trials]]
      trials = as.numeric(trials)
    }
  }else{
    trials = rep(1.0,N)
  }
  
  
  # Assign new values
  ## Get estimates
  yhats = 
    trials * 
    object$family$linkinv(eta = 
                            drop(X %*% 
                                   object$summary$`Post Mean`) + 
                            os)
  
  
  ## Get CI bounds
  if("posterior_covariance" %in% names(object)){
    
    if(object$family$family == "poisson"){
      grad_ginv_xbeta = 
        drop(exp(X %*% object$summary$`Post Mean` + os)) * X
    }
    if(object$family$family == "binomial"){
      probs = 
        1.0 / (1.0 + drop(exp(-X %*% object$summary$`Post Mean` - os)))
      grad_ginv_xbeta = 
        trials * probs * (1.0 - probs) * X
    }
    
    yhats_covar = 
      tcrossprod(grad_ginv_xbeta %*% object$posterior_covariance,
                 grad_ginv_xbeta)
    
    newdata =
      newdata |> 
      mutate(`Post Mean` = yhats,
             CI_lower = 
               qnorm(alpha / 2.0,
                     yhats,
                     sqrt(diag(yhats_covar))),
             CI_upper = 
               qnorm(1.0 - alpha / 2.0,
                     yhats,
                     sqrt(diag(yhats_covar))))
    
    yhats_sds = 
      sqrt(diag(yhats_covar))
    
    
    if(object$family$family == "poisson"){
      y_draws =
        future_sapply(1:5e3, # Might set this as an argument later.
                      function(i){
                        rpois(nrow(newdata),
                              pmax(rnorm(nrow(newdata),
                                         newdata$`Post Mean`,
                                         yhats_sds),
                                   .Machine$double.eps)
                        )
                      },
                      future.seed = seed)
      
    }
    if(object$family$family == "binomial"){
      y_draws =
        future_sapply(1:5e3, # Might set this as an argument later.
                      function(i){
                        rbinom(nrow(newdata),
                               trials,
                               pmin(
                                 pmax(
                                   rnorm(nrow(newdata),
                                         newdata$`Post Mean`,
                                         yhats_sds),
                                   .Machine$double.eps),
                                 1.0 - .Machine$double.eps)
                               )
                      },
                      future.seed = seed)
      
    }
    
    PI_bounds = 
      y_draws |> 
      future_apply(1,quantile,probs = c(0.5 * alpha_pi,
                                        1.0 - 0.5 * alpha_pi))
    newdata$PI_lower = 
      PI_bounds[1,]
    newdata$PI_upper = 
      PI_bounds[2,]
    
  }else{#End: if asymptotic approx was used.
    
    yhat_draws = 
      trials * 
      object$family$linkinv(os + tcrossprod(X, object$proposal_draws))
    
    CI_from_weighted_sample = function(x,w){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      LB = max(which(w <= 0.5 * alpha_ci))
      UB = min(which(w >= 1.0 - 0.5 * alpha_ci))
      return(c(lower = x[LB],
               upper = x[UB]))
    }
    CI_bounds = 
      apply(yhat_draws,1,
            CI_from_weighted_sample,
            w = object$importance_sampling_weights)
    
    newdata =
      newdata |> 
      as_tibble() |> 
      mutate(`Post Mean` = yhats,
             CI_lower = 
               CI_bounds["lower",],
             CI_upper = 
               CI_bounds["upper",])
    
    if(object$family$family == "poisson"){
      y_draws = 
        future_sapply(1:nrow(yhat_draws),
                      function(i){
                        rpois(ncol(yhat_draws),yhat_draws[i,])
                      },
                      future.seed = seed)
      PI_bounds = 
        y_draws |> 
        future_apply(2,quantile,probs = c(0.5 * alpha_pi,
                                          1.0 - 0.5 * alpha_pi))
      newdata$PI_lower = 
        PI_bounds[1,]
      newdata$PI_upper = 
        PI_bounds[2,]
    }
    
    if(object$family$family == "binomial"){
      y_draws = 
        future_sapply(1:nrow(yhat_draws),
                      function(i){
                        rbinom(ncol(yhat_draws),
                               trials,
                               yhat_draws[i,])
                      },
                      future.seed = seed)
      PI_bounds = 
        y_draws |> 
        future_apply(2,quantile,probs = c(0.5 * alpha_pi,
                                          1.0 - 0.5 * alpha_pi))
      newdata$PI_lower = 
        PI_bounds[1,]
      newdata$PI_upper = 
        PI_bounds[2,]
    }
    
  }
  
  
  return(newdata)
}
