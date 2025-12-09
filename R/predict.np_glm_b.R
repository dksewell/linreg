#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param CI_level numeric. Credible interval level.
#' @param ... optional arguments.
#' 
#' @return tibble with estimate, prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @exportS3Method predict np_glm_b


predict.np_glm_b = function(object,
                            newdata,
                            trials,
                            CI_level = 0.95,
                            ...){
  
  
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
    
    if(object$family$family == "gaussian"){
      grad_ginv_xbeta = 
        X
    }
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
    
  }else{
    yhat_draws = 
      trials * 
      object$family$linkinv(os + tcrossprod(X, object$posterior_draws))
    
    newdata =
      newdata |> 
      as_tibble() |> 
      mutate(`Post Mean` = yhats,
             CI_lower = 
               yhat_draws |> 
               apply(1,quantile, probs = alpha / 2.0),
             CI_upper = 
               yhat_draws |> 
               apply(1,quantile, probs = 1.0 - alpha / 2.0))
  }
  
  
  return(newdata)
}
