#' Predict method for bma model fits
#' 
#' 
#' @param object Object of class bma
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level numeric. Credible interval level.
#' @param seed integer. Always set your seed!!!
#' 
#' @return tibble with estimate, prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @export predict.bma
#' @export


predict.bma = function(object,
                       newdata,
                       CI_level = 0.95,
                       seed = 1){
  
  if(missing(newdata)){
    newdata = object$data
  }
  
  # Extract 
  X = model.matrix(object$formula,newdata)
  N = nrow(X)
  p = ncol(X)
  alpha = 1 - CI_level
  
  # Get means of E(y|X)
  mu_draws = 
    tcrossprod(as.matrix(object$posterior_draws[,1:ncol(X)]),
               X)
  # Convert from variance to sd
  object$posterior_draws$s2 = sqrt(object$posterior_draws$s2)
  
  # Get draws of y|X
  y_new_draws = 
    future_sapply(1:ncol(mu_draws),
                  function(i){
                    rnorm(nrow(mu_draws),
                          mu_draws[,i],
                          object$posterior_draws$s2)
                  },
                  future.seed = seed)
    
  # Compile results
  newdata %<>% 
    mutate(`Post Mean` = 
             colMeans(mu_draws),
           CI_Lower = 
             apply(mu_draws,2,quantile,probs = 0.5 * alpha),
           CI_upper = 
             apply(mu_draws,2,quantile,probs = 1.0 - 0.5 * alpha),
           PI_Lower = 
             apply(y_new_draws,2,quantile,probs = 0.5 * alpha),
           PI_upper = 
             apply(y_new_draws,2,quantile,probs = 1.0 - 0.5 * alpha))
  
  return(list(newdata = newdata,
              posterior_draws = 
                list(mean_of_ynew = mu_draws,
                     ynew = y_new_draws)))
}
