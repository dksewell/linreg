#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param n_draws If desired, the number of posterior samples drawn.
#' 
#' @return tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @export predict.lm_b
#' @export


predict.lm_b = function(object,
                        newdata,
                        CI_level = 0.95,
                        PI_level = 0.95,
                        n_draws = 0){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  if(missing(newdata)){
    newdata = object$data
  }
  X = model.matrix(as.formula(paste(as.character(object$formula)[c(1,3)],
                                    collapse = "")),
                   data = newdata)
  
  V_eig = eigen(object$posterior_parameters$V_tilde)
  XULneghalf = X %*% V_eig$vectors %*% diag(1 / sqrt(V_eig$values))
  X_VInverse_X = 
    apply(XULneghalf,1,crossprod)
  
  newdata %<>%
    mutate(`Post Mean` = drop(X %*% object$posterior_parameters$mu_tilde)) %>% 
    mutate(PI_lower = 
             qlst(alpha_pi/2.0,
                  df = object$posterior_parameters$a_tilde,
                  mu = `Post Mean`,
                  sigma = sqrt(object$posterior_parameters$b_tilde / 
                                 object$posterior_parameters$a_tilde * 
                                 (X_VInverse_X + 1.0) ) ),
           PI_upper = 
             qlst(1.0 - alpha_pi/2.0,
                  df = object$posterior_parameters$a_tilde,
                  mu = `Post Mean`,
                  sigma = sqrt(object$posterior_parameters$b_tilde / 
                                 object$posterior_parameters$a_tilde * 
                                 (X_VInverse_X + 1.0) ) ),
           CI_lower = 
             qlst(alpha_ci/2.0,
                  df = object$posterior_parameters$a_tilde,
                  mu = `Post Mean`,
                  sigma = sqrt(object$posterior_parameters$b_tilde / 
                                 object$posterior_parameters$a_tilde * 
                                 (X_VInverse_X) ) ),
           CI_upper = 
             qlst(1.0 - alpha_ci/2.0,
                  df = object$posterior_parameters$a_tilde,
                  mu = `Post Mean`,
                  sigma = sqrt(object$posterior_parameters$b_tilde / 
                                 object$posterior_parameters$a_tilde * 
                                 (X_VInverse_X) ) )
    )
  
  # Get new draws
  if(n_draws > 0){
    new_draws = matrix(0.0,
                       nrow(newdata),
                       n_draws,
                       dimnames = list(NULL,
                                       paste("y_new",1:n_draws,sep="")))
    for(it in 1:n_draws){
      new_draws[,it] = 
        rlst(nrow(newdata),
             df = object$posterior_parameters$a_tilde,
             mu = newdata$`Post Mean`,
             sigma = sqrt(object$posterior_parameters$b_tilde / 
                            object$posterior_parameters$a_tilde * 
                            (X_VInverse_X + 1.0) ) )
    }
    
    newdata %<>%
      bind_cols(new_draws %>% 
                  as_tibble())
  }
  
  return(newdata)
}
