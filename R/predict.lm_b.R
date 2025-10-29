#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level numeric. Credible interval level.
#' 
#' @return tibble with estimate, prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @export predict.lm_b
#' @export


predict.lm_b = function(object,
                        newdata,
                        CI_level = 0.95,
                        n_draws = 0){
  
  alpha = 1.0 - CI_level
  
  if(missing(newdata)){
    newdata = object$data
  }
  X = model.matrix(as.formula(paste(as.character(object$formula)[c(1,3)],
                                    collapse = "")),
                   data = newdata)
  
  V_eig = eigen(object$post_parms$V_tilde)
  XULneghalf = X %*% V_eig$vectors %*% diag(1 / sqrt(V_eig$values))
  X_VInverse_X = 
    apply(XULneghalf,1,crossprod)
  
  newdata %<>%
    mutate(Estimate = drop(X %*% object$post_parms$mu_tilde)) %>% 
    mutate(PI_lower = 
             qlst(alpha/2.0,
                  df = object$post_parms$a_tilde,
                  mu = Estimate,
                  sigma = sqrt(object$post_parms$b_tilde / 
                                 object$post_parms$a_tilde * 
                                 (X_VInverse_X + 1.0) ) ),
           PI_upper = 
             qlst(1.0 - alpha/2.0,
                  df = object$post_parms$a_tilde,
                  mu = Estimate,
                  sigma = sqrt(object$post_parms$b_tilde / 
                                 object$post_parms$a_tilde * 
                                 (X_VInverse_X + 1.0) ) ),
           CI_lower = 
             qlst(alpha/2.0,
                  df = object$post_parms$a_tilde,
                  mu = Estimate,
                  sigma = sqrt(object$post_parms$b_tilde / 
                                 object$post_parms$a_tilde * 
                                 (X_VInverse_X) ) ),
           CI_upper = 
             qlst(1.0 - alpha/2.0,
                  df = object$post_parms$a_tilde,
                  mu = Estimate,
                  sigma = sqrt(object$post_parms$b_tilde / 
                                 object$post_parms$a_tilde * 
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
             df = object$post_parms$a_tilde,
             mu = newdata$Estimate,
             sigma = sqrt(object$post_parms$b_tilde / 
                            object$post_parms$a_tilde * 
                            (X_VInverse_X + 1.0) ) )
    }
    
    newdata %<>%
      bind_cols(new_draws %>% 
                  as_tibble())
  }
  
  return(newdata)
}
