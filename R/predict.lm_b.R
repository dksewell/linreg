#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level numeric. Credible interval level.
#' 
#' @return matrix with estimate and prediction intervals.
#' 
#' @export predict.lm_b
#' @export


predict.lm_b = function(object,
                        newdata,
                        CI_level = 0.95){
  
  alpha = 1.0 - CI_level
  
  if(missing(newdata)){
    newdata = object$data
  }
  X = model.matrix(object$formula,
                   data = newdata)
  
  V_eig = eigen(object$post_parms$V_tilde)
  XULneghalf = X %*% V_eig$vectors %*% diag(1 / sqrt(V_eig$values))
  X_VInverse_X = 
    apply(XULneghalf,1,crossprod)
  
  ret = 
    cbind(Estimate = drop(object$fitted),
          Lower = 
            qlst(alpha/2.0,
                 df = object$post_parms$a_tilde,
                 mu = object$fitted,
                 sigma = sqrt(object$post_parms$b_tilde / 
                                object$post_parms$a_tilde * 
                                (X_VInverse_X + 1.0) ) ),
          Upper = 
            qlst(1.0 - alpha/2.0,
                 df = object$post_parms$a_tilde,
                 mu = object$fitted,
                 sigma = sqrt(object$post_parms$b_tilde / 
                                object$post_parms$a_tilde * 
                                (X_VInverse_X + 1.0) ) )
          )
  
  return(ret)
}