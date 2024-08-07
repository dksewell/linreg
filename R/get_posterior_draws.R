#' Get posterior samples from lm_b object
#' 
#' @param object Object of class lm_b
#' @param n_draws integer.  Number of posterior draws to obtain.
#' 
#' @export

get_posterior_draws = function(object, n_draws = 1e4){
  y = model.frame(object$formula,
                  object$data)[,1]
  X = model.matrix(object$formula,
                   object$data)
  
  p = nrow(object$summary)
  
  V_tilde_eig = eigen(object$post_parms$V_tilde)
  Vinv_sqrt = tcrossprod(diag(1 / sqrt(V_tilde_eig$values)),
                         V_tilde_eig$vectors)
  
  post_draws = 
    matrix(0.0,n_draws,p + 1,
           dimnames = list(NULL,
                           c(object$summary$Variable,"s2")))
  post_draws[,"s2"] = 
    rinvgamma(n_draws,
              0.5 * object$post_parms$a_tilde,
              0.5 * object$post_parms$b_tilde)
  post_draws[,1:p] = 
    matrix(1.0,n_draws,1) %*% matrix(object$summary$Post.Mean,nr=1) +
    matrix(rnorm(n_draws*p,
                 sd = sqrt(rep(post_draws[,"s2"],p))),n_draws,p) %*% Vinv_sqrt
  
  
  return(post_draws)
}