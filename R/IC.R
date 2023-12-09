#' @name IC
#' @aliases AIC
#' @aliases BIC
#' @aliases DIC
#' 
#' @title Compute AIC, BIC, or DIC for aov_b or lm_b objects.  (Lower is better.)  
#' 
#' @export


DIC = function(object){
  UseMethod("DIC")
}

#' @rdname IC
#' @export
BIC.lm_b = function(object){
  y = model.frame(object$formula,
                  object$data)[,1]
  
  llik = 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$post_parms$b_tilde / 
                      (0.5 * object$post_parms$a_tilde + 1.0) ),
          log = TRUE) %>% 
    sum()
  
  -2.0 * llik + log(nrow(object$data)) * (length(object$post_parms$mu_tilde) + 1.0)
}

#' @rdname IC
#' @export
AIC.lm_b = function(object){
  y = model.frame(object$formula,
                  object$data)[,1]
  
  llik = 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$post_parms$b_tilde / 
                      (0.5 * object$post_parms$a_tilde + 1.0) ),
          log = TRUE) %>% 
    sum()
  
  -2.0 * llik + 2 * (length(object$post_parms$mu_tilde) + 1.0)
}

#' @rdname IC
#' @export
DIC.lm_b = function(object){
  y = model.frame(object$formula,
                  object$data)[,1]
  X = model.matrix(object$formula,
                   object$data)
  
  n_draws = 1e4
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
  
  
  D_E = 
    -2.0 * 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$post_parms$b_tilde / 
                      (0.5 * object$post_parms$a_tilde + 1.0)),
          log = TRUE) %>% 
    sum()
  
  
  llik = 
    sapply(1:nrow(object$data),
           function(i){
             dnorm(y[i],
                   mean = drop(tcrossprod(post_draws[,1:p],X[i,,drop=FALSE])),
                   sd = sqrt(post_draws[,p + 1]),
                   log = TRUE)
           })
  E_D = -2 * mean(rowSums(llik))
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}



#' @rdname IC
#' @export
BIC.aov_b = function(object){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    llik = 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) %>% 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    llik = 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) %>% 
      sum()
  }
  
  -2.0 * llik + log(nrow(object$data)) * nparms
}


#' @rdname IC
#' @export
AIC.aov_b = function(object){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    llik = 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) %>% 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    llik = 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) %>% 
      sum()
  }
  
  -2.0 * llik + 2 * nparms
}

#' @rdname IC
#' @export
DIC.aov_b = function(object){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    D_E = 
      -2.0 * 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) %>% 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    D_E = 
      -2.0 * 
      dnorm(object$data$y,
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) %>% 
      sum()
  }
  
  
  if(nparms == G+1){
    llik = 
      sapply(1:nrow(object$data),
             function(i){
               dnorm(object$data$y[i],
                     mean = object$posterior_draws[,as.integer(object$data$group)[i]],
                     sd = sqrt(object$posterior_draws[,G + 1]),
                     log = TRUE)
               })
  }else{
    llik = 
      sapply(1:nrow(object$data),
             function(i){
               dnorm(object$data$y[i],
                     mean = object$posterior_draws[,as.integer(object$data$group)[i]],
                     sd = sqrt(object$posterior_draws[,G + as.integer(object$data$group)[i]]),
                     log = TRUE)
             })
  }
  E_D = -2 * mean(rowSums(llik))
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}
