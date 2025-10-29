#' @name coef
#' 
#' @title Coefficient extraction for linreg objects
#' 
#' @param object linreg object
#' 
#' 

#' @rdname coef
#' @export
coef.lm_b = function(object){
  object$summary$Post.Mean
}

#' @rdname coef
#' @export
coef.aov_b = function(object){
  object$posterior_parameters$mu_g
}

#' @rdname coef
#' @export
coef.np_lm_b = function(object){
  ret = object$summary$`Post Mean`
  names(ret) = object$summary$Variable
  ret
}
