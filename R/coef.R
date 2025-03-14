#' @name coef
#' 
#' @title Coefficient extraction for lm_b or aov_b object
#' 
#' @param object lm_b or aov_b object
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