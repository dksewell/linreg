#' Coefficient extraction for lm_b or aov_b object
#' 
#' @param object lm_b or aov_b object
#' 
#' 
#' @export coef.lm_b
#' @export coef.aov_b
#' @exportS3Method linreg::coef


coef.lm_b = function(object){
  object$summary$Post.Mean
}

coef.aov_b = function(object){
  object$posterior_parameters$mu_g
}