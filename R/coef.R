#' @name coef
#' 
#' @title Coefficient extraction for linreg objects
#' 
#' @param object linreg object
#' @param ... optional arguments.
#' 

#' @rdname coef
#' @export
coef.lm_b = function(object, ...){
  object$summary$`Post Mean`
}

#' @rdname coef
#' @export
coef.aov_b = function(object, ...){
  object$posterior_parameters$mu_g
}

#' @rdname coef
#' @export
coef.np_glm_b = function(object, ...){
  ret = object$summary$`Post Mean`
  names(ret) = object$summary$Variable
  ret
}

#' @rdname coef
#' @export
coef.glm_b = function(object, ...){
  ret = object$summary$`Post Mean`[1:(nrow(object$summary) -
                                        (object$family$family == "negbinom") )]
  names(ret) = object$summary$Variable[1:(nrow(object$summary) -
                                            (object$family$family == "negbinom") )]
  ret
}

#' @rdname coef
#' @export
coef.lm_b_bma = function(object, ...){
  ret = object$summary$`Post Mean`
  names(ret) = object$summary$Variable
  ret[-length(ret)]
}
