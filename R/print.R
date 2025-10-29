#' @name print
#' 
#' @title Print linreg objects.
#' 
#' @export

#' @rdname print
#' @export
print.aov_b = function(object){
  cat("\n----------\n\nAnalysis of Variance fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(object$formula)
  cat("\n----------\n\n")
  print(object$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @export
print.lm_b = function(object){
  cat("\n----------\n\nLinear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(object$formula)
  cat("\n----------\n\n")
  print(object$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @export
print.np_lm_b = function(object){
  cat("\n----------\n\nNon-parametric linear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(object$formula)
  cat("\n----------\n\n")
  print(object$summary)
  cat("\n----------\n")
}
