#' @name print
#' 
#' @title Print linreg objects.
#' 
#' @param x an object used to select a method.
#' 

#' @rdname print
#' @export
print.aov_b = function(x){
  cat("\n----------\n\nAnalysis of Variance fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @export
print.lm_b = function(x){
  cat("\n----------\n\nLinear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @export
print.np_lm_b = function(x){
  cat("\n----------\n\nNon-parametric linear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @export
print.bma = function(x){
  cat("\n----------\n\nBayesian model averaging for linear regression models\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}
