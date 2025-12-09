#' @name print
#' 
#' @title Print linreg objects.
#' 
#' @param x an object used to select a method.
#' @param ... optional arguments.

#' @rdname print
#' @method print aov_b 
#' @export
print.aov_b = function(x, ...){
  cat("\n----------\n\nAnalysis of Variance fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n") 
  if("BF_for_different_vs_same_means" %in% names(x)){
    cat(paste0(
      "Bayes factor in favor of the full vs. null model: ",
      format(signif(x$BF_for_different_vs_same_means, 3), 
             scientific = 
               (x$BF_for_different_vs_same_means > 1e3) | 
               (x$BF_for_different_vs_same_means < 1e-3))))
    
    cat("\n\n----------\n\n")
  }
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @method print lm_b 
#' @export
print.lm_b = function(x, ...){
  cat("\n----------\n\nLinear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @method print np_glm_b
#' @export
print.np_glm_b = function(x, ...){
  cat("\n----------\n\nNon-parametric linear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @method print lm_b_bma
#' @export
print.lm_b_bma = function(x, ...){
  cat("\n----------\n\nBayesian model averaging for linear regression models\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}

#' @rdname print
#' @method print glm_b
#' @export
print.glm_b = function(x, ...){
  cat("\n----------\n\nGeneralized linear regression fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
}
