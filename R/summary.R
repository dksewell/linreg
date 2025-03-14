#' @name summary
#' 
#' @title Summary functions for lm_b and aov_b objects
#' 
#' @param object lm_b or aov_b object
#' 
#' @export

#' @rdname summary
#' @export
summary.lm_b = function(object){
  object$summary
}

#' @rdname summary
#' @export
summary.aov_b = function(object){
  cat("\n--- Summary of factor level means ---\n")
  print(object$summary)
  
  cat("\n\n\n\n--- Summary of pairwise differences ---\n")
  
  print(object$pairwise_summary)
  
  invisible(list(summary = object$summary,
                 pairwise = object$pairwise_summary))
}