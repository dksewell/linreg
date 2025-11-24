#' t-test
#' 
#' One and two sample t-tests on vectors of data
#' 
#' @param x a (non-empty) numeric vector of data values
#' @param y an optional (non-empty) numeric vector of data values
#' @param mu optional.  If supplied, \code{t_test_b} will return the 
#' posterior probabilty that the population mean (ignored in 2 sample inference) 
#' is less than this value.
#' @param paired logical.  If TRUE, provide both x and y as vectors.
#' @param plot description
#' @param ... arguments to be passed to \code{\ling{aov_b}}.
#' 
#' @export


t_test_b = function(x,
                    y,
                    mu,
                    paired = FALSE,
                    data,
                    ...){
  
  if(is.formula(x) & missing(data)) stop("If formula is given, data must also be given.")
  if(missing(y) & paired) stop("Cannot have paired data without y.")
  
  if(is.numeric(x)){
    if(missing(y) | paired){
      if(paired) x = x - y
      # Copy lots over from aov_b for one sample inference.  Do it on x.
      
      
    }else{
      ttest_data = 
        tibble(group = rep(c("x","y"),
                           c(length(x),
                             length(y)))) %>% 
        mutate(y = c(x,y))
      
      return(aov_b(y ~ group,
                   data = ttest_data,
                   ...))
    }
  }else{
    return(aov_b(x,
                 data = data,
                 ...))
  }
  
  
  
}