#' Find parameters for Inverse gamma prior based on prior mean and one quantile
#' 
#' @param median numeric greater than 0 giving the prior median
#' @param quantile numeric greater than 0 giving the quantile lying at left_tail_prob
#' @param left_tail_prob numeric between 0 and 1 giving the prior probability of theta 
#' being less than or equal to quantile
#' @param search_bounds bounds with which to search.  Sometimes you need to adjust this to get a good solution.
#' 
#' @examples
#' find_invgamma_parms(10 / 4,4.1,0.9)
#'
#' @import extraDistr
#' @export 


find_invgamma_parms = function(median,quantile,left_tail_prob, plot_results = TRUE){
  if(missing(median)) stop("Must provide median, quantile, and the CDF at that quantile")
  helper = function(x){
    a = exp(x[1])
    b = exp(x[2])
    (pinvgamma(median,a,b) - 0.5)^2 +
      (pinvgamma(quantile,a,b) - left_tail_prob)^2
  }
  opt = optim(numeric(2),
              f = helper)
  a = exp(opt$par[1])
  b = exp(opt$par[2])
  
  if(plot_results){
    curve(dinvgamma(x,a,b),
          from= 0, 
          to = qinvgamma(0.99,a,b),
          lwd = 3,
          yaxt = "n",
          xlab = expression(theta))
  }
  
  
  return(c(shape1 = a, shape2 = b) * 2)
}

