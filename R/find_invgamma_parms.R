#' Find parameters for Inverse gamma prior based on prior mean and one quantile
#' 
#' @param mean numeric greater than 0 giving the prior mean
#' @param quantile numeric greater than 0 giving the quantile lying at left_tail_prob
#' @param left_tail_prob numeric between 0 and 1 giving the prior probability of theta 
#' being less than or equal to quantile
#' 
#' @examples
#' find_invgamma_parms(10 / 4,4.1,0.9)
#'
#' @import extraDistr
#' @export 


find_invgamma_parms = function(mean,quantile,left_tail_prob, plot_results = TRUE){
  if(missing(mean)) stop("Must provide mean, quantile, and the CDF at that quantile")
  helper = function(a){
    b = mean * (a - 1)
    (pinvgamma(quantile,a,b) - left_tail_prob)^2
  }
  opt = optimize(f = helper,
                 interval = c(1 + 1e-3,5e2))
  a = opt$min
  b = mean * (a - 1)
  
  if(plot_results){
    curve(dinvgamma(x,a,b),
          from= 0, 
          to = qinvgamma(0.99,a,b),
          lwd = 3,
          yaxt = "n",
          xlab = expression(theta))
  }
  
  
  return(c(shape1 = a, shape2 = b))
}
