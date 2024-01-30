#' Find parameters for Beta prior based on prior mean and one quantile
#' 
#' @param mean numeric between 0 and 1 giving the prior mean
#' @param quantile numeric between 0 and 1 giving the quantile lying at left_tail_prob
#' @param left_tail_prob numeric between 0 and 1 giving the prior probability of theta 
#' being less than or equal to quantile
#' 
#' @examples
#' find_beta_parms(2/5,0.68,0.9)
#'
#' @export 


find_beta_parms = function(mean,quantile,left_tail_prob, plot_results = TRUE){
  if(missing(mean)) stop("Must provide mean, quantile, and the CDF at that quantile")
  helper = function(a){
    b = a * (1.0 - mean) / mean
    (pbeta(quantile,a,b) - left_tail_prob)^2
  }
  opt = optimize(f = helper,
                 interval = c(1e-3,100))
  a = opt$min
  b = a * (1.0 - mean) / mean
  
  if(plot_results){
    curve(dbeta(x,a,b),
          from= 0, 
          to = 1,
          lwd = 3,
          yaxt = "n",
          xlab = expression(theta))
  }
  
  
  return(c(shape1 = a, shape2 = b))
}
