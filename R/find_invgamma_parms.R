#' Find parameters for Inverse gamma prior based on prior mean and one quantile
#' 
#' @param lower_quantile lower quantile desired
#' @param upper_quantile upper quantile desired
#' @param probability prior probability to be contained within the lower and upper quantiles
#' @param search_bounds bounds with which to search.  Sometimes you need to adjust this to get a good solution.
#' 
#' @returns twice the shape and rate of the inverse gamma distribution.
#' 
#' @examples
#' find_invgamma_parms(10 / 4,4.1,0.9)
#' extraDistr::qinvgamma(c(0.5,0.9),
#'                         16.25/2,
#'                         39/2)
#'
#' @import extraDistr
#' @export 


find_invgamma_parms = function(lower_quantile,
                               upper_quantile,
                               response_variance,
                               lower_R2,
                               upper_R2,
                               probability,
                               plot_results = TRUE){
  tail_prob = 0.5 * (1.0 - probability)
  
  if(missing(lower_quantile) | 
     missing(upper_quantile)){
    
    if(missing(response_variance) |
       missing(lower_R2) |
       missing(upper_R2))
      stop("Must provide EITHER\nboth lower_quantile and upper_quantile, \nOR\nprovide response_variance, lower_R2, and upper_R2.")
    
    if( (lower_R2 < 0) | (lower_R2 >= 1) )
      stop("lower_R2 must be between 0 and 1")
    
    if( (upper_R2 <= 0) | (upper_R2 > 1) )
      stop("upper_R2 must be between 0 and 1")
    
    if(lower_R2 >= upper_R2)
      stop("lower_R2 must be less than upper_R2")
    
    
    lower_quantile = (1.0 - upper_R2) * response_variance
    upper_quantile = (1.0 - lower_R2) * response_variance
  }
  
  if(missing(probability))
    stop("Must provide the probability to be contained within the two quantiles")
  
  if(!missing(response_variance) & 
     !missing(lower_R2) & 
     !missing(upper_R2))
  
  helper = function(x){
    a = exp(x[1])
    b = exp(x[2])
    (pinvgamma(upper_quantile,a,b) - (1.0 - tail_prob))^2 + 
      (pinvgamma(lower_quantile,a,b) - tail_prob)^2
  }
  opt = optim(numeric(2),
              f = helper)
  a = exp(opt$par[1])
  b = exp(opt$par[2])
  
  if(plot_results){
    curve(dinvgamma(x,a,b),
          from= qinvgamma(0.01,a,b),
          to = qinvgamma(0.99,a,b),
          lwd = 3,
          yaxt = "n",
          xlab = expression(sigma^2))
  }
  
  
  return(c(shape = a, rate = b) * 2)
}

