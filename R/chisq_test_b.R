#' @name chisq_test_b
#' @aliases homogeneity_b
#' @aliases independence_b
#' 
#' Homogeneity independence analyses for count data
#' 
#' @details
#' 
#' \strong{Homogeneity Analysis}
#' For a 2-way contingency table with 2 rows and C columns, compare the 
#' C-dimensional probability vectors corresponding to each row.
#' 
#' \strong{Indepdence analysis}
#' For a 2-way contingency table with R rows and C columns, evaluate 
#' the probability that the joint probabilities \eqn{p_{rc}} are all 
#' within the ROPE of \eqn{p_{r\cdot}\times p_{\cdot c}}.
#' 
#' 
#' @import extraDistr

#' @export
homogeneity_b = function(x,
                         y,
                         ROPE,
                         prior = "uniform",
                         prior_shapes,
                         CI_level = 0.95,
                         mc_relative_error =
                         seed = 1,
                         mc_relative_error = 0.01){
  
  # Get data
  x = as.matrix(x)
  if(dim(x)[2] == 1){
    x = as.vector(x)
    if(missing(y))
      stop("if x is a vector, y must also be supplied")
    if (length(x) != length(y)) 
      stop("'x' and 'y' must have the same length")
    x = rbind(x,c(y))
  }
  if(nrow(x) != 2)
    stop("x must have two rows.")
  
  J = ncol(x)
  alpha_ci = 1.0 - CI_level
  
  
  # Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "jeffreys")[pmatch(tolower(prior),
                                 c("uniform",
                                   "jeffreys"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,J)
    }
    if(prior == "jeffreys"){
      message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
      prior_shapes = rep(0.5,J)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
    if(length(prior_shapes) == 1)
      prior_shapes == rep(prior_shapes,J)
    if(length(prior_shapes) != J)
      stop("Length of prior_shapes should either match the number of categories or be of length 1")
  }
  
  # Get ROPE
  if(missing(ROPE)){
    ROPE = c(1.0 / 1.125, 1.125)
    # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
    #   Use the same thing for odds ratios.
  }else{
    if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
    if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
    if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
  }
  
  # Get posterior parameters
  post_shapes = 
    tcrossprod(matrix(1.0,2,1), prior_shapes) + x
  
  # Get posterior draws
  set.seed(seed)
  ## Get preliminary draws
  p1_draws = 
    rdirichlet(500,
               post_shapes[1,])
  p2_draws = 
    rdirichlet(500,
               post_shapes[2,])
  ## Use CLT for empirical quantiles:
  #     A Central Limit Theorem For Empirical Quantiles in the Markov Chain Setting. Peter W. Glynn and Shane G. Henderson
  #     With prob 0.99 we will be within mc_relative_error of the alpha_ci/2 quantile
  fhat = 
    lapply(1:J,
           function(j){
             density(p1_draws[,j] - p2_draws[,j],
                     from = -1.0 + .Machine$double.eps,
                     to = 1.0 - .Machine$double.eps)
           })
    
  n_draws = 
    sapply(1:J,
           function(j){
             0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
               (
                 qnorm(0.5 * (1.0 - 0.99)) / 
                   mc_relative_error /
                   quantile(p1_draws[,j] - p2_draws[,j], 0.5 * alpha_ci) /
                   fhat[[j]]$y[which.min(abs(fhat[[j]]$x - 
                                          quantile(p1_draws[,j] - p2_draws[,j], 0.5 * alpha_ci)))]
               )^2 |> 
               round()
           })
    
  
  
  
}


