#' Bayesian Wilcoxon Rank Sum and Signed Rank Analyses
#' 
#' 
#' @param x numeric vector of data values. Non-finite (e.g., infinite or 
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite values will be omitted.
#' @param p 
#' \itemize{
#'  \item Signed rank: If provided, \code{wilcox_test_b} will return the 
#' posterior probability that the population proportion is less than this value.
#'  \item Rank sum: 
#' }
#' @param ROPE 
#' \itemize{
#'  \item Signed rank: If a single number, ROPE will be \code{p}\eqn{\pm}\code{ROPE}. 
#'  If a vector of length 2, these will serve as the ROPE bounds. Defaults to 
#'  \eqn{\pm 0.05}.
#'  \item Rank sum: 
#' }
#' @param prior Prior used on the probability that x > y.  Either 
#' "jeffreys" (Beta(1/2,1/2)), "uniform" (Beta(1,1)), or "centered" (Beta(2,2)). 
#' This is ignored if prior_shapes is provided.
#' @param prior_shapes Vector of length two, giving the shape parameters for the 
#' beta distribution that will act as the prior on the population proportions.
#' @param CI_level The posterior probability to be contained in the credible interval.
#' @param plot logical.  Should a plot be shown?
#' @param seed Always set your seed! (Unused for \eqn{\geq} 20 observations.)
#' 
#' @import DFBA
#' @export

wilcox_test_b = function(x,
                         y,
                         paired = FALSE,
                         p = 0.5,
                         ROPE,
                         prior = "centered",
                         prior_shapes,
                         CI_level = 0.95,
                         plot = TRUE,
                         seed = 1){
  if( !("numeric" %in% class(x)) ) 
    stop("x must be numeric")
  if( !missing(y) && !("numeric" %in% class(y)) )
    stop("y must be numeric")
  
  
  if(paired | missing(y)){
    
    # Do signed rank test
    
    ## Get outcome vector of differences
    if(!missing(y)){
      x = na.omit(x); y = na.omit(y)
      if(length(x) != length(y))
        stop("If paired = TRUE, then x and y must be of the same length.")
      x = x - y
    }
    
    ## Get ROPE
    if(missing(ROPE)){
      ROPE = 
        min(c(0.05,
              p / 2,
              (1.0 - p) / 2))
      ROPE_bounds = p + c(-1,1) * ROPE
    }else{
      if(!(length(ROPE) %in% 1:2))
        stop("ROPE must be of length 1 or 2")
      if(length(ROPE) == 1){
        if(ROPE < 0)
          stop("ROPE must be positive")
        if(ROPE > min(c(p,1.0 - p) / 2))
          stop("ROPE extends from p0 beyond 0 or 1")
        ROPE_bounds = p + c(-1,1) * ROPE
      }else{
        if(ROPE[1] < 0)
          stop("Lower bound of ROPE cannot be less than 0")
        if(ROPE[2] < 0)
          stop("Upper bound of ROPE cannot be greater than 1")
        ROPE_bounds = ROPE
      }
    }
    
    ## Prior distribution
    if(missing(prior_shapes)){
      prior = c("uniform",
                "jeffreys",
                "centered")[pmatch(tolower(prior),
                                   c("uniform",
                                     "jeffreys",
                                     "centered"))]
      
      if(prior == "uniform"){
        message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
        prior_shapes = rep(1.0,2)
      }
      if(prior == "jeffreys"){
        message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
        prior_shapes = rep(0.5,2)
      }
      if(prior == "centered"){
        message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
        prior_shapes = rep(2.0,2)
      }
    }else{
      if(any(prior_shapes <= 0))
        stop("Prior shape parameters must be positive.")
    }
    
    
    ## Use DFBA package to do analysis
    dfba_object = 
      DFBA::dfba_wilcoxon(Y1 = x,
                          Y2 = numeric(length(x)),
                          a0 = prior_shapes[1],
                          b0 = prior_shapes[2],
                          prob_interval = CI_level)
    
    # Compute results
    results = 
      list(posterior_mean = 
             dfba_object$a_post / 
             sum( dfba_object$a_post + dfba_object$b0),
           CI = 
             c(dfba_object$eti_lower,
               dfba_object$eti_upper),
           Pr_less_than_p = 
             pbeta(p,
                   dfba_object$a_post,
                   dfba_object$b_post),
           Pr_in_ROPE = 
             pbeta(ROPE_bounds[2],
                   dfba_object$a_post,
                   dfba_object$b_post) - 
             pbeta(ROPE_bounds[1],
                   dfba_object$a_post,
                   dfba_object$b_post)
      )
    
    # Print results
    cat("\n----------\n\nWilcoxon signed-rank analysis using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: ", 
               format(signif(results$posterior_mean, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: (", 
               format(signif(results$CI[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that Pr(x > y) > ",
               format(signif(p, 3), 
                      scientific = FALSE),
               " = ",
               format(signif(1.0 - results$Pr_less_than_p, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0("Probability that Pr(x > y) is in the ROPE, defined to be (",
               format(signif(ROPE_bounds[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE_bounds[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      
      results$prob_plot = 
        tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                       l = 50)) |> 
        ggplot(aes(x=x)) +
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                prior_shapes[1],
                                prior_shapes[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                dfba_object$a_post,
                                dfba_object$b_post)
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      print(results$prob_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = 
      c(shape_1 = dfba_object$a_post,
        shape_2 = dfba_object$b_post)
    
    invisible(results)
    
  }else{#End: Wilcoxon signed rank analysis
    # Do rank sum test
    
    
    
  }
  
}


.wilcox_signed_rank = function(x,
                               y,
                               asdf){
  
}

.wilcox_rank_sum = function(x,y,
                            asdf){
  
}