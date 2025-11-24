#' t-test
#' 
#' One and two sample t-tests on vectors of data
#' 
#' @param x Either a (non-empty) numeric vector of data values, or a formula 
#' of the form outcome ~ grouping variable.
#' @param y an optional (non-empty) numeric vector of data values
#' @param mu optional.  If supplied, \code{t_test_b} will return the 
#' posterior probabilty that the population mean (ignored in 2 sample inference) 
#' is less than this value.
#' @param paired logical.  If TRUE, provide both x and y as vectors.
#' @param data logical.  Only used if x is a formula.
#' @param heterscedastic logical.  Set to FALSE to assume all groups have 
#' equal variance.
#' @param n_joint_draws non-negative integer.  Number of posterior draws of the 
#' variance(s) and group means.
#' @param prior_mean_mu numeric. Hyperparameter for the a priori mean of the 
#' group means.
#' @param prior_mean_nu numeric. Hyperparameter which scales the precision of 
#' the group means.
#' @param prior_var_shape numeric. Twice the shape parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param prior_var_rate  numeric. Twice the rate parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param CI_level numeric. Credible interval level.
#' @param ROPE numeric.  Used to compute posterior probability that Cohen's D +/- ROPE
#' @param contrasts numeric/matrix. Either vector of length equal to the number of 
#' levels in the grouping variable, or else a matrix where each row is a separate 
#' contrast.
#' @param n_joint_draws integer. Number of posterior draws to obtain.
#' @param improper logical.  Should we use an improper prior that is proportional 
#' to the inverse of the variance?
#' @param seed integer.  Always set your seed!!!
#' 
#' @returns Either an aov_b object, if two samples are being compared,
#' or a list with the following elements:
#' \itemize{
#'  \item Variable
#'  \item Post Mean
#'  \item Lower (bound of credible interval)
#'  \item Upper (bound of credible interval)
#'  \item ProbDir (Probability of Direction)
#' }
#' 
#' @export


t_test_b = function(x,
                    y,
                    mu,
                    paired = FALSE,
                    data,
                    heteroscedastic = TRUE,
                    prior_mean_mu = 0,
                    prior_mean_nu = 0.001,
                    prior_var_shape = 0.001,
                    prior_var_rate = 0.001,
                    CI_level = 0.95,
                    ROPE = 0.1,
                    contrasts,
                    improper = FALSE,
                    seed = 1,
                    mc_relative_error = 0.01){
  
  outcome_name = NULL
  if((class(x) == "formula")){
    outcome_name = all.vars(x)[1]
    if(length(all.vars(x)) == 1){#Intercept only model, i.e., one population
      x = data[[outcome_name]]
    }
  }
  
  if((class(x) == "formula") & missing(data)) stop("If formula is given, data must also be given.")
  if(missing(y) & paired) stop("Cannot have paired data without y.")
  
  if(is.numeric(x)){
    if(missing(y) | paired){
      if(!missing(y) && paired && (length(x) != length(y)) ) stop("Length of x must equal that of y.")
      if(paired){
        x = x - y
        outcome_name = "x minus y"
      }
      
      if(is.null(outcome_name)) outcome_name = "x"
      
      # Set alpha lv
      a = 1 - CI_level
      
      # Check if improper prior \propto 1/\sigma^2 is requested
      if(improper){
        prior_mean_mu = 0.0
        prior_mean_nu = 0.0
        prior_var_shape = -1.0
        prior_var_rate = 0.0
      }
      
      # Get summary stats
      data_quants = 
        tibble(n = NROW(x),
               ybar = mean(x),
               y2 = sum(x^2),
               sample_var = var(x)) |> 
        mutate(s2 = (n - 1.0) / n * sample_var)
      
      # Get posterior parameters
      nu_g = 
        prior_mean_nu + data_quants$n
      mu_g =
        (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
        nu_g
      a_G =
        prior_var_shape + sum(data_quants$n)
      b_G =
        prior_var_rate +
        sum(
          data_quants$n * data_quants$s2 +
            prior_mean_nu * data_quants$n / (nu_g + data_quants$n) * (prior_mean_mu - data_quants$ybar)^2
        )
      
      # Return a summary including the posterior mean, credible intervals, and probability of direction
      ret = 
        tibble(Variable = 
                 c(outcome_name,
                   "Var"),
               `Post Mean` = c(mu_g, b_G/2 / (a_G/2 - 1.0)),
               Lower = c(extraDistr::qlst(a/2, 
                                          df = a_G,
                                          mu = mu_g,
                                          sigma = sqrt(b_G / nu_g / a_G)),
                         extraDistr::qinvgamma(a/2, alpha = a_G/2, beta = b_G/2)),
               Upper = c(extraDistr::qlst(1 - a/2, 
                                          df = a_G,
                                          mu = mu_g,
                                          sigma = sqrt(b_G / nu_g / a_G)),
                         extraDistr::qinvgamma(1 - a/2, alpha = a_G/2, beta = b_G/2)),
               ProbDir = c(extraDistr::plst(0, 
                                            df = a_G,
                                            mu = mu_g,
                                            sigma = sqrt(b_G / nu_g / a_G)),
                           NA))
      ret$ProbDir = 
        sapply(ret$ProbDir, function(x) max(x,1-x))
      
      return(ret)
      
    }else{#End: one sample inference
      
      ttest_data = 
        tibble(group = rep(c("x","y"),
                           c(length(x),
                             length(y)))) %>% 
        mutate(y = c(x,y))
      
      if(missing(prior_mean_mu)){
        return(aov_b(y ~ group,
                     data = ttest_data,
                     heteroscedastic = TRUE,
                     prior_mean_nu = prior_mean_nu,
                     prior_var_shape = prior_var_shape,
                     prior_var_rate = prior_var_rate,
                     CI_level = CI_level,
                     ROPE = ROPE,
                     improper = improper,
                     seed = seed,
                     mc_relative_error = mc_relative_error))
      }else{
        return(aov_b(y ~ group,
                     data = ttest_data,
                     heteroscedastic = TRUE,
                     prior_mean_mu = prior_mean_mu,
                     prior_mean_nu = prior_mean_nu,
                     prior_var_shape = prior_var_shape,
                     prior_var_rate = prior_var_rate,
                     CI_level = CI_level,
                     ROPE = ROPE,
                     improper = improper,
                     seed = seed,
                     mc_relative_error = mc_relative_error))
      }
      
    }
  }else{
    if(missing(prior_mean_mu)){
      return(aov_b(x,
                   data = data,
                   heteroscedastic = TRUE,
                   prior_mean_nu = prior_mean_nu,
                   prior_var_shape = prior_var_shape,
                   prior_var_rate = prior_var_rate,
                   CI_level = CI_level,
                   ROPE = ROPE,
                   improper = improper,
                   seed = seed,
                   mc_relative_error = mc_relative_error))
    }else{
      return(aov_b(x,
                   data = data,
                   heteroscedastic = TRUE,
                   prior_mean_mu = prior_mean_mu,
                   prior_mean_nu = prior_mean_nu,
                   prior_var_shape = prior_var_shape,
                   prior_var_rate = prior_var_rate,
                   CI_level = CI_level,
                   ROPE = ROPE,
                   improper = improper,
                   seed = seed,
                   mc_relative_error = mc_relative_error))
    }
  }
  
}