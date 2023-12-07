#' Analysis of Variance using Bayesian methods
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.
#' @param heterscedastic logical.  Set to FALSE to assume all groups have 
#' equal variance.
#' @param n_joint_draws non-negative integer.  Number of posterior draws of the 
#' variance(s) and group means.
#' @param prior_mean_mu numeric. Hyperparameter for the a priori mean of the 
#' group means.
#' @param prior_mean_nu numeric. Hyperparameter which scales the precision of 
#' the group means.
#' @param prior_var_shape numeric. Shape parameter for the inverse gamma prior on
#' the variance(s).
#' @param prior_var_rate  numeric. Rate parameter for the inverse gamma prior on
#' the variance(s).
#' @param CI_level numeric. Credible interval level.
#' 
#' @return asdfasdfa
#' 
#' @import magrittr
#' @import dplyr
#' @import coda
#' @import extraDistr

if(FALSE){
  library(magrittr)
  library(dplyr)
  library(coda)
  library(extraDistr)
  formula = mpg ~ gear
  data = mtcars# %>% mutate(gear = factor(gear))
  prior_mean_mu = 15
  prior_mean_nu = 0.5
  prior_var_shape = 0.01
  prior_var_rate = 0.01
  CI_level = 0.95
  
}

aov_b = function(formula,
                 data,
                 heteroscedastic = TRUE,
                 prior_mean_mu,
                 prior_mean_nu,
                 prior_var_shape,
                 prior_var_rate,
                 CI_level = 0.95,
                 ROPE = 0.1,
                 n_joint_draws = 1e3){
  
  # Set alpha lv
  a = 1 - CI_level
  
  # Extract variable names
  variables = all.vars(formula)
  if(!("factor" %in% class(data[[variables[2]]]))){
    data[[variables[2]]] = 
      factor(data[[variables[2]]])
  }
  data %<>%
    dplyr::rename(y = variables[1],
                  group = variables[2])
  
  
  
  if(heteroscedastic){
    
    # Get summary stats
    data_quants = 
      data %>% 
      group_by(group) %>% 
      summarize(n = n(),
                ybar = mean(y),
                y2 = sum(y^2),
                sample_var = var(y)) %>% 
      mutate(s2 = (n - 1) / n * sample_var)
    
    # Get posterior parameters
    nu_g = 
      prior_mean_nu + data_quants$n
    mu_g =
      (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
      (nu_g)
    a_g =
      prior_var_shape + data_quants$n
    b_g =
      prior_var_rate +
      data_quants$n * data_quants$s2 +
      prior_mean_nu * data_quants$n / (nu_g + data_quants$n) * (prior_mean_mu - data_quants$ybar)^2
    G = length(nu_g)
    
    # Return a summary including the posterior mean, credible intervals, and probability of direction
    ret = list()
    ret$summary = 
      data.frame(Variable = 
                   paste(rep(c("Mean","Var"),each = G),
                         rep(variables[2],2*G),
                         rep(levels(data$group),2),
                           sep = " : "),
                 Estimate = c(mu_g, b_g/2 / (a_g/2 - 1.0)),
                 Lower = c(extraDistr::qlst(a/2, 
                                            df = a_g,
                                            mu = mu_g,
                                            sigma = sqrt(b_g / nu_g / a_g)),
                           extraDistr::qinvgamma(a/2, alpha = a_g/2, beta = b_g/2)),
                 Upper = c(extraDistr::qlst(1 - a/2, 
                                            df = a_g,
                                            mu = mu_g,
                                            sigma = sqrt(b_g / nu_g / a_g)),
                           extraDistr::qinvgamma(1 - a/2, alpha = a_g/2, beta = b_g/2)),
                 ProbDir = c(extraDistr::plst(0, 
                                              df = a_g,
                                              mu = mu_g,
                                              sigma = sqrt(b_g / nu_g / a_g)),
                             rep(NA,G)))
    ret$summary$ProbDir = 
      sapply(ret$summary$ProbDir, function(x) max(x,1-x))
    
    
    # Get posterior samples
    mu_g_draws = 
      matrix(0.0,
             n_joint_draws,
             G,
             dimnames = list(NULL,
                             paste("mean",levels(data$group),sep="_")))
    s2_g_draws = 
      matrix(0.0,
             n_joint_draws,
             G,
             dimnames = list(NULL,
                             paste("variance",levels(data$group),sep="_")))
             
             
    for(g in 1:G){
      s2_g_draws[,g] =
        extraDistr::rinvgamma(n_joint_draws,
                              alpha = a_g[g]/2,
                              beta = b_g[g]/2)
      mu_g_draws[,g] = 
        rnorm(n_joint_draws,
              mean = mu_g[g],
              sd = sqrt(s2_g_draws[,g] / nu_g[g]))
    }
    
    ret$posterior_draws = 
      cbind(mu_g_draws,
            s2_g_draws) %>% 
      as.mcmc()
    
    # Get pairwise comparisons
    temp = 
      combn(1:length(levels(data$group)),2)
    ret$pairwise_summary = 
      data.frame(Comparison = 
                   apply(combn(levels(data$group),2),
                         2,
                         function(x) paste(x[1],x[2],sep="-")),
                 Estimate = mu_g[temp[1,]] - mu_g[temp[2,]],
                 Lower = 0.0,
                 Upper = 0.0)
    for(i in 1:nrow(ret$pairwise_summary)){
      ret$pairwise_summary[i,c("Lower","Upper")] = 
        quantile(ret$posterior_draws[,temp[1,i]] - 
                   ret$posterior_draws[,temp[2,i]],
                 probs = c(a/2, 1 - a/2))
      
      
      
    }
    
    
    
    return(ret)
      
  }else{# start homoscedastic approach
    
    
  }
  
}
