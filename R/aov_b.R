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
#' @param ROPE numeric.  Used to compute posterior probability that Cohen's D +/- ROPE
#' @param contrasts numeric/matrix. Either vector of length equal to the number of 
#' levels in the grouping variable, or else a matrix where each row is a separate 
#' contrast.
#' @param n_joint_draws integer. Number of posterior draws to obtain.
#' @param improper logical.  Should we use an improper prior that is proportional 
#' to the inverse of the variance?
#' 
#' @return Object of class "aov_b" with the following elements:
#' \itemize{
#'  \item summary - data.frame giving the summary of the model parameters
#'  \item pairwise summary - data.frame giving the summary comparing all 
#'  factor level means
#'  \item contrasts (if provided) - list with named elements L (the contrasts provided 
#'  by the user) and summary.
#'  \item posterior_draws - mcmc object (see coda package) giving the posterior draws
#'  \item formula, data - input by user
#'  \item posterior_parameters - 
#'  \itemize{
#'    \item mu_g - the post. means of the group means
#'    \item nu_g - the post. scalars of the precision
#'    \item a_g - the post. shape of the inv. gamma for the group variances
#'    \item b_g - the post. rate of the inv. gamma for the group variances.
#'  }
#' }
#' 
#' @import magrittr
#' @import dplyr
#' @import coda
#' @import extraDistr
#' @export
#' @exportClass aov_b

aov_b = function(formula,
                 data,
                 heteroscedastic = TRUE,
                 prior_mean_mu = 0,
                 prior_mean_nu = 0.001,
                 prior_var_shape = 0.001,
                 prior_var_rate = 0.001,
                 CI_level = 0.95,
                 ROPE = 0.1,
                 contrasts,
                 n_joint_draws = 1e3,
                 improper = FALSE){
  
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
  
  # Check if improper prior \propto 1/\sigma^2 is requested
  if(improper){
    prior_mean_mu = 0.0
    prior_mean_nu = 0.0
    prior_var_shape = -1.0
    prior_var_rate = 0.0
  }
  
  # Get summary stats
  data_quants = 
    data %>% 
    group_by(group) %>% 
    summarize(n = n(),
              ybar = mean(y),
              y2 = sum(y^2),
              sample_var = var(y)) %>% 
    mutate(s2 = (n - 1) / n * sample_var)
  
  
  if(heteroscedastic){
    
    # Get posterior parameters
    nu_g = 
      prior_mean_nu + data_quants$n
    mu_g =
      (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
      nu_g
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
                 Upper = 0.0,
                 ROPE = 0.0,
                 EPR = 0.0)
    for(i in 1:nrow(ret$pairwise_summary)){
      ## Get CI for D(g,h)
      ret$pairwise_summary[i,c("Lower","Upper")] = 
        quantile(ret$posterior_draws[,temp[1,i]] - 
                   ret$posterior_draws[,temp[2,i]],
                 probs = c(a/2, 1 - a/2))
      
      # Get ROPE for D(g,h) based on Cohen's D
      ret$pairwise_summary$ROPE[i] = 
        mean(
          abs(mu_g_draws[,temp[1,i]] - 
                mu_g_draws[,temp[2,i]]) / 
            sqrt(
              (data_quants$n[temp[1,i]] * s2_g_draws[,temp[1,i]] + 
                 data_quants$n[temp[2,i]] * s2_g_draws[,temp[2,i]]) /
                sum(data_quants$n[temp[,i]])
            ) <
            ROPE
          )
      
      # Get EPR(g,h)
      ret$pairwise_summary$EPR[i] = 
        mean(
          1 - pnorm((mu_g_draws[,temp[1,i]] - mu_g_draws[,temp[2,i]]) / 
                      sqrt(s2_g_draws[,temp[1,i]] + s2_g_draws[,temp[2,i]]))
        )
      
    }
    colnames(ret$pairwise_summary)[ncol(ret$pairwise_summary) - 1] = 
      paste0("ROPE (",ROPE,")")
    
    
    
    # Compute contrasts if requested
    if(!missing(contrasts)){
      
      if("matrix" %in% class(contrasts)){
        L = tcrossprod(mu_g_draws, contrasts)
        rownames(contrasts) = paste("contrast",1:nrow(contrasts),sep="_")
      }else{
        L = mu_g_draws %*% contrasts
        contrasts = 
          matrix(contrasts,nr=1,dimnames = list("contrasts_1",NULL))
      }
      
      ret$contrasts = 
        list(L = contrasts,
             summary = data.frame(contrast = 1:nrow(contrasts),
                                  Estimate = colMeans(L),
                                  Lower = apply(L,2,quantile,probs = a/2),
                                  Upper = apply(L,2,quantile,probs = 1 - a/2)))
      
    }
    
    ret$formula = formula
    ret$data = data
    ret$posterior_parameters = 
      list(mu_g = mu_g,
           nu_g = nu_g,
           a_g = a_g,
           b_g = b_g)
    class(ret) = "aov_b"
    return(ret)
      
  }else{# start homoscedastic approach
    
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
    G = length(nu_g)
    
    # Return a summary including the posterior mean, credible intervals, and probability of direction
    ret = list()
    ret$summary = 
      data.frame(Variable = 
                   c(paste(rep("Mean",G),
                           rep(variables[2],G),
                           levels(data$group),
                           sep = " : "),
                     "Var"),
                 Estimate = c(mu_g, b_G/2 / (a_G/2 - 1.0)),
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
    ret$summary$ProbDir = 
      sapply(ret$summary$ProbDir, function(x) max(x,1-x))
    
    
    # Get posterior samples
    mu_g_draws = 
      matrix(0.0,
             n_joint_draws,
             G,
             dimnames = list(NULL,
                             paste("mean",levels(data$group),sep="_")))
    s2_G_draws =
      extraDistr::rinvgamma(n_joint_draws,
                            alpha = a_G/2,
                            beta = b_G/2)
    
    for(g in 1:G){
      mu_g_draws[,g] = 
        rnorm(n_joint_draws,
              mean = mu_g[g],
              sd = sqrt(s2_G_draws / nu_g[g]))
    }
    
    ret$posterior_draws = 
      cbind(mu_g_draws,
            Var = s2_G_draws) %>% 
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
                 Upper = 0.0,
                 ROPE = 0.0)
    for(i in 1:nrow(ret$pairwise_summary)){
      ## Get CI for D(g,h)
      ret$pairwise_summary[i,c("Lower","Upper")] = 
        quantile(ret$posterior_draws[,temp[1,i]] - 
                   ret$posterior_draws[,temp[2,i]],
                 probs = c(a/2, 1 - a/2))
      
      # Get ROPE for D(g,h) based on Cohen's D
      ret$pairwise_summary$ROPE[i] = 
        mean(
          abs(mu_g_draws[,temp[1,i]] - 
                mu_g_draws[,temp[2,i]]) / 
            sqrt(s2_G_draws) <
            ROPE
        )
      
      # Get EPR(g,h)
      ret$pairwise_summary$EPR[i] = 
        mean(
          1 - pnorm((mu_g_draws[,temp[1,i]] - mu_g_draws[,temp[2,i]]) / 
                      sqrt(2 * s2_G_draws))
        )
      
    }
    colnames(ret$pairwise_summary)[ncol(ret$pairwise_summary)] = 
      paste0("ROPE (",ROPE,")")
    
    
    # Compute contrasts if requested
    if(!missing(contrasts)){
      
      if("matrix" %in% class(contrasts)){
        L = tcrossprod(mu_g_draws, contrasts)
        rownames(contrasts) = paste("contrast",1:nrow(contrasts),sep="_")
      }else{
        L = mu_g_draws %*% contrasts
        contrasts = 
          matrix(contrasts,nr=1,dimnames = list("contrasts_1",NULL))
      }
      
      ret$contrasts = 
        list(L = contrasts,
             summary = data.frame(contrast = 1:nrow(contrasts),
                                  Estimate = colMeans(L),
                                  Lower = apply(L,2,quantile,probs = a/2),
                                  Upper = apply(L,2,quantile,probs = 1 - a/2)))
      
    }
    
    ret$formula = formula
    ret$data = data
    ret$posterior_parameters = 
      list(mu_g = mu_g,
           nu_g = nu_g,
           a_g = a_G,
           b_g = b_G)
    class(ret) = "aov_b"
    return(ret)
    
  }
  
}
