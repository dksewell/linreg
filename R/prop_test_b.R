#' Bayesian test of Equal or Given Proportions
#' 
#' 


prop_test_b = function(n_successes,
                       n_failures,
                       n_total,
                       p,
                       predict_for_n,
                       ROPE,
                       prior = "jeffreys",
                       prior_shapes,
                       CI_level = 0.95,
                       PI_level = 0.95,
                       plot = TRUE,
                       seed = 1,
                       n_draws = 5e4){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  # Preliminary checks
  if(missing(n_successes)) 
    stop("Must provide the number of successes")
  if(missing(n_failures) & missing(n_total)) 
    stop("Must provide either number of failures or the total number of trials.")
  if( !missing(n_failures) & !missing(n_total) ) 
    if(!isTRUE(all.equal(n_failures,n_total - n_successes)))
      stop("Number of failures provided is inconsistent with total number of trials")
  if(missing(n_failures))
    n_failures = n_total - n_successes
  if(missing(n_total))
    n_total = n_failures + n_successes
  if(missing(predict_for_n)){
    predict_for_n = n_total
  }
  if(length(n_successes) != length(n_failures))
    stop("Length of n_successes must be the same length as n_successes/n_failures")
  if(length(n_successes) > 2)
    stop("The number of populations must be either 1 or 2.")
  
  
  # Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
             "jeffreys")[pmatch(tolower(prior),
                               c("uniform",
                                 "jeffreys"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "jeffreys"){
      message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
      prior_shapes = rep(0.5,2)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  # One sample inference 
  if(length(n_successes) == 1){
    
    ## Get posterior parameters
    post_shapes = 
      prior_shapes + 
      c(n_successes,
        n_failures)
    
    ## Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      pbbinom(0:predict_for_n,
              predict_for_n,
              post_shapes[1],
              post_shapes[2])
    
    ## Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean = 
             post_shapes[1] / sum(post_shapes),
           CI_lower = 
             qbeta(0.5 * alpha_ci,
                   post_shapes[1],
                   post_shapes[2]),
           CI_upper = 
             qbeta(1.0 - 0.5 * alpha_ci,
                   post_shapes[1],
                   post_shapes[2]),
           PI_lower = 
             max(which(cdf_probs <= 0.5 * alpha_pi)),
           PI_upper = 
             min(predict_for_n,
                 1 + min(which(cdf_probs > 1.0 - 0.5 * alpha_pi)))
      )
    
    ## Print results
    cat("\n----------\n\nAnalysis of a single population proportion using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: ", n_successes,"\n\n"))
    cat(paste0("Number of failures: ", n_failures,"\n\n"))
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
               format(signif(results$CI_lower, 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_upper, 3), 
                      scientific = FALSE),
               ")\n\n"))
    if(!missing(p)){
      results$Pr_less_than_p = 
        pbeta(p,post_shapes[1],post_shapes[2])
      cat(paste0("Probability that p < ",
                 format(signif(p, 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(results$Pr_less_than_p, 3), 
                        scientific = FALSE),
                 "\n\n"))
    }
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n == n_total,
                      "another ",
                      ""),
               predict_for_n,
               " trials: (", 
               results$PI_lower,
               ", ",
               results$PI_upper,
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    ## Plot (if requested)
    if(plot){
      results$prop_plot = 
        tibble(x = seq(.Machine$double.eps,1.0 - .Machine$double.eps,
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
                                post_shapes[1],
                                post_shapes[2])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        ylim(0,
             ifelse(sum(post_shapes > 1) == 2,
                    1.05 * dbeta((post_shapes[1] - 1) / (sum(post_shapes) - 2),
                                 post_shapes[1],
                                 post_shapes[2]),
                    2)) + 
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      results$prop_plot
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = 
      c(shape1 = post_shapes[1],
        shape2 = post_shapes[2])
    
    invisible(results)
  }else{#End: One sample inference
    # Two sample inference
    
    ## Get posterior parameters
    post_shapes = 
      rbind(prior_shapes,
            prior_shapes) + 
      cbind(n_successes,
            n_failures)
    
    ## Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      list(pop1 = 
             pbbinom(0:predict_for_n[1],
                     predict_for_n,
                     post_shapes[1,1],
                     post_shapes[1,2]),
           pop2 = 
             pbbinom(0:predict_for_n[2],
                     predict_for_n,
                     post_shapes[2,1],
                     post_shapes[2,2])
      )
    
    ## Compute post prob that p1 < p2
    # joint_post = function(p1,p2){
    #   dbeta(p1,
    #         post_shapes[1,1],
    #         post_shapes[1,2]) * 
    #     dbeta(p2,
    #           post_shapes[2,1],
    #           post_shapes[2,2])
    # }
    # # cubature package seems a bit slow...
    # # prob_p1_less_p2 = function(x){
    # #   if(x[1] < x[2]){
    # #     joint_post(x[1],x[2])
    # #   }else{
    # #     return(0.0)
    # #   }
    # # }
    # # 
    # # test1 = 
    # #   cubature::adaptIntegrate(prob_p1_less_p2,
    # #                          lowerLimit = c(0.0,0.0),
    # #                          upperLimit = c(1.0,1.0))
    # results$prob_p1_less_p2 = 
    #   integrate(function(p2) {
    #     inner = sapply(p2,function(x){
    #       integrate(function(p1) joint_post(p1, x),
    #                        lower = 0, upper = x)$value
    #     })
    #     inner
    #   }, lower = 0, upper = 1)$value
    p1_draws = 
      rbeta(n_draws,
            post_shapes[1,1],
            post_shapes[1,2])
    p2_draws = 
      rbeta(n_draws,
            post_shapes[2,1],
            post_shapes[2,2])
    
    results$prob_p1_less_p2 = 
      mean(p1_draws < p2_draws)
    
    ## Find CI for difference in prob
    # fun_to_optim = function(bounds){
    #   abs(
    #     bounds[1] - 
    #       integrate(function(p2) {
    #         inner = sapply(p2,function(x){
    #           integrate(function(p1) joint_post(p1, x),
    #                     lower = max(0.0,),
    #                     upper = x)$value
    #         })
    #         inner
    #       }, lower = 0, upper = 1)$value
    #   ) + 
    #     abs(
    #       asdf
    #     )
    # }
    CI_bounds = 
      quantile(p1_draws - p2_draws,c(0.5 * alpha_ci,
                                     1.0 - 0.5 * alpha_ci))
    results$CI_lower_diff_in_props = 
      CI_bounds[1]
    results$CI_upper_diff_in_props = 
      CI_bounds[2]
    stop("left off here")
    
           
    ## Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean = 
             post_shapes[1] / sum(post_shapes),
           CI_lower_pop1 = 
             qbeta(0.5 * alpha_ci,
                   post_shapes[1,1],
                   post_shapes[2,1]),
           CI_upper_pop1 = 
             qbeta(1.0 - 0.5 * alpha_ci,
                   post_shapes[1,1],
                   post_shapes[2,1]),
           PI_lower_pop1 = 
             max(which(cdf_probs$pop1 <= 0.5 * alpha_pi)),
           PI_upper_pop1 = 
             min(predict_for_n,
                 1 + min(which(cdf_probs$pop1 > 1.0 - 0.5 * alpha_pi))),
           CI_lower_pop2 = 
             qbeta(0.5 * alpha_ci,
                   post_shapes[1,2],
                   post_shapes[2,2]),
           CI_upper_pop2 = 
             qbeta(1.0 - 0.5 * alpha_ci,
                   post_shapes[1,2],
                   post_shapes[2,2]),
           PI_lower_pop2 = 
             max(which(cdf_probs$pop2 <= 0.5 * alpha_ci)),
           PI_upper_pop2 = 
             min(predict_for_n,
                 1 + min(which(cdf_probs$pop2 > 1.0 - 0.5 * alpha_ci)))
      )
    
    ## Print results
    cat("\n----------\n\nAnalysis of a single population proportion using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: ", n_successes,"\n\n"))
    cat(paste0("Number of failures: ", n_failures,"\n\n"))
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
               format(signif(results$CI_lower, 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_upper, 3), 
                      scientific = FALSE),
               ")\n\n"))
    if(!missing(p)){
      results$Pr_less_than_p = 
        pbeta(p,post_shapes[1],post_shapes[2])
      cat(paste0("Probability that p < ",
                 format(signif(p, 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(results$Pr_less_than_p, 3), 
                        scientific = FALSE),
                 "\n\n"))
    }
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n == n_total,
                      "another ",
                      ""),
               predict_for_n,
               " trials: (", 
               results$PI_lower,
               ", ",
               results$PI_upper,
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    ## Plot (if requested)
    if(plot){
      results$prop_plot = 
        tibble(x = seq(.Machine$double.eps,1.0 - .Machine$double.eps,
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
                                post_shapes[1],
                                post_shapes[2])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        ylim(0,
             ifelse(sum(post_shapes > 1) == 2,
                    1.05 * dbeta((post_shapes[1] - 1) / (sum(post_shapes) - 2),
                                 post_shapes[1],
                                 post_shapes[2]),
                    2)) + 
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      results$prop_plot
      
    }
    
    
    invisible(results)
  }#End: 2 sample inference
}

