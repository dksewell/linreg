#' @name summary
#' 
#' @title Summary functions for linreg objects
#' 
#' @param object linreg object
#' @param CI_level Posterior probability covered by credible interval
#' 
#' 
#' @export

#' @rdname summary
#' @export
summary.lm_b = function(object,
                        CI_level = 0.95){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  if(object$prior != "improper"){
    summ$Lower = 
      qlst(alpha/2,
           object$posterior_parameters$a_tilde,
           object$posterior_parameters$mu_tilde,
           sqrt(object$posterior_parameters$b_tilde/object$posterior_parameters$a_tilde * 
                  diag(qr.solve(object$posterior_parameters$V_tilde))))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           object$posterior_parameters$a_tilde,
           object$posterior_parameters$mu_tilde,
           sqrt(object$posterior_parameters$b_tilde/object$posterior_parameters$a_tilde * 
                  diag(qr.solve(object$posterior_parameters$V_tilde))))
  }else{
    summ$Lower = 
      qlst(alpha/2,
           nrow(object$data) - length(object$posterior_parameters$mu_tilde),
           object$posterior_parameters$mu_tilde,
           sqrt(diag(object$posterior_parameters$Sigma)))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           nrow(object$data) - length(object$posterior_parameters$mu_tilde),
           object$posterior_parameters$mu_tilde,
           sqrt(diag(object$posterior_parameters$Sigma)))
  }
  
  summ
}

#' @rdname summary
#' @export
summary.aov_b = function(object,
                         CI_level = 0.95){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  pw_summ = 
    object$pairwise_summary %>% 
    as.data.frame()
  
  
  cat("\n--- Summary of factor level means ---\n")
  summ$Lower = 
    c(extraDistr::qlst(alpha/2, 
                       df = object$posterior_parameters$a_g,
                       mu = object$posterior_parameters$mu_g,
                       sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
      extraDistr::qinvgamma(alpha/2, 
                            alpha = object$posterior_parameters$a_g/2, 
                            beta = object$posterior_parameters$b_g/2))
  summ$Upper = 
    c(extraDistr::qlst(1 - alpha/2, 
                       df = object$posterior_parameters$a_g,
                       mu = object$posterior_parameters$mu_g,
                       sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
      extraDistr::qinvgamma(1 - alpha/2, 
                            alpha = object$posterior_parameters$a_g/2, 
                            beta = object$posterior_parameters$b_g/2))
  print(summ)
  
  
  cat("\n\n\n\n--- Summary of pairwise differences ---\n")
  temp = 
    combn(1:length(levels(object$data$group)),2)
  for(i in 1:nrow(pw_summ)){
    pw_summ[i,c("Lower","Upper")] = 
      quantile(object$posterior_draws[,temp[1,i]] - 
                 object$posterior_draws[,temp[2,i]],
               probs = c(alpha/2, 
                         1 - alpha/2))
  }
  pw_summ %<>% as_tibble()
  print(pw_summ)
  cat("\n\n   *Note: EPR (Exceedence in Pairs Rate) for a Comparison of g-h = Pr(Y_(gi) > Y_(hi)|parameters) ")
  
  if(is.null(object$contrasts)){
    invisible(list(summary = object$summary,
                   pairwise = object$pairwise_summary))
  }else{
    
    cat("\n\n\n\n--- Summary of Contrasts ---\n")
    csumm = object$contrasts$summary
    contrast_draws = 
      tcrossprod(object$posterior_draws[,grep("mean_",colnames(object$posterior_draws))],
                 object$contrasts$L)
    csumm$Lower = 
      apply(contrast_draws,2,quantile,probs = alpha/2)
    csumm$Upper = 
      apply(contrast_draws,2,quantile,probs = 1 - alpha/2)
    
    print(csumm)
    
    invisible(list(summary = object$summary,
                   pairwise = object$pairwise_summary,
                   contrasts = 
                     list(L = object$contrasts$L,
                          summary = csumm)))
  }
}

#' @rdname summary
#' @export
summary.np_lm_b = function(object,
                           CI_level = 0.95){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  if("posterior_covariance" %in% names(object)){
    summ$Lower = 
      qnorm(alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
    summ$Upper = 
      qnorm(1 - alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
  }else{
    summ$Lower = 
      object$posterior_draws %>% 
      apply(2,quantile,prob = alpha / 2)
    summ$Upper = 
      object$posterior_draws %>% 
      apply(2,quantile,prob = 1.0 - alpha / 2)
  }
  
  summ
}

#' @rdname summary
#' @export
summary.lm_b_bma = function(object,
                       CI_level = 0.95,
                       interpretable_scale = TRUE){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  summ$Lower = 
    apply(object$posterior_draws,2,quantile,probs = alpha/2)
  summ$Upper =
    apply(object$posterior_draws,2,quantile,probs = 1.0 - alpha/2)
  
  
  
  # Exponentiate
  if(( (object$family$family == "binomial") & 
       (object$family$link != "logit") ) | 
     ( (object$family$family == "poisson") & 
       (object$family$link != "log") ) | 
     (object$family$familt == "gaussian") ){
    interpretable_scale = FALSE
  }
  if(interpretable_scale){
    paste0("\n----------\n\nValues given in terms of ",
           ifelse(object$family$family == "binomial",
                  "odds ratios",
                  "rate ratios")
    ) %>% 
      cat()
    cat("\n\n----------\n\n")
    summ = summ[-1,]
    summ[,c("Post Mean","Lower","Upper")] %<>%
      exp()
    summ[,"ROPE bounds"] = 
      paste("(",
            round(exp(-object$ROPE[-1]),3),
            ",",
            round(exp(object$ROPE[-1]),3),
            ")",
            sep="")
  }
  
  summ
}


#' @rdname summary
#' @export
summary.glm_b = function(object,
                         CI_level = 0.95,
                         interpretable_scale = TRUE){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  if("posterior_covariance" %in% names(object)){
    summ$Lower = 
      qnorm(alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
    summ$Upper = 
      qnorm(1 - alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
  }else{
    # Get CI bounds
    CI_from_weighted_sample = function(x,w){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      LB = max(which(w <= 0.5 * alpha))
      UB = min(which(w >= 1.0 - 0.5 * alpha))
      return(c(lower = x[LB],
               upper = x[UB]))
    }
    CI_bounds = 
      apply(object$proposal_draws,2,
            CI_from_weighted_sample,
            w = object$importance_sampling_weights)
    summ$Lower = 
      CI_bounds["lower",]
    summ$Upper = 
      CI_bounds["upper", ]
  }
  
  # Exponentiate
  if(( (object$family$family == "binomial") & 
       (object$family$link != "logit") ) | 
     ( (object$family$family == "poisson") & 
       (object$family$link != "log") )){
    interpretable_scale = FALSE
  }
  
  if(interpretable_scale){
    paste0("\n----------\n\nValues given in terms of ",
           ifelse(object$family$family == "binomial",
                  "odds ratios",
                  "rate ratios")
    ) %>% 
    cat()
    cat("\n\n----------\n\n")
    summ = summ[-1,]
    summ[,c("Post Mean","Lower","Upper")] %<>%
      exp()
    summ[,"ROPE bounds"] = 
      paste("(",
            round(exp(-object$ROPE[-1]),3),
            ",",
            round(exp(object$ROPE[-1]),3),
            ")",
            sep="")
  }
  
  summ
}
