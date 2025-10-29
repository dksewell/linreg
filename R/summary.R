#' @name summary
#' 
#' @title Summary functions for lm_b and aov_b objects
#' 
#' @param object lm_b or aov_b object
#' 
#' @export

#' @rdname summary
#' @export
summary.lm_b = function(object,
                        CI_level = 0.95){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  if(object$prior == "improper"){
    summ$Lower = 
      qlst(alpha/2,
           object$post_parms$a_tilde,
           object$post_parms$mu_tilde,
           sqrt(object$post_parms$b_tilde/object$post_parms$a_tilde * 
                  diag(qr.solve(object$post_parms$V_tilde))))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           object$post_parms$a_tilde,
           object$post_parms$mu_tilde,
           sqrt(object$post_parms$b_tilde/object$post_parms$a_tilde * 
                  diag(qr.solve(object$post_parms$V_tilde))))
  }else{
    summ$Lower = 
      qlst(alpha/2,
           nrow(object$data) - length(object$post_parms$mu_tilde),
           object$post_parms$mu_tilde,
           sqrt(diag(object$post_parms$Sigma)))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           nrow(object$data) - length(object$post_parms$mu_tilde),
           object$post_parms$mu_tilde,
           sqrt(diag(object$post_parms$Sigma)))
  }
  
  summ
}

#' @rdname summary
#' @export
summary.aov_b = function(object,
                         CI_level = 0.95){
  alpha = alpha = 1 - CI_level
  summ = object$summary
  pw_summ = object$pairwise_summary
  
  
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
  print(pw_summ)
  
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
summary.np_lm_b = function(object){
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
