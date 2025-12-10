#' Mediation using Bayesian methods
#' 
#' 
#' Mediation analysis done in the framework of Imai et al. (2010).  Currently 
#' only applicable to linear models.
#' 
#' @details
#' The model is the same as that of Imai et al. (2010):
#' \deqn{
#'  M_i(X) =  w_i'\alpha_m + X\beta_m + \epsilon_{m,i}, \\
#'  y_i(X, M(\tilde X)) = w_i'\alpha_y + X\beta_y + M(\tilde X)\gamma + \epsilon_{y,i}, \\
#'  \epsilon_{m,i} \overset{iid}{\sim} N(0,\sigma^2_m), \\
#'  \epsilon_{y,i} \overset{iid}{\sim} N(0,\sigma^2_y), \\
#' }
#' where \eqn{M_i(X)} is the mediator as a function of the treatment variable
#' \eqn{X}, and \eqn{w_i} are confounder covariates.
#' 
#' 
#' Unlike the \code{mediation} R package, the estimation in \code{mediate_b} 
#' is fully Bayesian (as opposed to "quasi-Bayesian").
#' 
#' 
#' @references 
#' 
#' Imai, Kosuke, et al. 
#' “A General Approach to Causal Mediation Analysis.” Psychological Methods, 
#' vol. 15, no. 4, 2010, pp. 309–34, https://doi.org/10.1037/a0020761.
#' 
#' @param model_m a fitted model object of class lm_b for mediator.
#' @param model_y a fitted model object of class lm_b for outcome.
#' @param sims number of posterior draws for causal quantities
#' @param treat a character string indicating the name of the 
#' treatment variable used in the models.  NOTE: Treatment variable must be
#' numeric.
#' @param control_value value of the treatment variable used as the 
#' control condition. Default is the 1st quintile of the treat variable.
#' @param treat_value value of the treatment variable used as the treatment condition. 
#' Default is the 4th quintile of the treat variable.
#' @param CI_level numeric. Credible interval level.
#' @param simple logical.  Set to TRUE only if there are no interactions.
#' @param show_progress logical.
#' 
#' @return A list with two elements:
#' \itemize{
#'  \item summary - Posterior mean and credible intervals for causal quantities
#'  \item counterfactual_draws (if \code{simple = FALSE}) - list with the following elements:
#'  \itemize{
#'    \item E_y_00 - The expected value of the counterfactual y(control_value,M(control_value))
#'    \item E_y_01 - The expected value of the counterfactual y(control_value,M(treat_value))
#'    \item E_y_10 - The expected value of the counterfactual y(treat_value,M(control_value))
#'    \item E_y_11 - The expected value of the counterfactual y(treat_value,M(treat_value)) 
#'  }
#'  \item control_value, treat_value
#' }
#' 
#' @export

mediate_b = function(model_m,
                     model_y,
                     sims = 1e3,
                     treat,
                     control_value,
                     treat_value,
                     CI_level = 0.95,
                     simple = TRUE,
                     show_progress = TRUE){
  
  a = 1 - CI_level
  
  if(!all.equal(model_m$data,model_y$data)){
    stop("Data in model_m and model_y must match.")
  }
  
  mediator = as.character(model_m$formula)[[2]]
  
  if(missing(control_value)){
    warning(paste0("control_value missing; set to be the 1st quintile of ",
                   treat))
    control_value = quantile(model_m$data[[treat]],probs = 0.2)
  }
  if(missing(treat_value)){
    warning(paste0("treat_value missing; set to be the 4th quintile of ",
                   treat))
    treat_value = quantile(model_m$data[[treat]],probs = 0.8)
  }
  
  
  if(simple){
    
    # Get posterior draws for ACME and ADE
    mediator_draws = 
      get_posterior_draws(model_m,
                          n_draws = sims)
    outcome_draws = 
      get_posterior_draws(model_y,
                          n_draws = sims)
    
    acme_draws = 
      (treat_value - control_value) * 
      mediator_draws[,treat] *
      outcome_draws[,all.vars(model_m$formula)[1]]
    
    ade_draws = 
      (treat_value - control_value) * 
      outcome_draws[,treat]
    
    ate_draws = 
      acme_draws + ade_draws
    
    ret = list()
    ret$summary = 
      data.frame(Estimand = c("ACME",
                              "ADE",
                              "Total Effect",
                              "Prop. Mediated"),
                 Estimate = c(mean(acme_draws),
                              mean(ade_draws),
                              mean(ate_draws),
                              mean(acme_draws / ate_draws)),
                 Lower = c(quantile(acme_draws,probs = 0.5 * a),
                           quantile(ade_draws,probs = 0.5 * a),
                           quantile(ate_draws,probs = 0.5 * a),
                           quantile(acme_draws / ate_draws,probs = 0.5 * a)),
                 Upper = c(quantile(acme_draws,probs = 1.0 - 0.5 * a),
                           quantile(ade_draws,probs = 1.0 - 0.5 * a),
                           quantile(ate_draws,probs = 1.0 - 0.5 * a),
                           quantile(acme_draws / ate_draws,probs = 1.0 - 0.5 * a))
      )
    
  }else{
    
    # Get mediator counterfactual posterior draws
    counterfactual_data0 = 
      counterfactual_data1 = 
      model_m$data
    counterfactual_data0[[treat]] = control_value
    counterfactual_data1[[treat]] = treat_value
    
    M_0 =
      predict(model_m,
              newdata = counterfactual_data0,
              n_draws = sims)
    
    M_1 =
      predict(model_m,
              newdata = counterfactual_data1,
              n_draws = sims)
    
    
    # Get outcome counterfactual posterior draws
    y_00 = y_01 = y_10 = y_11 = model_y$data
    if(show_progress) pb = txtProgressBar(0,sims,style=3)
    for(iter in 1:sims){
      
      counterfactual_data0[[mediator]] = M_0[[paste("y_new",iter,sep="")]]
      y_00[[paste("y_new",iter,sep="")]] = 
        predict(model_y,
                newdata = counterfactual_data0,
                n_draws = 1)$y_new1
      
      counterfactual_data0[[mediator]] = M_1[[paste("y_new",iter,sep="")]]
      y_01[[paste("y_new",iter,sep="")]] = 
        predict(model_y,
                newdata = counterfactual_data0,
                n_draws = 1)$y_new1
      
      counterfactual_data1[[mediator]] = M_0[[paste("y_new",iter,sep="")]]
      y_10[[paste("y_new",iter,sep="")]] = 
        predict(model_y,
                newdata = counterfactual_data1,
                n_draws = 1)$y_new1
      
      counterfactual_data1[[mediator]] = M_1[[paste("y_new",iter,sep="")]]
      y_11[[paste("y_new",iter,sep="")]] = 
        predict(model_y,
                newdata = counterfactual_data1,
                n_draws = 1)$y_new1
        
      if(show_progress) setTxtProgressBar(pb,iter)
    }
    
    # Now get means over the empirical distribution of X
    E_y_00 = 
      y_00[,ncol(y_00) + 1 - sims:1] |> 
      as.matrix() |> 
      colMeans()
    E_y_01 = 
      y_01[,ncol(y_00) + 1 - sims:1] |> 
      as.matrix() |> 
      colMeans()
    E_y_10 = 
      y_10[,ncol(y_00) + 1 - sims:1] |> 
      as.matrix() |> 
      colMeans()
    E_y_11 = 
      y_11[,ncol(y_00) + 1 - sims:1] |> 
      as.matrix() |> 
      colMeans()
      
    # Get causal quantities
    tot_eff = mean(E_y_11 - E_y_00)
    
    ACME = c(control = mean(E_y_01 -E_y_00),
             treat = mean(E_y_11 -E_y_10))
    
    ADE = c(control = mean(E_y_10 -E_y_00),
            treat = mean(E_y_11 -E_y_01))
    
    
    # Put it together to return
    ret = list()
    ret$summary = 
      tibble(Estimand = c("ACME (Control)",
                              "ACME (Treat)",
                              "ADE (Control)",
                              "ADE (Treat)",
                              "Total Effect",
                              "Prop. Mediated"),
                 Estimate = c(ACME["control"],
                              ACME["treat"],
                              ADE["control"],
                              ADE["treat"],
                              tot_eff,
                              mean(0.5 * (E_y_01 - E_y_00 + E_y_11 -E_y_10) / 
                                     (E_y_11 - E_y_00))),
                 Lower = c(quantile(E_y_01 -E_y_00,probs = 0.5 * a),
                           quantile(E_y_11 -E_y_10,probs = 0.5 * a),
                           quantile(E_y_10 -E_y_00,probs = 0.5 * a),
                           quantile(E_y_11 -E_y_01,probs = 0.5 * a),
                           quantile(E_y_11 - E_y_00,probs = 0.5 * a),
                           quantile(0.5 * (E_y_01 - E_y_00 + E_y_11 -E_y_10) / 
                                      (E_y_11 - E_y_00),
                                    probs = 0.5 * a)),
                 Upper = c(quantile(E_y_01 -E_y_00,probs = 1.0 - 0.5 * a),
                           quantile(E_y_11 -E_y_10,probs = 1.0 - 0.5 * a),
                           quantile(E_y_10 -E_y_00,probs = 1.0 - 0.5 * a),
                           quantile(E_y_11 -E_y_01,probs = 1.0 - 0.5 * a),
                           quantile(E_y_11 - E_y_00,probs = 1.0 - 0.5 * a),
                           quantile(0.5 * (E_y_01 - E_y_00 + E_y_11 -E_y_10) / 
                                      (E_y_11 - E_y_00),
                                    probs = 1.0 - 0.5 * a))
      )
    
    ret$counterfactual_draws = 
      list(E_y_00 = E_y_00,
           E_y_01 = E_y_01,
           E_y_10 = E_y_10,
           E_y_11 = E_y_11
           )
  }
  
  
  ret$treat_value = treat_value
  ret$control_value = control_value
  
  return(ret)
}


