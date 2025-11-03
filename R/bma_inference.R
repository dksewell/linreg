#' Estimates and CI from BMA
#' 
#' @param formula A formula specifying the model.
#' @param data Data used in linear regression model
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  
#' @param CI_level Level for credible interval
#' @param mc_draws Integer.  Number of draws in Monte Carlo integration.
#' @param future.seed Integer. (Optional) If parallelizing through future, 
#' make sure to set your seed for reproducibility!
#' @param ... Other arguments for BMS::bms().
#' 
#' @return A list with the following elements:
#' \itemize{
#'  \item summary Tibble with point and interval estimates
#'  \item lm_b_fits A list of lm_b fits using zellner's g prior for
#'  all the top models from BMS::bms()
#'  \item hyperparameters A named list with the user-specified zellner's g value.
#'  \item posterior_draws matrix of posterior draws of the regression parameters, 
#'  marginalizing out the model
#' }
#' 
#' @importFrom BMS bms
#' @export


bma_inference = function(formula,
                         data,
                         zellner_g = nrow(data),
                         CI_level = 0.95,
                         mc_draws = 5e4,
                         future.seed = 1,
                         ...){
  
  # Use BMS package to get top (a posteriori) models
  X.data = 
    cbind(
      y = data[[all.vars(formula)[1]]],
      model.matrix(formula,data)[,-1]
    )
  
  bms_fit = 
    BMS::bms(X.data,
             g = zellner_g,
             iter = mc_draws)
  
  
  # Get the posterior probabilities of the models
  model_post_probs = 
    bms_fit$topmod$ncount()
  model_post_probs =
    model_post_probs / sum(model_post_probs)
  
  # Get number of posterior samples per model
  mc_draws_by_model = 
    round(mc_draws * model_post_probs)
  
  # Extract a matrix showing us which variables are used 
  # in each of our top-most models
  var_inclusion = 
    bms_fit$topmod$bool_binary()
  if(any(mc_draws_by_model == 0)) var_inclusion = var_inclusion[,-which(mc_draws_by_model == 0)]
  
  # Fit top models
  X.data = as.data.frame(X.data)
  full_fits = 
    future_lapply(1:ncol(var_inclusion),
                  function(i){
                    lm_b(paste0(all.vars(formula)[1], " ~ ", 
                                paste(colnames(X.data)[-1][as.logical(var_inclusion[,i])],
                                      collapse = " + ")) %>% 
                           as.formula(),
                         data = X.data,
                         prior = "zellner",
                         zellner_g = zellner_g)
                  })
  
  # Get posterior samples
  post_samples = 
    future_lapply(1:length(full_fits),
                  function(i){
                    samples = 
                      get_posterior_draws(full_fits[[i]],
                                          n_draws = mc_draws_by_model[i]) %>% 
                      as_tibble()
                    if(ncol(samples) < ncol(X.data)){ # Yes X.data includes y, but the samples also ought to include s2
                      for(j in setdiff(colnames(X.data)[-1],
                                       colnames(samples))) samples[[j]] = 0.0
                    }
                    return(samples)
                  },
                  future.seed = future.seed)
  post_samples = 
    do.call(bind_rows,post_samples)
  
  post_samples %<>% 
    relocate(all_of(c(colnames(X.data)[-1],"s2")),
             .after = "(Intercept)")
  
  # Summarize results
  alpha = 1.0 - CI_level
  ## Get values for ROPE
  s_y = sd(X.data[,1])
  s_j = apply(X.data[,-1,drop = FALSE],2,sd)
  boundaries = 
    matrix(0.05 * s_y / c(NA,s_j,NA),
           nr = nrow(post_samples),
           nc = ncol(post_samples),
           byrow = TRUE)
  
  ## Compile results
  results = 
    tibble(Variable = colnames(post_samples),
           `Post Mean` = colMeans(post_samples),
           Lower = 
             apply(post_samples,2,quantile,probs = alpha/2),
           Upper =
             apply(post_samples,2,quantile,probs = 1.0 - alpha/2),
           ROPE = 
             colMeans(
               (-boundaries < post_samples) & 
                 (boundaries > post_samples)
             ),
           `Prob Dir` = 
             c(sapply(apply(post_samples[,-ncol(post_samples)],
                            2,
                            function(x) mean(x < 0)),
                    function(x) max(x,1-x)),
               NA)
    )
  
  # Get fitted values and residuals
  fitted = 
    drop(
      cbind(1.0,as.matrix(X.data[,-1])) %*% results$`Post Mean`[-nrow(results)]
    )
  residuals = 
    drop(
      X.data[,1] - fitted
    )
  
  return_object = 
    list(summary = results,
         lm_b_fits = full_fits,
         bms_fit = bms_fit,
         hyperparms = list(zellner_g = zellner_g),
         posterior_draws = post_samples,
         fitted = fitted,
         residuals = residuals,
         formula = formula,
         data = data)
  
  class(return_object) = "bma"
}
