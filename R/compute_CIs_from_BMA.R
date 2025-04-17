#' Estimates and CI from BMA
#' 
#' @param object An object returned from BMS::bms().  Make sure you return 
#' a reasonably large number of models (nmodel)
#' @param data Data used in linear regression model
#' @param y_name character for response variable.
#' @param covar_names character vector for covariates in fullest model
#' @param CI_level Level for credible interval
#' @param mc_draws Integer.  Number of draws in Monte Carlo integration.
#' 
#' @return A list with the following elements:
#' \itemize{
#'  \item results Dataframe with point and interval estimates
#'  \item full_fits a list of lm_b fits using zellner's g prior for
#'  all the top models from BMS::bms()
#'  \item mc_draws matrix of posterior draws of the regression coefficients, 
#'  marginalizing out the model
#' }
#' 
#' @importFrom mvtnorm rmvt
#' @export


get_est_and_ci = function(object,
                          data,
                          y_name,
                          covar_names,
                          CI_level = 0.95,
                          mc_draws = 5e4){
  
  # Extract a matrix showing us which variables are used 
  # in each of our top-most models
  var_inclusion = 
    bma_fits$uniform_on_models$topmod$bool_binary()
  dim(var_inclusion)
  
  full_fits = 
    lapply(1:ncol(var_inclusion),
           function(i){
             lm_b(paste0(y_name,
                         " ~ ",
                         paste(covar_names[as.logical(var_inclusion[,i])],
                               collapse = " + ")) %>% 
                    as.formula(),
                  data = data,
                  prior = "zellner")
           })
  
  # Get the posterior probabilities of the models
  model_post_probs = 
    object$topmod$ncount()
  model_post_probs =
    model_post_probs / sum(model_post_probs)
  
  full_fits_summary = 
    do.call(bind_rows,
            lapply(1:ncol(var_inclusion),
                   function(i){
                     full_fits[[i]]$summary %>% 
                       select(-Prob.Dir) %>% 
                       bind_rows(tibble(Variable = covar_names[!as.logical(var_inclusion[,i])],
                                        Post.Mean = 0,
                                        Lower = 0,
                                        Upper = 0,
                                        ROPE = 1.0)) %>% 
                       mutate(model_prob = model_post_probs[i])
                   })) %>% 
    select(Variable,Post.Mean,ROPE,model_prob) %>% 
    mutate(across(Post.Mean:ROPE,
                  function(x) x * model_prob)) %>% 
    select(-model_prob) %>% 
    group_by(Variable) %>% 
    summarize(across(Post.Mean:ROPE,sum)) %>% 
    rename(`Post Mean` = Post.Mean)
  
  # Get rmixt for MC integration
  model_draw = sample(length(model_post_probs),
                      mc_draws,
                      replace = TRUE,
                      prob = model_post_probs)
  ctable <- tabulate(model_draw, nbins = length(model_post_probs))
  t_draws = matrix(0.0,mc_draws,length(covar_names) + 1,
                   dimnames = list(NULL,c("(Intercept)",covar_names)))
  for(m in 1:length(model_post_probs)){
    if(ctable[m] > 0){
      t_draws[which(model_draw == m),c(T,as.logical(var_inclusion[,m]))] = 
        rmvt(ctable[m],
             df = full_fits[[m]]$post_parms$a_tilde,
             delta = full_fits[[m]]$post_parms$mu_tilde[,1],
             sigma = 
               full_fits[[m]]$post_parms$b_tilde / 
               full_fits[[m]]$post_parms$a_tilde * 
               chol2inv(chol(full_fits[[m]]$post_parms$V_tilde)))
    }
  }
  
  alpha = 1.0 - CI_level
  results = 
    full_fits_summary %>% 
    left_join(
      tibble(Variable = colnames(t_draws),
             `Post Median` = apply(t_draws,2,quantile,probs = 0.5),
             Lower = apply(t_draws,2,quantile,probs = alpha/2.0),
             Upper = apply(t_draws,2,quantile,probs = 1.0 - alpha/2.0),
             `Pr(<0)` = apply(t_draws,2, function(x) mean(x < 0)),
             `Pr(=0)` = apply(t_draws,2, function(x) mean(x == 0)),
             `Pr(>0)` = apply(t_draws,2, function(x) mean(x > 0))),
      by = "Variable") %>% 
    relocate(ROPE, .after = `Pr(>0)`)
      
  
  return(list(results = results,
              full_fits = full_fits,
              mc_draws = t_draws))
}
