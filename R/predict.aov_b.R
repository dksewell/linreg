#' Predict method for aov_b model fits
#' 
#' 
#' @param object Object of class aov_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. If set to NULL, prediction intervals will only be given on a group
#' level basis.
#' @param CI_level numeric. Credible interval level.
#' 
#' @return list with the following elements.
#' \itemize{
#'  \item newdata - list with posterior_draws and summary
#'  \item by_group - list with posterior_draws and summary
#' }
#' 
#' @export predict.aov_b
#' @export



predict.aov_b = function(object, newdata, CI_level = 0.95){
  
  ret = list()
  n_joint_draws = nrow(object$posterior_draws)
  a = 1 - CI_level
  G = length(levels(object$data$group))
  
  if(missing(newdata)){
    newdata = object$data
  }
  
  # Get prediction estimates and intervals for new data
  if(!is.null(newdata)){
    ## Get design matrix for new predictions
    X = model.matrix(~ -1 + group,data = newdata)
    
    ## Get predictions with noise
    ret$newdata = list()
    if(ncol(object$posterior_draws) == G+1){
      ### Get mean
      ret$newdata$posterior_draws = 
        tcrossprod(X,
                   object$posterior_draws[,1:G])
      ### Add noise
      for(i in 1:nrow(X)){
        ret$newdata$posterior_draws[i,] = 
          ret$newdata$posterior_draws[i,] + 
          rnorm(n_joint_draws,
                sd = sqrt(object$posterior_draws[,G + 1]))
      }
    }else{
      ### Get mean
      ret$newdata$posterior_draws = 
        tcrossprod(X,
                   object$posterior_draws[,1:G])
      ### Add noise
      for(i in 1:nrow(X)){
        ret$newdata$posterior_draws[i,] = 
          ret$newdata$posterior_draws[i,] + 
          rnorm(n_joint_draws,
                sd = sqrt(object$posterior_draws[,G + as.integer(newdata$group)]))
          
      }
    }
    ret$newdata$summary = 
      cbind(Estimate = rowMeans(ret$newdata$posterior_draws),
            Lower = apply(ret$newdata$posterior_draws,1,quantile,probs = a/2),
            Upper = apply(ret$newdata$posterior_draws,1,quantile,probs = 1.0 - a/2))
  }
  
  
  # Get prediction intervals for each group
  ret$by_group = list()
  ret$by_group$posterior_draws = 
    object$posterior_draws[,1:G]
  
  for(g in 1:G){
    if(ncol(object$posterior_draws) == G+1){
      ret$by_group$posterior_draws[,g] = 
        ret$by_group$posterior_draws[,g] +
        rnorm(n_joint_draws,
              sd = sqrt(object$posterior_draws[,G + 1]))
    }else{
      ret$by_group$posterior_draws[,g] = 
        ret$by_group$posterior_draws[,g] +
        rnorm(n_joint_draws,
              sd = sqrt(object$posterior_draws[,G + g]))
    }
  }
  ret$by_group$summary = 
    data.frame(Variable = object$summary$Variable[1:G],
               Estimate = object$summary$Estimate[1:G],
               Lower = apply(ret$by_group$posterior_draws,2,quantile,probs = a/2),
               Upper = apply(ret$by_group$posterior_draws,2,quantile,probs = 1.0 - a/2))
  rownames(ret$by_group$summary) = NULL
  
  return(ret)
}




