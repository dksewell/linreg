#' Create survival curves
#' 
#' Use the semi-parametric piecewise exponential survival model 
#' to fit a survival curve
#' 
#' 
#' @references 
#' Qing Y, Thall PF, Yuan Y. A Bayesian piecewise exponential phase II design for monitoring a time-to-event endpoint. Pharm Stat. 2023 Jan;22(1):34-44. doi: 10.1002/pst.2256. Epub 2022
#' 
#' 
#' @importFrom survival Surv
#' @export




survfit_b = function(Surv_object,
                     prior_shape,
                     prior_rate,
                     max_n_time_bins,
                     n_time_bins){
  
  
  if(!inherits(Surv_object, "Surv")){
    stop("Surv_object must be a survival::Surv object.")
  }
  
  # Extract
  time = Surv_object[,1]
  status = Surv_object[,2]
  max_time = max(time)
  
  
  # If missing prior hyperparameters, use data to get reasonable guess based on t_i\sim exp(\lambda)
  if(missing(prior_shape) | missing(prior_rate)){
    # Based on posterior updating with unit information gain
    prior_shape = 
      0.001 + mean(status)
    prior_rate = 
      0.001 + mean(time)
  }
  
  
  # Create helper for computing posterior parameters
  get_post_parms = function(breakpoints){
    bpt_assignment = 
      cut(time,
          breakpoints) |> 
      as.integer()
    n_j = 
      bpt_assignment[which(dplyr::near(status,1))] |> 
      table() |> 
      as.vector()
    bpt_diffs = diff(breakpoints)
    
    e_ij = matrix(0.0,length(time),length(breakpoints) - 1)
    for(j in 1:ncol(e_ij)){
      e_ij[which(bpt_assignment > j),j] = 
        bpt_diffs[j]
      
      happened_in_this_bin = 
        which(bpt_assignment == j)
      e_ij[happened_in_this_bin,j] = 
        time[happened_in_this_bin] -
        breakpoints[j]
    }
    
    cbind(a_j_tilde =
            prior_shape + n_j,
          b_j_tilde = 
            prior_rate + colSums(e_ij)
    ) 
  }
  
  # Try out multiple breakpoints, select the optimal via marginal likelihood
  if(missing(n_time_bins)){
    trimmed_time = c(time[-which(time == max_time)],max_time)
    
    if(missing(max_n_time_bins))
      max_n_time_bins = max(2,(NROW(time)/5))
    
    ml_values = c(NA,numeric(max_n_time_bins - 1))
    for(J in 2:max_n_time_bins){
      bpts = quantile(trimmed_time,seq(0,1,l = 1 + J))
      bpts[1] = 0.0
      a_b_tilde = 
        get_post_parms(bpts)
      
      ml_values[J] = 
        sum(lgamma(a_b_tilde[,1])) -
        length(time) * lgamma(prior_shape) +
        length(time) * prior_shape * log(prior_rate) - 
        sum(a_b_tilde[,1] * log(a_b_tilde[,2]))
    }
    
    J_opt = which.max(ml_values)
  }else{
    J_opt = n_time_bins
  }
  
  # Get best (or pre-specified) breakpoints
  ret = list()
  bpts = quantile(trimmed_time,seq(0,1,l = 1 + J_opt))
  bpts[1] = 0.0
  ret$posterior_parameters = 
    get_post_parms(bpts)
  
  # Get intervals
  if(J_opt == 2){
    ret$intervals = 
      rep(bpts, c(1,2,1)) |> 
      matrix(ncol = 2,
             byrow = TRUE)
  }else{
    ret$intervals = 
      rep(bpts, 
          c(1,rep(2,length(bpts) - 2),1)) |> 
      matrix(ncol = 2,
             byrow = TRUE)
  }
  
  # Return data also
  ret$data = Surv_object
  
  return(structure(ret,
                   class = "survfit_b"))
}
