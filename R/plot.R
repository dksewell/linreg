#' @name plot
#' 
#' @title Print linreg objects.
#' 
#' @param x A linreg object
#' @param type character. Select any of "diagnostics" ("dx" is also allowed),
#'  "pdp" (partial dependence plot), "ci band", and/or "pi band".  NOTE: the CI 
#'  and PI bands only work for numeric variables.
#' @param variable character. If type = "pdp" , which variable should be plotted?
#' @param exemplar_covariates data.frame or tibble with exactly one row.  
#' Used to fix other covariates while varying the variable of interest for the plot.
#' @param combine_pi_ci logical. If type includes both "ci band" and "pi band", 
#' should the credible band be superimposed on the prediction band or 
#' plotted separately?
#' @param variable_seq_length integer. Number of points used to draw pdp.
#' @param return_as_list logical.  If TRUE, a list of ggplots will be returned, 
#' rather than a single plot produced by the patchwork package.
#' 
#' @import ggplot2
#' @import patchwork
#' @importFrom cluster pam
#' 

#' @rdname plot
#' @export
print.lm_b = function(x,
                      type = c("diagnostics",
                               "pdp",
                               "ci band",
                               "pi band")[1],
                      variable,
                      exemplar_covariates,
                      combine_pi_ci = TRUE,
                      variable_seq_length = 100,
                      return_as_list = FALSE){
  
  type = c("diagnostics",
           "diagnostics",
           "pdp",
           "ci band",
           "pi band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "pdp",
                               "ci band",
                               "pi band"))]
  
  if(missing(variable)) variable = all.vars(x$formula)[-1]
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    dx_data = 
      tibble(yhat = x$fitted,
             epsilon = x$residuals)
    
    plot_list[["fitted_vs_residuals"]] =
      dx_data |>
      ggplot(aes(y = epsilon,x = yhat)) +
      geom_hline(yintercept = 0,
                 linetype = 2,
                 color = "gray35") +
      geom_point(alpha = 0.6) +
      xlab(expression(hat(y))) +
      ylab(expression(hat(epsilon))) +
      theme_classic() +
      ggtitle("Fitted vs. Residuals")
      
    plot_list[["qqnorm"]] =
      dx_data |>
      ggplot(aes(sample = yhat)) +
      geom_qq(alpha = 0.3) + 
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
    
  }# End: diagnostics
  
  
  # Get unique values and x sequences for plots
  if( length(intersect(c("pdp","ci band","pi band"),
                type)) > 0){
    
    x_unique = 
      lapply(variable,
             function(v) unique(x$data[[v]]))
    x_seq = 
      lapply(x_unique,
             function(xvals){
               if(length(xvals) > variable_seq_length){
                 return( 
                   seq(min(xvals),
                       max(xvals),
                       l = variable_seq_length)
                 )
               }else{
                 if(is.numeric(xvals)){
                   return(sort(xvals))
                 }else{
                   if(is.character(xvals)){
                     return(
                       factor(sort(xvals),
                              levels = sort(xvals))
                     )
                   }else{
                     return(xvals)
                   }
                 }
               }
             })
    
    names(x_unique) = 
      names(x_seq) = variable
  }# End: Get unique values and x_seq
  
  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble(var_of_interest = x_seq[[v]],
               y = 0.0)
      for(i in 1:length(x_seq[[v]])){
        temp_preds = 
          predict(x,
                  newdata = 
                    x$data |>
                    dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]))
        newdata$y[i] = mean(temp_preds$Estimate)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        x$data |>
        ggplot(aes(x = !!sym(variable[v]),
                   y = !!sym(all.vars(x$formula)[1]))) + 
        geom_point(alpha = 0.2)
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_line(data = newdata,
                    aes(x = var_of_interest,
                        y = y))
      }else{
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                    aes(x = var_of_interest,
                        y = y),
                    size = 5)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        plot_list[[paste0("pdp_",variable[v])]] + 
        xlab(variable[v]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }# End: PDP
  
  
  # If drawing CI/PI bands, get reference covariate values and prediction/CIs
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) %>% 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
    
    # Get CI and PI values
    newdata = list()
    for(v in variable){
      newdata[[v]] = 
        tibble(!!v := x_seq[[v]])
      for(j in setdiff(names(exemplar_covariates),v)){
        if(is.character(exemplar_covariates[[j]])){
          newdata[[v]][[j]] = 
            factor(exemplar_covariates[[j]],
                   levels = unique(x$data[[j]]))
        }else{
          newdata[[v]][[j]] = exemplar_covariates[[j]]
        }
      }
      
      newdata[[v]] = 
        predict(x,
                newdata = newdata[[v]])
    }
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pi band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    for(v in variable){
      if(is.numeric(x$data[[v]])){
        plot_list[[paste0("band_",v)]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[paste0("band_",v)]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_violin(alpha = 0.2)
      }
    }
    
    for(v in variable){
      
      if(class(x) == "lm_b"){
        
        if(is.numeric(x_seq[[v]])){
          plot_list[[paste0("band_",v)]] =
            plot_list[[paste0("band_",v)]] +
            geom_ribbon(data = newdata[[v]],
                        aes(ymin = PI_lower,
                            ymax = PI_upper),
                        fill = "lightsteelblue3",
                        alpha = 0.5) +
            geom_line(data = newdata[[v]],
                      aes(x = !!sym(v),
                          y = Estimate))
        }else{
          plot_list[[paste0("band_",v)]] =
            plot_list[[paste0("band_",v)]] +
            geom_errorbar(data = newdata[[v]],
                          aes(x = !!sym(v),
                              ymin = PI_lower,
                              ymax = PI_upper),
                          color = "lightsteelblue3") +
            geom_point(data = newdata[[v]],
                       aes(x = !!sym(v),
                           y = Estimate),
                       size = 5)
        }
      }
      
      if(class(x) == "aov_b"){
        
      }
      
    }
      
      
    if(class(x) == "np_lm_b") stop("Prediction interval band not applicable for np_lm_b objects.")
      
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if(!combine_pi_ci){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_violin(alpha = 0.2)
        }
      }
    }
    
    for(v in variable){
      
      if(class(x) == "lm_b"){
        
        if(is.numeric(x_seq[[v]])){
          plot_list[[paste0("band_",v)]] =
            plot_list[[paste0("band_",v)]] +
            geom_ribbon(data = newdata[[v]],
                        aes(ymin = CI_lower,
                            ymax = CI_upper),
                        fill = "steelblue4",
                        alpha = 0.5) +
            geom_line(data = newdata[[v]],
                      aes(x = !!sym(v),
                          y = Estimate))
        }else{
          plot_list[[paste0("band_",v)]] =
            plot_list[[paste0("band_",v)]] +
            geom_errorbar(data = newdata[[v]],
                          aes(x = !!sym(v),
                              ymin = CI_lower,
                              ymax = CI_upper),
                          color = "steelblue4") +
            geom_point(data = newdata[[v]],
                      aes(x = !!sym(v),
                          y = Estimate),
                      size = 5)
        }
      }
      
      if(class(x) == "aov_b"){
        
      }
      
      if(class(x) == "np_lm_b"){
        
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    for(v in variable){
      plot_list[[paste0("band_",v)]] =
        plot_list[[paste0("band_",v)]] +
        theme_classic() +
        ggtitle(paste0(
          ifelse( ("pi band" %in% type) & ("ci band" %in% type),
                  paste0("CI and PI bands for ",v),
                  ifelse("pi band" %in% type,
                         paste0("PI band for ",v),
                         paste0("CI band for ",v)))
        ))
        
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
      # eval(parse(text = paste(
      #   paste("plot_list[['",
      #         names(plot_list),
      #         "']]",sep=''),
      # collapse = "+")))
      )
  }
  
  
}

