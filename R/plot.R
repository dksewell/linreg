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
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' 
#' @import ggplot2
#' @import patchwork
#' @importFrom cluster pam
#' 

#' @rdname plot
#' @export
plot.lm_b = function(x,
                      type = c("diagnostics",
                               "pdp",
                               "ci band",
                               "pi band"),
                      variable,
                      exemplar_covariates,
                      combine_pi_ci = TRUE,
                      variable_seq_length = 100,
                      return_as_list = FALSE,
                     CI_level = 0.95,
                     PI_level = 0.95){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
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
      ggplot(aes(sample = epsilon)) +
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
                    dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]),
                  CI_level = CI_level,
                  PI_level = PI_level)
        newdata$y[i] = mean(temp_preds$`Post Mean`)
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
                newdata = newdata[[v]],
                CI_level = CI_level,
                PI_level = PI_level)
    }
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pi band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("ci band" %in% type),
                      "pi_band_","band_"),v)
      
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_violin(alpha = 0.2)
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("ci band" %in% type),
                      "pi_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = PI_lower,
                          ymax = PI_upper),
                      fill = "lightsteelblue3",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = !!sym(v),
                        y = `Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]],
                        aes(x = !!sym(v),
                            ymin = PI_lower,
                            ymax = PI_upper),
                        color = "lightsteelblue3") +
          geom_point(data = newdata[[v]],
                     aes(x = !!sym(v),
                         y = `Post Mean`),
                     size = 5)
      }
    
      
    }
      
      
      
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("ci_band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("ci_band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_violin(alpha = 0.2)
        }
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("pi band" %in% type),
                      "ci_band_","band_"),v)
      
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = CI_lower,
                          ymax = CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = !!sym(v),
                        y = `Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]],
                        aes(x = !!sym(v),
                            ymin = CI_lower,
                            ymax = CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                    aes(x = !!sym(v),
                        y = `Post Mean`),
                    size = 5)
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    for(v in variable){
      
      for(j in names(plot_list)[grepl("band",names(plot_list)) & grepl(v,names(plot_list))]){
        plot_list[[j]] =
          plot_list[[j]] +
          theme_classic() +
          ggtitle(
            paste0(
              ifelse(
                grepl("pi_",j),
                paste0("PI band for ",v),
                ifelse(grepl("ci_",j),
                       paste0("CI band for ",v),
                       paste0("CI and PI bands for ",v)
                )
              )
            )
          )
      }
        
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
      )
  }
  
  
}



#' @rdname plot
#' @export
plot.aov_b = function(x,
                      type = c("diagnostics",
                               "ci band",
                               "pi band"),
                      combine_pi_ci = TRUE,
                      return_as_list = FALSE,
                      CI_level = 0.95,
                      PI_level = 0.95){
  
  type = c("diagnostics",
           "diagnostics",
           "ci band",
           "pi band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "ci band",
                               "pi band"))]
  
  
  plot_list = list()
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    dx_data =
      tibble(group = x$data$group,
             yhat = x$fitted,
             epsilon = x$residuals)
    
    plot_list[["residuals_by_group"]] =
      dx_data |>
      ggplot(aes(y = epsilon,x = group)) +
      geom_hline(yintercept = 0,
                 linetype = 2,
                 color = "gray35") +
      geom_violin(alpha = 0.6) +
      xlab(all.vars(x$formula)[2]) +
      ylab(expression(hat(epsilon))) +
      theme_classic() +
      ggtitle("Residual plot by group")
    
    plot_list[["qqnorm"]] =
      dx_data |>
      ggplot(aes(sample = epsilon)) +
      geom_qq(alpha = 0.3) +
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
    
  }# End: diagnostics
  
  
  # If drawing CI/PI bands, get newdata for prediction/CIs
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    
    # Get CI and PI values
    newdata =
      predict(x,
              CI_level = CI_level,
              PI_level = PI_level)
    
  }# End: Get newdata
  
  
  # Prediction Band plots
  if("pi band" %in% type){
    
    # Get starter plots
    plot_name_v =
      ifelse((!combine_pi_ci) | !("ci band" %in% type),
             "pi_intervals","intervals")
    
    plot_list[[plot_name_v]] =
      x$data %>%
      ggplot(aes(x = group,
                 y = !!sym(all.vars(x$formula)[1]))) +
      geom_violin(alpha = 0.2) +
      geom_errorbar(data = newdata,
                    aes(x = !!sym(all.vars(x$formula)[2]),
                        y = `Post Mean`,
                        ymin = PI_lower,
                        ymax = PI_upper),
                    color = "lightsteelblue3") +
      geom_point(data = newdata,
                 aes(x = !!sym(all.vars(x$formula)[2]),
                     y = `Post Mean`),
                 size = 5)
    
    
    
    
    
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      plot_list[["ci_intervals"]] =
        x$data %>%
        ggplot(aes(x = group,
                   y = !!sym(all.vars(x$formula)[1]))) +
        geom_violin(alpha = 0.2)
    }
    
    plot_name_v =
      ifelse((!combine_pi_ci) | !("pi band" %in% type),
             "ci_intervals","intervals")
    
    plot_list[[plot_name_v]] =
      plot_list[[plot_name_v]] +
      geom_errorbar(data = newdata,
                    aes(x = !!sym(all.vars(x$formula)[2]),
                        y = `Post Mean`,
                        ymin = CI_lower,
                        ymax = CI_upper),
                    color = "steelblue4") +
      geom_point(data = newdata,
                 aes(x = !!sym(all.vars(x$formula)[2]),
                     y = `Post Mean`),
                 size = 5)
    
  }
  
  
  # Polish up plots
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    
    for(j in names(plot_list)[grepl("intervals",names(plot_list))]){
      plot_list[[j]] =
        plot_list[[j]] +
        theme_classic() +
        xlab(all.vars(x$formula)[2]) +
        ggtitle(
          paste0(
            ifelse(
              grepl("pi_",j),
              "PI intervals",
              ifelse(grepl("ci_",j),
                     "CI intervals",
                     "CI and PI intervals"
              )
            )
          )
        )
    }
    
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
}






#' @rdname plot
#' @export
plot.lm_b_bma = function(x,
                         type = c("diagnostics",
                                  "pdp",
                                  "ci band",
                                  "pi band"),
                         variable,
                         exemplar_covariates,
                         combine_pi_ci = TRUE,
                         bayes_pvalues_quantiles = c(0.01,1:19/20,0.99),
                         variable_seq_length = 30,
                         return_as_list = FALSE,
                         CI_level = 0.95,
                         PI_level = 0.95){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
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
    
    message("Bayesian p-values measure GOF via \nPr(T(y_obs) - T(y_pred) > 0 | y_obs).\nThus values close to 0.5 are ideal.  Be concerned if values are near 0 or 1.\nThese Bayesian p-values correspond to quantiles of the distribution of y.")
    
    bayes_pvalues_quantiles = sort(bayes_pvalues_quantiles)
    
    preds = predict(x)
    
    T_pred = 
      preds$posterior_draws$ynew %>% 
      apply(1,quantile,probs = bayes_pvalues_quantiles)
    
    T_obs = quantile(x$data[[ all.vars(x$formula)[1] ]],
                     bayes_pvalues_quantiles)
    
    bpvals = 
      rowMeans(T_obs - T_pred > 0)
    
    plot_list$bpvals = 
      tibble(quants = bayes_pvalues_quantiles,
             bpvals  = bpvals) %>% 
      ggplot(aes(x = quants,
                 y = bpvals)) + 
      geom_line() + 
      geom_point() + 
      geom_polygon(data = tibble(x = c(0,1,1,0,0),
                                 y = c(0,0,0.05,0.05,0)),
                   aes(x=x,y=y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.25) +
      geom_polygon(data = tibble(x = c(0,1,1,0,0),
                                 y = c(1,1,0.95,0.95,1)),
                   aes(x=x,y=y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.25) + 
      geom_polygon(data = tibble(x = c(0,1,1,0,0),
                                 y = c(0,0,0.025,0.025,0)),
                   aes(x=x,y=y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.5) +
      geom_polygon(data = tibble(x = c(0,1,1,0,0),
                                 y = c(1,1,0.975,0.975,1)),
                   aes(x=x,y=y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.5) + 
      xlab("Quantiles of outcome") +
      ylab("Bayesian p-values") + 
      theme_minimal()
    
    rm(preds,T_pred)
    
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
    
    message("Partial dependence plots typically require long run times.  Plan accordingly.")
    
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
        newdata$y[i] = mean(temp_preds$newdata$`Post Mean`)
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
                newdata = newdata[[v]],
                CI_level = CI_level,
                PI_level = PI_level)
    }
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pi band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("ci band" %in% type),
                      "pi_band_","band_"),v)
      
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data %>% 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_violin(alpha = 0.2)
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("ci band" %in% type),
                      "pi_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]]$newdata,
                      aes(ymin = PI_lower,
                          ymax = PI_upper),
                      fill = "lightsteelblue3",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]]$newdata,
                    aes(x = !!sym(v),
                        y = `Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]]$newdata,
                        aes(x = !!sym(v),
                            ymin = PI_lower,
                            ymax = PI_upper),
                        color = "lightsteelblue3") +
          geom_point(data = newdata[[v]]$newdata,
                     aes(x = !!sym(v),
                         y = `Post Mean`),
                     size = 5)
      }
      
      
    }
    
    
    
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("ci_band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("ci_band_",v)]] =
            x$data %>% 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_violin(alpha = 0.2)
        }
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pi_ci) | !("pi band" %in% type),
                      "ci_band_","band_"),v)
      
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]]$newdata,
                      aes(ymin = CI_lower,
                          ymax = CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]]$newdata,
                    aes(x = !!sym(v),
                        y = `Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]]$newdata,
                        aes(x = !!sym(v),
                            ymin = CI_lower,
                            ymax = CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]]$newdata,
                     aes(x = !!sym(v),
                         y = `Post Mean`),
                     size = 5)
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    for(v in variable){
      
      for(j in names(plot_list)[grepl("band",names(plot_list)) & grepl(v,names(plot_list))]){
        plot_list[[j]] =
          plot_list[[j]] +
          theme_classic() +
          ggtitle(
            paste0(
              ifelse(
                grepl("pi_",j),
                paste0("PI band for ",v),
                ifelse(grepl("ci_",j),
                       paste0("CI band for ",v),
                       paste0("CI and PI bands for ",v)
                )
              )
            )
          )
      }
      
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
  
}


