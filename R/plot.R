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
#' @param ... optional arguments.
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
                     variable_seq_length = 30,
                     return_as_list = FALSE,
                     CI_level = 0.95,
                     PI_level = 0.95,
                     ...){
  
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
  
  if(missing(variable)){
    variable = 
      x$formula |> 
      terms() |> 
      delete.response() |> 
      all.vars() |> 
      unique()
  }
  
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
                    size = 3)
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
                     x$data) |> 
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
          x$data |> 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data |> 
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
                     size = 3)
      }
    
      
    }
      
      
      
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("ci_band_",v)]] =
            x$data |> 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("ci_band_",v)]] =
            x$data |> 
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
                    size = 3)
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
                      PI_level = 0.95,
                      ...){
  
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
      x$data |>
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
                 size = 3)
    
    
    
    
    
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      plot_list[["ci_intervals"]] =
        x$data |>
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
                 size = 3)
    
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
                         PI_level = 0.95,
                         seed = 1,
                         ...){
  
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
  
  if(missing(variable)) variable = attr(terms(x$formula),"term.labels")
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    set.seed(seed)
    
    message("Bayesian p-values measure GOF via \nPr(T(y_obs) - T(y_pred) > 0 | y_obs).\nThus values close to 0.5 are ideal.  Be concerned if values are near 0 or 1.\nThese Bayesian p-values correspond to quantiles of the distribution of y.")
    
    bayes_pvalues_quantiles = sort(bayes_pvalues_quantiles)
    
    preds = predict(x)
    
    T_pred = 
      preds$posterior_draws$ynew |> 
      apply(1,quantile,probs = bayes_pvalues_quantiles)
    
    T_obs = quantile(x$data[[ all.vars(x$formula)[1] ]],
                     bayes_pvalues_quantiles)
    
    bpvals = 
      rowMeans(T_obs - T_pred > 0)
    
    plot_list$bpvals = 
      tibble(quants = bayes_pvalues_quantiles,
             bpvals  = bpvals) |> 
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
                     size = 3)
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
                     x$data) |> 
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
          x$data |> 
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data |> 
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
                     size = 3)
      }
      
      
    }
    
    
    
  }
  
  if("ci band" %in% type){
    
    # Get starter plots if !combine_pi_ci
    if( (!combine_pi_ci) | !("pi band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("ci_band_",v)]] =
            x$data |> 
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("ci_band_",v)]] =
            x$data |> 
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
                     size = 3)
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
plot.glm_b = function(x,
                      type,
                      variable,
                      exemplar_covariates,
                      combine_pi_ci = TRUE,
                      variable_seq_length = 30,
                      return_as_list = FALSE,
                      CI_level = 0.95,
                      PI_level = 0.95,
                      seed = 1,
                      ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  if(missing(type)){
    type = 
      c("diagnostics",
        "pdp",
        "ci band")
    if(x$family$family != "binomial") type = c(type,"pi band")
  }
  
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
  
  if( (x$family$family == "binomial") & 
      ("pi band" %in% type) ){
    type = setdiff(type,"pi band")
    if(length(type) == 0){
      warning("Prediction band cannot be supplied for a binomial outcome.\nResults shown will be credible band instead.")
      type = "ci band"
    }
  }
  
  if(missing(variable)) variable = attr(terms(x$formula),"term.labels")
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    # Extract 
    mframe = model.frame(x$formula, x$data)
    y = model.response(mframe)
    X = model.matrix(x$formula,x$data)
    os = model.offset(mframe)
    N = nrow(X)
    p = ncol(X)
    if(is.null(os)) os = numeric(N)
    
    
    message("Bayesian p-values measure GOF via \nPr(T(y_obs) - T(y_pred) > 0 | y_obs).\nThus values close to 0.5 are ideal.  Be concerned if values are near 0 or 1.\nThis Bayesian p-value corresponds to the deviance.")
    
    if("posterior_covariance" %in% names(x)){
      
      # Get posterior draws of E(y)
      yhat_draws = 
        x$trials * 
        x$family$linkinv(os + tcrossprod(X, rmvnorm(5e3,
                                                    x$summary$`Post Mean`,
                                                    x$posterior_covariance)))
      
      # Get posterior draws of y
      if(x$family$family == "binomial"){
        y_draws = 
          future_sapply(1:nrow(yhat_draws),
                        function(i){
                          rbinom(ncol(yhat_draws),
                                 x$trials[i],
                                 yhat_draws[i,])
                        },
                        future.seed = seed)
        
        deviances_pred = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dbinom(y_draws[draw,],
                                       x$trials,
                                       yhat_draws[,draw],
                                       log = TRUE))
                        })
        deviances_obs = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dbinom(y,
                                       x$trials,
                                       yhat_draws[,draw],
                                       log = TRUE))
                        })
        
        
      }
      if(x$family$family == "poisson"){
        y_draws = 
          future_sapply(1:nrow(yhat_draws),
                        function(i){
                          rpois(ncol(yhat_draws),yhat_draws[i,])
                        },
                        future.seed = seed)
        deviances_pred = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dpois(y_draws[draw,],
                                      yhat_draws[,draw],
                                      log = TRUE))
                        })
        deviances_obs = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dpois(y,
                                      yhat_draws[,draw],
                                      log = TRUE))
                        })
      }
      
      
    }else{#End: Getting pvals for large sample approx
      
      # Get posterior draws of E(y)
      yhat_draws = 
        x$trials * 
        x$family$linkinv(os + tcrossprod(X, x$proposal_draws))
      
      # Get posterior draws of y
      if(x$family$family == "binomial"){
        y_draws = 
          future_sapply(1:nrow(yhat_draws),
                        function(i){
                          rbinom(ncol(yhat_draws),
                                 x$trials[i],
                                 yhat_draws[i,])
                        },
                        future.seed = seed)
        
        deviances_pred = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dbinom(y_draws[draw,],
                                       x$trials,
                                       yhat_draws[,draw],
                                       log = TRUE))
                        })
        deviances_obs = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dbinom(y,
                                       x$trials,
                                       yhat_draws[,draw],
                                       log = TRUE))
                        })
      }
      if(x$family$family == "poisson"){
        y_draws = 
          future_sapply(1:nrow(yhat_draws),
                        function(i){
                          rpois(ncol(yhat_draws),yhat_draws[i,])
                        },
                        future.seed = seed)
        deviances_pred = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dpois(y_draws[draw,],
                                      yhat_draws[,draw],
                                      log = TRUE))
                        })
        deviances_obs = 
          future_sapply(1:nrow(y_draws),
                        function(draw){
                          -2.0 * 
                            sum(dpois(y,
                                      yhat_draws[,draw],
                                      log = TRUE))
                        })
        
      }
      
      resample_index = 
        sample(1:length(deviances_obs),length(deviances_obs),TRUE,x$importance_sampling_weights)
      deviances_obs = deviances_obs[resample_index]
      deviances_pred = deviances_pred[resample_index]
      
    }#End: Getting pvals for IS
    
    
    dx_data = 
      tibble(T_obs = deviances_obs,
             T_pred = deviances_pred) |> 
      mutate(obs_gr_pred = T_obs > T_pred)
    
    plot_list$bpvals = 
      dx_data |> 
      ggplot(aes(x = T_pred,
                 y = T_obs,
                 color = obs_gr_pred)) + 
      geom_point() + 
      geom_abline(intercept = 0,
                  slope = 1) + 
      xlab(bquote(T(y[pred] * "," * beta))) +
      ylab(bquote(T(y[obs] * "," * beta))) +
      theme_classic() +
      scale_color_viridis_d() +
      ggtitle(paste0("Bayesian p-value based on deviance = ",
                     round(mean(deviances_obs > deviances_pred),3))) + 
      theme(legend.position = "none")
    
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
      suppressMessages({
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
      })
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          x$data |>
          ggplot(aes(x = !!sym(variable[v]),
                     y = !!sym(all.vars(x$formula)[1]))) + 
          geom_point(alpha = 0.2) +
          geom_line(data = newdata,
                    aes(x = var_of_interest,
                        y = y))
      }else{
        if(x$family$family == "poisson"){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            ggplot(aes(x = !!sym(variable[v]),
                       y = !!sym(all.vars(x$formula)[1]))) + 
            geom_violin(alpha = 0.2)
        }
        if(x$family$family == "binomial"){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            group_by(get(variable[v])) |> 
            summarize(prop1 = mean(near(!!sym(all.vars(x$formula)[1]), 1))) |> 
            ggplot(aes(x = `get(variable[v])`,
                       y = prop1)) + 
            geom_col(fill="gray70")
        }
        
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = var_of_interest,
                         y = y),
                     size = 3)
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
                     x$data) |> 
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
      
      suppressMessages({
        newdata[[v]] = 
          predict(x,
                  newdata = newdata[[v]],
                  CI_level = CI_level,
                  PI_level = PI_level)
      })
      
      
      # Get starter plots
      two_plots = 
        (!combine_pi_ci) &
        ( ("ci band" %in% type) & ("pi band" %in% type) )
      if(two_plots){
        plot_name_v1 = 
            paste0("pi_band_",v)
        plot_name_v2 = 
          paste0("ci_band_",v)
        
        if(is.numeric(x$data[[v]])){
          plot_list[[plot_name_v1]] =
            plot_list[[plot_name_v2]] =
            x$data |>
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) + 
            geom_point(alpha = 0.2)
        }else{
          if(x$family$family == "poisson"){
            plot_list[[plot_name_v1]] =
              plot_list[[plot_name_v2]] =
              x$data |>
              ggplot(aes(x = !!sym(v),
                         y = !!sym(all.vars(x$formula)[1]))) + 
              geom_violin(alpha = 0.2)
          }
          if(x$family$family == "binomial"){
            plot_list[[plot_name_v1]] =
              plot_list[[plot_name_v2]] =
              x$data |>
              group_by(get(v)) |> 
              summarize(prop1 = mean(near(!!sym(all.vars(x$formula)[1]), 1))) |> 
              rename(!!v := `get(v)`) |> 
              ggplot(aes(x = !!sym(v),
                         y = prop1)) + 
              geom_col(fill="gray70")
          }
        }
        
      }else{#End: starting two plots for bands/intervals
        
        band_to_plot = 
          paste0(gsub("\ ","_",type[grep("band",type)]),
                 "_",v)
        if(length(band_to_plot) == 2) band_to_plot = paste0("band_",v)
        
        
        for( plot_name_v in band_to_plot){
          if(is.numeric(x$data[[v]])){
            plot_list[[plot_name_v]] =
              x$data |>
              ggplot(aes(x = !!sym(v),
                         y = !!sym(all.vars(x$formula)[1]))) + 
              geom_point(alpha = 0.2)
          }else{
            if(x$family$family == "poisson"){
              plot_list[[plot_name_v]] =
                x$data |>
                ggplot(aes(x = !!sym(v),
                           y = !!sym(all.vars(x$formula)[1]))) + 
                geom_violin(alpha = 0.2)
            }
            if(x$family$family == "binomial"){
              plot_list[[plot_name_v]] =
                x$data |>
                group_by(get(v)) |> 
                summarize(prop1 = mean(near(!!sym(all.vars(x$formula)[1]), 1))) |> 
                rename(!!v := `get(v)`) |> 
                ggplot(aes(x = !!sym(v),
                           y = prop1)) + 
                geom_col(fill="gray70")
            }
          }
        }
        
      }
      
      
    }#End: for loop through variables
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pi band" %in% type){
    
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
                     size = 3)
      }
      
      
    }
    
    
    
  }
  
  if("ci band" %in% type){
    
    
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
          geom_errorbar(data = 
                          newdata[[v]] |> 
                          mutate(prop1 = 0.0), # Stupid hack to make ggplot work right.
                        aes(x = !!sym(v),
                            ymin = CI_lower,
                            ymax = CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                     aes(x = !!sym(v),
                         y = `Post Mean`),
                     size = 3)
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
plot.np_glm_b = function(x,
                         type,
                         variable,
                         exemplar_covariates,
                         variable_seq_length = 30,
                         return_as_list = FALSE,
                         CI_level = 0.95,
                         seed = 1,
                         ...){
  
  alpha_ci = 1.0 - CI_level
  
  if(missing(type)){
    type = 
      c("pdp",
        "ci band")
  }
  
  type = c("pdp",
           "ci band")[pmatch(tolower(type),
                             c("pdp",
                               "ci band"))]
  
  
  if(missing(variable)) variable = attr(terms(x$formula),"term.labels")
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  
  # Get unique values and x sequences for plots
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

  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble(var_of_interest = x_seq[[v]],
               y = 0.0)
      suppressMessages({
        for(i in 1:length(x_seq[[v]])){
          temp_preds = 
            predict(x,
                    newdata = 
                      x$data |>
                      dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]),
                    CI_level = CI_level)
          newdata$y[i] = mean(temp_preds$`Post Mean`)
        }
      })
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          x$data |>
          ggplot(aes(x = !!sym(variable[v]),
                     y = !!sym(all.vars(x$formula)[1]))) + 
          geom_point(alpha = 0.2) +
          geom_line(data = newdata,
                    aes(x = var_of_interest,
                        y = y))
      }else{
        if(x$family$family == "binomial"){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            group_by(get(variable[v])) |> 
            summarize(prop1 = mean(near(!!sym(all.vars(x$formula)[1]), 1))) |> 
            ggplot(aes(x = `get(variable[v])`,
                       y = prop1)) + 
            geom_col(fill="gray70")
        }else{
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            ggplot(aes(x = !!sym(variable[v]),
                       y = !!sym(all.vars(x$formula)[1]))) + 
            geom_violin(alpha = 0.2)
        }
        
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = var_of_interest,
                         y = y),
                     size = 3)
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
  if("ci band" %in% type){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) |> 
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
      
      suppressMessages({
        newdata[[v]] = 
          predict(x,
                  newdata = newdata[[v]],
                  CI_level = CI_level)
      })
      
      
      # Get starter plot
      plot_name_v = paste0("ci_band_",v)
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data |>
          ggplot(aes(x = !!sym(v),
                     y = !!sym(all.vars(x$formula)[1]))) + 
          geom_point(alpha = 0.2) +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = CI_lower,
                          ymax = CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = !!sym(v),
                        y = `Post Mean`))
      }else{
        if(x$family$family == "binomial"){
          plot_list[[plot_name_v]] =
            x$data |>
            group_by(get(v)) |> 
            summarize(prop1 = mean(near(!!sym(all.vars(x$formula)[1]), 1))) |> 
            rename(!!v := `get(v)`) |> 
            ggplot(aes(x = !!sym(v),
                       y = prop1)) + 
            geom_col(fill="gray70")
        }else{
          plot_list[[plot_name_v]] =
            x$data |>
            ggplot(aes(x = !!sym(v),
                       y = !!sym(all.vars(x$formula)[1]))) + 
            geom_violin(alpha = 0.2)
        }
        
        
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = 
                          newdata[[v]] |> 
                          mutate(prop1 = 0.0), # Stupid hack to make ggplot work right.
                        aes(x = !!sym(v),
                            ymin = CI_lower,
                            ymax = CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                     aes(x = !!sym(v),
                         y = `Post Mean`),
                     size = 3)
        
      }
    
      plot_list[[plot_name_v]] =
        plot_list[[plot_name_v]] +
        theme_classic() +
        ggtitle(paste0("CI band for ",v))
      
    }#End: loop through variables
    
  }#End: CI band code
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
  
}



