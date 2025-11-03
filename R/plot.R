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
  
  plot_list = list()
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
    
  }
  if(type == "pdp"){
    
    if(missing(variable)) variable = all.vars(x$formula)[-1]
    
    for(v in 1:length(variable)){
      N = nrow(x$data)
      x_unique = 
        unique(x$data[[variable[v]]])
      if(length(x_unique) > variable_seq_length){
        variable_seq_length_v = variable_seq_length
        xseq = 
          seq(min(x$data[[variable[v]]]),
              max(x$data[[variable[v]]]),
              l = variable_seq_length_v)
      }else{
        xseq = factor(sort(x_unique),levels = sort(x_unique))
        variable_seq_length_v = length(x_unique)
      }
      
      newdata = 
        tibble(var_of_interest = xseq,
               y = numeric(variable_seq_length_v))
      for(i in 1:variable_seq_length_v){
        temp_preds = 
          predict(x,
                  newdata = 
                    x$data |>
                    dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]))
        newdata$y[i] = mean(temp_preds$Estimate)
      }
      
      plot_list[[variable[v]]] = 
        x$data |>
        ggplot(aes(x = !!sym(variable[v]),
                   y = !!sym(all.vars(x$formula)[1]))) + 
        geom_point(alpha = 0.2)
      if(variable_seq_length_v == variable_seq_length){
        plot_list[[variable[v]]] = 
          plot_list[[variable[v]]] + 
          geom_line(data = newdata,
                    aes(x = var_of_interest,
                        y = y))
      }else{
        plot_list[[variable[v]]] = 
          plot_list[[variable[v]]] + 
          geom_point(data = newdata,
                    aes(x = var_of_interest,
                        y = y),
                    size = 5)
      }
      
      plot_list[[variable[v]]] = 
        plot_list[[variable[v]]] + 
        xlab(variable[variable[v]]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }
  
  
  if( ("pi band" %in% type) | ("ci band" %in% type) ){
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) %>% 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
  }
  
  
  if("pi band" %in% type){
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      
    }
    
  }
  
  if("ci band" %in% type){
    if(missing(exemplar_covariates)) stop("Must provide other covariate values.")
    
  }
  
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      eval(parse(text = paste(
        paste("plot_list[['",
              names(plot_list),
              "']]"),
      collapse = "+")))
      )
  }
  
  
}

