% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lm_b.R
\name{lm_b}
\alias{lm_b}
\title{Bayesian Linear Models}
\usage{
lm_b(
  formula,
  data,
  prior = c("zellner", "conjugate", "improper")[1],
  zellner_g,
  prior_beta_mean,
  prior_beta_precision,
  prior_var_shape = 0.001,
  prior_var_rate = 0.001,
  CI_level = 0.95
)
}
\arguments{
\item{formula}{A formula specifying the model.}

\item{data}{A data frame in which the variables specified in the formula
will be found. If missing, the variables are searched for in the standard way.}

\item{prior}{character.  One of "zellner", "conjugate", or "improper", giving
the type of prior used on the regression coefficients.}

\item{zellner_g}{numeric.  Positive number giving the value of "g" in Zellner's
g prior.  Ignored unless prior = "zellner".}

\item{prior_beta_mean}{numeric vector of same length as regression coefficients
(denoted p). Unless otherwise specified, automatically set to rep(0,p).  Ignored
unless prior = "conjugate".}

\item{prior_beta_precision}{pxp matrix giving a priori precision matrix to be
scaled by the residual precision.}

\item{prior_var_shape}{numeric. Twice the shape parameter for the inverse gamma prior on
the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).}

\item{prior_var_rate}{numeric. Twice the rate parameter for the inverse gamma prior on
the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).}

\item{CI_level}{numeric. Credible interval level.}
}
\value{
lm_b returns an object of class "lm_b", which behaves as a list with
the following elements:
\itemize{
\item summary - data.frame giving results for regression coefficients
\item post_parms - list giving the posterior parameters
\item hyperparms - list giving the user input (or default) hyperparameters used
\item fitted - posterior mean of the individuals' means
\item residuals - posterior mean of the residuals
\item formula, data - input by user
}
}
\description{
lm_b is used to fit linear models.  It can be used to carry out
regression, single stratum analysis of variance and analysis of
covariance (although aov_b may provide a more convenient
interface for ANOVA.)
}
