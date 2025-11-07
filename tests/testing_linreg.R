# Add seeds to anything stochastic!
# Test bma
# Test mediation
# Add BF to determine hetero or homo for aov_b
# Add glm's via IS
# Add SUBSET
# Maybe add AR(p) 

library(linreg)
library(future)


# Helper functions --------------------------------------------------------

beta_parms = find_beta_parms(2/5,0.68,0.9)
beta_parms[1] / sum(beta_parms)
pbeta(0.68,beta_parms[1],beta_parms[2])

igamma_parms = find_invgamma_parms(10 / 4,5,0.75)
extraDistr::pinvgamma(10/4,igamma_parms[1] * 0.5,igamma_parms[2] * 0.5)
extraDistr::pinvgamma(4.1,igamma_parms[1] * 0.5,igamma_parms[2] * 0.5)

rm(list=ls())



# Linear regression -------------------------------------------------------

# Create data

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$y = 
  rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )


# Conjugate prior
## Make sure CI_level works (and print.lm_b works)
fita = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "conj")
fitb = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "conj",
       CI_level = 0.9)
fita
fitb

## Make sure summary.lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prior hyperparameters work
fitc = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "conj",
       prior_beta_mean = rep(2,7),
       prior_beta_precision = diag(7) * 0.1,
       prior_var_shape = 2,
       prior_var_rate = 2,
       CI_leve = 0.9)
fita
fitc

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = PI_lower, 
                  ymax = PI_upper), 
              fill = "lightsteelblue") +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

## Make sure Savage-Dickey ratio works
SDratio(fita)

## Make sure plotting function works
plot(fita,
     type = "diagnostics")
plot(fita,
     type = "pdp")
test = 
  plot(fita,
       type = "pdp",
       return_as_list = TRUE)
patchwork::wrap_plots(test)
plot(fita,
     type = "pdp",
     variable = "x1")
plot(fita,
     type = c("ci","pi"),
     variable = "x1",
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE)
plot(fita,
     type = "pi",
     variable = "x1")
plot(fita,
     type = "pi",
     variable = "x1",
     PI_level = 0.5)
plot(fita,
     type = "pi")
plot(fita,
     type = "pi",
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = "ci",
     variable = "x1")
plot(fita,
     type = "ci")
plot(fita,
     type = "ci",
     CI_level = 0.999)
plot(fita)




# Zellner-g prior
## Make sure CI_level works (and print.lm_b works)
fita = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "zellner")
fitb = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "zellner",
       CI_level = 0.9)
fita
fitb

## Make sure summary.lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prior hyperparameters work
fitc = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "zellner",
       zellner_g = 10)
fita
fitc

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = PI_lower, 
                  ymax = PI_upper), 
              fill = "lightsteelblue") +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

## Make sure Savage-Dickey ratio works
SDratio(fita)

## Make sure plotting function works
plot(fita,
     type = "diagnostics")
plot(fita,
     type = "pdp")
test = 
  plot(fita,
       type = "pdp",
       return_as_list = TRUE)
patchwork::wrap_plots(test)
plot(fita,
     type = "pdp",
     variable = "x1")
plot(fita,
     type = c("ci","pi"),
     variable = "x1",
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE)
plot(fita,
     type = "pi",
     variable = "x1")
plot(fita,
     type = "pi")
plot(fita,
     type = "pi",
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = "ci",
     variable = "x1")
plot(fita,
     type = "ci")


rm(list = setdiff(ls(),"test_data"))


# Improper prior
## Make sure CI_level works (and print.lm_b works)
fita = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "impr")
fitb = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "impr",
       CI_level = 0.9)
fita
fitb

## Make sure summary.lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prior hyperparameters DO NOT work
fitc = 
  lm_b(y ~ x1 + x2 + x3,
       data = test_data,
       prior = "impr",
       zellner_g = 5,
       prior_beta_mean = rep(2,7),
       prior_beta_precision = diag(7) * 0.1,
       prior_var_shape = 2,
       prior_var_rate = 2,
       CI_leve = 0.95)
all.equal(summary(fita),
          summary(fitc))

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
predict(fita,
        CI_level = 0.8) |> 
  head()
predict(fita,
        PI_level = 0.8) |> 
  head()


preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = PI_lower, 
                  ymax = PI_upper), 
              fill = "lightsteelblue") +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

## Make sure Savage-Dickey ratio DOES NOT works
SDratio(fita)

## Make sure plotting function works
plot(fita,
     type = "diagnostics")
plot(fita,
     type = "pdp")
test = 
  plot(fita,
       type = "pdp",
       return_as_list = TRUE)
patchwork::wrap_plots(test)
plot(fita,
     type = "pdp",
     variable = "x1")
plot(fita,
     type = c("ci","pi"),
     variable = "x1",
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE,
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE)
plot(fita,
     type = "pi",
     variable = "x1")
plot(fita,
     type = "pi")
plot(fita,
     type = "pi",
     exemplar_covariates = fita$data[1,])
plot(fita,
     type = "ci",
     variable = "x1")
plot(fita,
     type = "ci")

rm(list=ls())




# Check BMA ---------------------------------------------------------------


# Create data
pacman::p_load(coda,
               dplyr,
               extraDistr,
               magrittr,
               mvtnorm,
               future,
               future.apply,
               ggplot2,
               patchwork)

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             x4 = rnorm(N),
             x5 = rnorm(N),
             x6 = rnorm(N),
             x7 = rnorm(N),
             x8 = rnorm(N),
             x9 = rnorm(N),
             x10 = rnorm(N))
test_data$y = 
  rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )

formula = y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
data = test_data
zellner_g = nrow(data)
CI_level = 0.95
mc_draws = 5e4






# ANOVA -------------------------------------------------------------------

# Create data

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rep(letters[1:5],N/5))
test_data$y = 
  rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )


# Proper, heteroscedastic

## Make sure CI_level works (and print.lm_b works)
fita = 
  aov_b(y ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01)

fitb = 
  aov_b(y ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        CI_level = 0.9)
fita
fitb

## Make sure summary.aov_b works
summary(fita)

## Make sure coef.aov_b works
coef(fita)

## Make sure prior hyperparameters work
fitc = 
  aov_b(y ~ x1,
        test_data,
        prior_mean_mu = 200,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        CI_level = 0.9)
fita
fitc

## Make sure prediction function works
predict(fita)
predict(fita,
        CI_level = 0.8)
predict(fita,
        PI_level = 0.8)

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

## Make sure plotting function works
plot(fita,
     type = "diagnostics")
plot(fita,
     type = c("ci","pi"),
     combine_pi_ci = TRUE)
plot(fita,
     type = "pi")
plot(fita,
     type = "pi",
     PI_level = 0.8)
plot(fita,
     type = "ci")
plot(fita,
     type = "ci",
     CI_level = 0.999)
plot(fita)


## Make sure parallelization works
plan(multisession,workers = 5)
fitd = 
  aov_b(y ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01)
fitd
fita
plan(sequential)

rm(list = setdiff(ls(),"test_data"))


# Proper, homoscedastic

## Make sure CI_level works (and print.lm_b works)
fita = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01)

fitb = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        CI_level = 0.8)
fita
fitb

## Make sure summary.aov_b works
summary(fita)

## Make sure coef.aov_b works
coef(fita)

## Make sure prior hyperparameters work
fitc = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        prior_mean_mu = 200,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        CI_level = 0.9)
fita
fitc

## Make sure prediction function works
predict(fita)
predict(fita,
        CI_level = 0.8)
predict(fita,
        PI_level = 0.8)

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

rm(list = setdiff(ls(),"test_data"))



# Improper, heteroscedastic

## Make sure CI_level works (and print.lm_b works)
fita = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = TRUE,
        improper = TRUE)

fitb = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = TRUE,
        improper = TRUE,
        CI_level = 0.8)
fita
fitb

## Make sure summary.aov_b works
summary(fita)

## Make sure coef.aov_b works
coef(fita)

## Make sure prior hyperparameters DO NOT work
fitc = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = TRUE,
        prior_mean_mu = 200,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        improper = TRUE)
all.equal(summary(fita)$summary,
          summary(fitc)$summary)

## Make sure prediction function works
predict(fita)
predict(fita,
        CI_level = 0.8)
predict(fita,
        PI_level = 0.8)


## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

rm(list = setdiff(ls(),"test_data"))



# Improper, homoscedastic

## Make sure CI_level works (and print.lm_b works)
fita = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        improper = TRUE)

fitb = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        improper = TRUE,
        CI_level = 0.8)
fita
fitb

## Make sure summary.aov_b works
summary(fita)

## Make sure coef.aov_b works
coef(fita)

## Make sure prior hyperparameters DO NOT work
fitc = 
  aov_b(y ~ x1,
        test_data,
        heteroscedastic = FALSE,
        prior_mean_mu = 200,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        improper = TRUE)
all.equal(summary(fita)$summary,
          summary(fitc)$summary)

## Make sure prediction function works
predict(fita)
predict(fita,
        CI_level = 0.8)
predict(fita,
        PI_level = 0.8)

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

rm(list = ls())





# Test loss-likelihood approach (Gaussian) --------------------------------

# Create data

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$y = 
  rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )


# Bootstrapping approach - sequential
plan(sequential)
## Make sure CI_level works (and print.lm_b works)
fita = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          n_draws = 100,
          seed = 2025)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          n_draws = 100,
          seed = 2025,
          CI_level = 0.8)
fita
fitb

## Make sure summary.np_lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()


rm(list = setdiff(ls(),"test_data"))




# Bootstrapping approach - parallel
plan(multisession, workers = 10)
## Make sure future.seed is reproducible, CI_level works (and print.lm_b works)
fita = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          n_draws = 100,
          seed = 2025)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          n_draws = 100,
          seed = 2025)
fitc = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          n_draws = 100,
          seed = 2025,
          CI_level = 0.8)
all.equal(fita$summary,
          fitb$summary)
fita
fitc

## Make sure summary.np_lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()

## Check that parallelization is speeding things up

system.time({
  fitd = 
    np_lm_b(y ~ x1 + x2 + x3,
            data = test_data,
            family = gaussian(),
            n_draws = 500,
            seed = 2025)
})
# user  system elapsed 
# 11.76    6.28   23.77 

plan(sequential)
system.time({
  fite = 
    np_lm_b(y ~ x1 + x2 + x3,
            data = test_data,
            family = gaussian(),
            n_draws = 500,
            seed = 2025)
})
# user  system elapsed 
# 6.84    0.07    7.00 # Eh?  Not sure what's happening.
fitd
fite

rm(list = setdiff(ls(),"test_data"))




# large sample approach
## Make sure CI_level works (and print.lm_b works)
fita = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          seed = 2025)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          seed = 2025,
          CI_level = 0.8)
fita
fitb

## Make sure summary.np_lm_b works
summary(fita)

## Make sure coef.lm_b works
coef(fita)

## Make sure prediction function works
preds0a = 
  predict(fita)
head(preds0a)
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()


rm(list = setdiff(ls(),"test_data"))
