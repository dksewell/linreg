library(linreg)
library(ggplot2)
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
  geom_line(aes(y = Estimate), 
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

rm(list = setdiff(ls(),"test_data"))


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
  geom_line(aes(y = Estimate), 
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
preds0a[order(preds0a$y),] |> 
  ggplot(aes(x = y)) +
  geom_ribbon(aes(ymin = PI_lower, 
                  ymax = PI_upper), 
              fill = "lightsteelblue") +
  geom_ribbon(aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  geom_line(aes(y = Estimate), 
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

rm(list=ls())


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
preds0a = 
  predict(fita)
names(preds0a);lapply(preds0a,names)
preds0a$by_group$summary

plot(preds0a$by_group$summary$Estimate,
     as.numeric(by(test_data$y,test_data$x1,mean)))
abline(0,1)

## Make sure information criteria work
AIC(fita)
BIC(fita)
DIC(fita)
WAIC(fita)

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
preds0a = 
  predict(fita)
lapply(preds0a,names)
preds0a$by_group$summary

plot(preds0a$by_group$summary$Estimate,
     as.numeric(by(test_data$y,test_data$x1,mean)))
abline(0,1)

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
preds0a = 
  predict(fita)
lapply(preds0a,names)
preds0a$by_group$summary

plot(preds0a$by_group$summary$Estimate,
     as.numeric(by(test_data$y,test_data$x1,mean)))
abline(0,1)

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
preds0a = 
  predict(fita)
lapply(preds0a,names)
preds0a$by_group$summary

plot(preds0a$by_group$summary$Estimate,
     as.numeric(by(test_data$y,test_data$x1,mean)))
abline(0,1)

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
  geom_line(aes(y = Estimate), 
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
  geom_line(aes(y = Estimate), 
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
  geom_line(aes(y = Estimate), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  theme_minimal()


rm(list = setdiff(ls(),"test_data"))
