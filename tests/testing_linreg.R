# Change all outcome names from y to something else.  Things will break.
# Set summary of glm_b and np_glm_b to report in exponentiated scale.
# Set summary of glm_b to report ESS
# Add generics for glm_b: plot
# Add ROPE functionality and add bounds for 
#   np_lm_b
# Add plot.np_glm_b
# Add IC for glm_b object
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

# Check to make sure print.lm_b_bma and coef.lm_b_bma works
fita = 
  bma_inference(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
                test_data)
fita
coef(fita)

# Make sure parallelization works. 
plan(multisession, workers = 10)
fitb = 
  bma_inference(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
                test_data,
                mc_draws = 5e3)
plan(sequential)
fitb

## Make sure summary.lm_b_bma works
summary(fita)
summary(fita,
        CI_level = 0.8)


# Test predictions
preds = 
  predict(fita)
colnames(preds$newdata)
dim(preds$posterior_draws$mean_of_ynew)
dim(preds$posterior_draws$ynew)


# Test plotting
plot(fita,
     "dx")

plot(fita,
     "pdp",
     variable = c("x1","x2","x3"))

plot(fita,
     "ci",
     variable = c("x1","x2","x3"))
plot(fita,
     "ci",
     variable = c("x1","x2","x3"),
     CI_level = 0.99)
plot(fita,
     "pi",
     variable = c("x1","x2","x3"))
plot(fita,
     "pi",
     variable = c("x1","x2","x3"),
     PI_level = 0.5)
plot(fita,
     c("ci","pi"),
     variable = c("x1","x2","x3"))
plot(fita,
     c("ci","pi"),
     combine_pi_ci = FALSE,
     variable = c("x1","x2","x3"))



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



# Mediation ---------------------------------------------------------------

# Binary treatment
set.seed(2025)
N = 500
test_data = 
  data.frame(tr = rep(0:1,N/2),
             x1 = rnorm(N))
test_data$m = 
  rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
test_data$y = 
  rnorm(N,-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1)

m1 = 
  lm_b(m ~ tr + x1,
       data = test_data)
m2 = 
  lm_b(y ~ m + tr + x1,
       data = test_data)
m3 = 
  mediate_b(m1,m2,
            treat = "tr",
            control_value = 0,
            treat_value = 1)


# Continuous treatment
set.seed(2025)
N = 500
test_data = 
  data.frame(tr = rnorm(N),
             x1 = rnorm(N))
test_data$m = 
  rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
test_data$y = 
  rnorm(N,-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1)

m1 = 
  lm_b(m ~ tr + x1,
       data = test_data)
m2 = 
  lm_b(y ~ m + tr + x1,
       data = test_data)
m3 = 
  mediate_b(m1,m2,
            treat = "tr",
            control_value = 0,
            treat_value = 1)
m3


# Test glm_b (binary) -----------------------------------------------------

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$y = 
  rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))

table(test_data$y,test_data$x3) |>
  prop.table(2)
boxplot(x1 ~ y,test_data)

# IS
## Zellner
fita =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        seed = 2025)
fitb =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        seed = 2025,
        CI_level = 0.8)
fita
fitb
coef(fita)
summary(fita)
summary(fita,
        CI_level = 0.8)
summary(fita,
        interpretable_scale = FALSE)
preds = predict(fita)
boxplot(`Post Mean` ~ y, data = preds)



## Normal
fitc =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        prior = "normal",
        seed = 2025)
fita
fitc
coef(fitc)
summary(fitc)

## improper
fitc =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        prior = "improper")
fita
fitc
coef(fitc)
summary(fitc)


# Large sample approx
## Zellner
fitb =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        n_draws = NA,
        seed = 2025)
fitc =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        n_draws = NA,
        seed = 2025,
        CI_level = 0.8)
fita
fitb
fitc
coef(fitb)
summary(fitb)
preds = predict(fitb)
boxplot(`Post Mean` ~ y, data = preds)


## Normal
fitd =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        prior = "normal",
        n_draws = NA,
        seed = 2025)
fitd
coef(fitb)
summary(fitb)

## improper
fite =
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        n_draws = NA,
        prior = "improper")
fite
coef(fite)
summary(fite)

rm(list=ls())

# Test glm_b (Poisson) ----------------------------------------------------

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$y = 
  rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)


# IS
## Zellner
fita =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        seed = 2025)
fitb =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        seed = 2025,
        CI_level = 0.8)
fita
fitb
coef(fita)
summary(fita)
summary(fita,
        CI_level = 0.8)
summary(fita,
        interpretable_scale = FALSE)
preds = predict(fita)
colnames(preds)
plot(`Post Mean` ~ y, data = preds %>% arrange(y))


## Normal
fitc =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        prior = "normal",
        seed = 2025)
fita
fitc
coef(fitc)
summary(fitc)

## improper
fitc =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        prior = "improper")
fita
fitc
coef(fitc)
summary(fitc)


# Large sample approx
## Zellner
fitb =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        n_draws = NA,
        seed = 2025)
fitc =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        n_draws = NA,
        seed = 2025,
        CI_level = 0.8)
fita
fitb
fitc
coef(fitb)
summary(fitb)
preds = predict(fitb)
str(preds)
plot(`Post Mean` ~ y, data = preds %>% arrange(y))



## Normal
fitd =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        prior = "normal",
        n_draws = NA,
        seed = 2025)
fitd
coef(fitb)
summary(fitb)

## improper
fite =
  glm_b(y ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        n_draws = NA,
        prior = "improper")
fite
coef(fite)
summary(fite)



# Test glm_b (Gaussian) ----------------------------------------------------


set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$y = 
  rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )

fita = 
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        prior = "normal",
        family = "gaussian")
fitb = 
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        prior = "zellner",
        family = "gaussian")
fitb = 
  glm_b(y ~ x1 + x2 + x3,
        data = test_data,
        prior = "zellner",
        zellner_g = 10,
        family = "gaussian")


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
  ggplot2::ggplot(ggplot2::aes(x = y)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_lower, 
                  ymax = CI_upper), 
              fill = "steelblue3") +
  ggplot2::geom_line(ggplot2::aes(y = `Post Mean`), 
            color = "steelblue4", 
            linewidth = 1, 
            linetype = "solid",
            alpha = 0.5) +
  ggplot2::theme_minimal()


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

plan(sequential)


rm(list = setdiff(ls(),c("fita","test_data")))




# large sample approach
## Make sure CI_level works (and print.lm_b works)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          seed = 2025)
fitc = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          seed = 2025,
          CI_level = 0.8)
fita
fitb
fitc

## Make sure summary.np_lm_b works
summary(fitb)

## Make sure coef.lm_b works
coef(fitb)

## Make sure prediction function works
preds0a = 
  predict(fitb)
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


# Test GLS
## Bootstrapping
fitd = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          loss = "gls",
          n_draws = 100,
          seed = 2025)
fitd
summary(fitd)
coef(fitd)
preds0d = 
  predict(fitd)
head(preds0d)

## large sample approx
fite = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = gaussian(),
          loss = "gls",
          seed = 2025)
fitd
fite
summary(fite)
coef(fite)
preds0e = 
  predict(fite)
head(preds0e)




rm(list = ls())





# Test loss-likelihood approach (Binomial) --------------------------------

# Create data

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$y = 
  rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
table(test_data$y,test_data$x3) |>
  prop.table(2)
boxplot(x1 ~ y,test_data)
# 
# formula = y ~ x1 + x2 + x3
# data = test_data
# family = binomial()
# n_draws = 10
# loss = "selfinformation"
# CI_level = 0.95

# Bootstrapping approach - sequential
plan(sequential)
## Make sure CI_level works (and print.lm_b works)
fita = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          n_draws = 500,
          seed = 2025)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          n_draws = 100,
          trials = rep(1,N),
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


rm(list = setdiff(ls(),c("fita","test_data")))



# large sample approach
## Make sure CI_level works (and print.lm_b works)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          seed = 2025)
fitc = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          seed = 2025,
          CI_level = 0.8)
fita
fitb
fitc

## Make sure summary.np_lm_b works
summary(fitb)

## Make sure coef.lm_b works
coef(fitb)

## Make sure prediction function works
preds0b = 
  predict(fitb)
head(preds0b)


# Make sure parallelization works
plan(multisession,
     workers = 10)
## Make sure CI_level works (and print.lm_b works)
fitd = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          n_draws = 500,
          seed = 2025)
fita
fitd
summary(fitd)
coef(fitd)
preds0d = 
  predict(fitd)
head(preds0d)


# Make sure gls works
## Bootstrapping - parallelized
fite = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          loss = "gls",
          n_draws = 250,
          seed = 2025)
fitd
fite
summary(fite)
coef(fite)
preds0e = 
  predict(fite)
head(preds0e)

## Bootstrapping - sequential
plan(sequential)
fitf = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          loss = "gls",
          n_draws = 250,
          seed = 2025)
fite
fitf
summary(fitf)
coef(fitf)
preds0f = 
  predict(fitf)
head(preds0f)

## Large sample approx
fitg = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          loss = "gls",
          seed = 2025)
fitf
fitg
summary(fitg)
coef(fitg)
preds0f = 
  predict(fitg)
head(preds0f)



rm(list = ls())




# Test loss-likelihood approach (Poisson) ---------------------------------

# Create data

set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$y = 
  rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)


# Bootstrapping approach - sequential
plan(sequential)
## Make sure CI_level works (and print.lm_b works)
fita = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          n_draws = 250,
          seed = 2025)
fitb = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = poisson(),
          n_draws = 100,
          trials = rep(1,N),
          seed = 2025,
          CI_level = 0.8)
fita
fitb

## Make sure summary.np_lm_b works
summary(fita)
summary(fita,
        CI_level = 0.8)

## Make sure coef.lm_b works
coef(fita)
glm(y ~ x1 + x2 + x3 + offset(log(time)),
    data = test_data,
    family = poisson())

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


# Bootstrapping approach - parallelized
plan(multisession, workers = 10)
## Make sure CI_level works (and print.lm_b works)
fitc = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          n_draws = 250,
          seed = 2025)
fita
fitc
summary(fitc)
coef(fitc)
head(predict(fitc))



# large sample approach
## Make sure CI_level works (and print.lm_b works)
fitd = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025)
fite = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025,
          CI_level = 0.8)
fita
fitd
fite

## Make sure summary.np_lm_b works
summary(fitd)
summary(fitd,
        CI_level = 0.8)

## Make sure coef.lm_b works
coef(fitd)

## Make sure prediction function works
head(predict(fitd))



# Make sure gls works
## Bootstrapping - parallelized
fite = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          loss = "gls",
          n_draws = 250,
          seed = 2025)
fitd
fite
summary(fite)
coef(fite)
head(predict(fite))


## Bootstrapping - sequential
plan(sequential)
fitf = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          loss = "gls",
          n_draws = 250,
          seed = 2025)
fite
fitf
summary(fitf)
coef(fitf)
head(predict(fitf))


## Large sample approx
fitg = 
  np_lm_b(y ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          loss = "gls",
          seed = 2025)
fitf
fitg
summary(fitg)
coef(fitg)
preds0f = 
  predict(fitg)
head(preds0f)




# Check if custom loss function will work.
fitg = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          loss = function(y,mu) (y-mu)^2,
          n_draws = 250,
          seed = 2025)
fitg
summary(fitg)
coef(fitg)
head(predict(fitg))
fith = 
  np_lm_b(y ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          loss = function(y,mu) (y-mu)^2,
          seed = 2025)
fith
summary(fith)
coef(fith)
head(predict(fith))