if(FALSE){
  library(devtools)
  
  document("C:/Users/dksewell/Documents/linreg/")
  # For some reason, roxygen2 is not including:
  # S3method(linreg::coef,aov_b)
  # S3method(linreg::summary,aov_b)
  # in the namespace.  Add manually.
  
  
  build(pkg = "C:/Users/dksewell/Documents/linreg/",
        path = "C:/Users/dksewell/Documents/")
  
  pacman::p_unload(all)
  
  install.packages("C:/Users/dksewell/Documents/linreg_0.0.0.9000.tar.gz",
                   type = "source",
                   repos = NULL)

}


# Test package ------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(linreg)

sim_data = 
  tibble(x = sample(letters[1:3],200,TRUE)) %>% 
  mutate(y = -1 + case_match(x,
                             "a" ~ -2,
                             "b" ~ 1.0,
                             "c" ~ 0.25) +
           rnorm(200))

aovfits = list()
aovfits$c = 
  aov_b(y~x,sim_data)
aovfits$i = 
  aov_b(y~x,sim_data,
        improper = TRUE)

aovfits$c
aovfits$i
coef(aovfits$c)
coef(aovfits$i)
summary(aovfits$c)
summary(aovfits$i)

sapply(aovfits,AIC)
sapply(aovfits,BIC)
sapply(aovfits,DIC)
sapply(aovfits,WAIC)

aovfits$c = 
  aov_b(y~x,sim_data,
        heteroscedastic = FALSE)
aovfits$i = 
  aov_b(y~x,sim_data,
        heteroscedastic = FALSE,
        improper = TRUE)

aovfits$c
aovfits$i
coef(aovfits$c)
coef(aovfits$i)
summary(aovfits$c)
summary(aovfits$i)

sapply(aovfits,AIC)
sapply(aovfits,BIC)
sapply(aovfits,DIC)
sapply(aovfits,WAIC)


sim_data = 
  tibble(x = rnorm(200,mean = 10)) %>% 
  mutate(y = -1 + 0.75 * x +
           rnorm(200))
fits = list()
fits$c = 
  lm_b(y~x,sim_data,"conjugate")
fits$z = 
  lm_b(y~x,sim_data,"zellner")
fits$i = 
  lm_b(y~x,sim_data,"improper")

fits$c
fits$z
fits$i
coef(fits$c)
coef(fits$z)
coef(fits$i)
summary(fits$c)
summary(fits$z)
summary(fits$i)

sapply(fits,AIC)
sapply(fits,BIC)
sapply(fits,DIC)
sapply(fits,WAIC)
