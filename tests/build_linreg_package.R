rm(list=ls())
pacman::p_unload(all)

devtools::document("C:/Users/dksewell/Documents/linreg")

remove.packages("linreg")

devtools::build(pkg = "C:/Users/dksewell/Documents/linreg",
                path = "C:/Users/dksewell/Downloads",
                vignettes = FALSE)

pacman::p_unload(linreg)
install.packages("C:/Users/dksewell/Downloads/linreg_2.0.0.tar.gz",
                 repos=NULL,type='source')

remotes::install_github("dksewell/linreg")
