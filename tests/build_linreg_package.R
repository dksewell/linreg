library(devtools)

pacman::p_unload(all)
rm(list=ls())

remove.packages("linreg")

devtools::build(pkg = "C:/Users/dksewell/Documents/linreg",
                path = "C:/Users/dksewell/Downloads",
                vignettes = FALSE)

pacman::p_unload(all)
install.packages("C:/Users/dksewell/Downloads/linreg_2.0.0.tar.gz",#'C:/Users/dksewell/Documents/linear_subspace_shrinkage/code/SUBSET_0.0.0.9000.tar.gz',
                 repos=NULL,type='source')

remotes::install_github("dksewell/linreg")
