
run <- FALSE

if (run) {
  # Clean session
  rm(list=ls())
  pacman::p_unload(all)
  gc()
  remove.packages("bayesics")
  
  # Document, including making NAMESPACE.  Safest to delete old one first.
  file.remove("C:/Users/dksewell/Documents/bayesics/NAMESPACE")
  devtools::document("C:/Users/dksewell/Documents/bayesics")
  
  # Build tarball
  devtools::build(pkg = "C:/Users/dksewell/Documents/bayesics",
                  path = "C:/Users/dksewell/Downloads",
                  vignettes = FALSE)
  
  # Install from tarball
  pacman::p_unload(bayesics)
  install.packages("C:/Users/dksewell/Downloads/bayesics_2.0.0.tar.gz",
                   repos=NULL,type='source')
  pacman::p_load(bayesics,future)
  beepr::beep(4)
  
  # Install from github
  remotes::install_github("dksewell/bayesics")
  
  
  pacman::p_load(coda,
                 dplyr,
                 extraDistr,
                 mvtnorm,
                 Matrix,
                 future,
                 future.apply,
                 ggplot2,
                 patchwork,
                 BMS,
                 cluster,
                 DFBA,
                 linreg,
                 future)
  
  beepr::beep(4)
  
}