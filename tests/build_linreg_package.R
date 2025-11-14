# Clean session
rm(list=ls())
pacman::p_unload(all)
gc()
remove.packages("linreg")

# Document, including making NAMESPACE.  Safest to manually delete old one first.
devtools::document("C:/Users/dksewell/Documents/linreg")

# Build tarball
devtools::build(pkg = "C:/Users/dksewell/Documents/linreg",
                path = "C:/Users/dksewell/Downloads",
                vignettes = FALSE)

# Install from tarball
pacman::p_unload(linreg)
install.packages("C:/Users/dksewell/Downloads/linreg_2.0.0.tar.gz",
                 repos=NULL,type='source')
pacman::p_load(linreg,future)
beepr::beep(4)

# Install from github
remotes::install_github("dksewell/linreg")
