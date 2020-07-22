## Code for building large 
##
##

library(filesstrings)

##
## convert the original Rmd to a pre-configured Rmd
##

knitr::knit("vignettes/pgLM.Rmd.orig", output = "vignettes/pgLM.Rmd")
knitr::knit("vignettes/pgSTLM.Rmd.orig", output = "vignettes/pgSTLM.Rmd.orig")
# knitr::knit("vignettes/pgLM.Rmd.orig", output = "vignettes/pgLM.R")

##
## extract vignette R code
##

knitr::purl("vignettes/pgLM.Rmd.orig", output = "vignettes/pgLM.R")
knitr::purl("vignettes/pgSTLM.Rmd.orig", output = "vignettes/pgSTLM.Rmd.orig")
# knitr::purl("vignettes/pgLM.Rmd.orig", output = "vignettes/pgLM.R")

##
## move images to the vignette folder
##

## list the images created by the vignette and move to the vignettes folder
pg_lm_images <- list.files(".")[grep("pglm-", list.files("."))]
file.move(pg_lm_images, destinations = "./vignettes/", overwrite = TRUE)

##
## build vignettes
##

## note: you can have build errors if the file names have spaces: https://stackoverflow.com/questions/27338970/vignette-creation-on-package-build-fails-with-the-error-failed-to-locate-the-w
devtools::build_vignettes()

##
## build the documentation site
##

pkgdown::build_site()
