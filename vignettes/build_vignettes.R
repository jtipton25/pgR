## Code for building long running vignettes
##
## From: https://ropensci.org/blog/2019/12/08/precompute-vignettes/
## From: https://blog.r-hub.io/2020/06/03/vignettes/#how-to-include-a-compute-intensive--authentication-dependent-vignette

library(filesstrings)
library(xfun)

##
## convert the original Rmd to a pre-configured Rmd
##

knitr::knit("vignettes/pg_lm.Rmd.orig", output = "vignettes/pg_lm.Rmd", envir = new.env())
## change the file path in the .Rmd
gsub_file("vignettes/pg_lm.Rmd", "figure/", "")

##
## move images to the vignette folder and delete cached files
##

## list the images created by the vignette and move to the vignettes folder
# mra_images <- list.files("figure/")[grep(".svg", list.files("figure/"))]
# file.move(paste0("figure/", mra_images), destinations = "./vignettes/", overwrite = TRUE)
pg_lm_images <- list.files("figure/")[grep(".png", list.files("figure/"))]
file.move(paste0("figure/", pg_lm_images), destinations = "./vignettes/", overwrite = TRUE)

if(length(dir("figure/", all.files = TRUE)) ==0)
    file.remove("./figure/*")

file.remove("./figure/")
cache_files <- list.files("cache/")
file.remove(paste0("cache/", cache_files))
file.remove("./cache")



# knitr::knit("vignettes/mra-simulation.Rmd.orig", output = "vignettes/mra-simulation.Rmd")



##
## extract vignette R code
##

knitr::purl("vignettes/pg_lm.Rmd.orig", output = "vignettes/pg_lm.R")
# knitr::purl("vignettes/mra-simulation.Rmd.orig", output = "vignettes/mra-simulation.R")



##
## build vignettes
##

## note: you can have build errors if the file names have spaces: https://stackoverflow.com/questions/27338970/vignette-creation-on-package-build-fails-with-the-error-failed-to-locate-the-w
devtools::build_vignettes()

##
## build the documentation site
##

# pkgdown::build_site()
