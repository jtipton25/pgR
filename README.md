# pgR

## Package Installation
To install the package, you can use 

```r
devtools::install_github("jtipton25/pgR")
```

or you can clone the project and install it locally using RStudio. You will need a personal access token -- see [here](https://happygitwithr.com/github-pat.html) for how to set this up. It's pretty simple, just use 

```r
usethis::browse_github_pat()
``` 

and copy the token into your .Renviron file as 

```r
GITHUB_PAT=XXXXX
```

where the XXXX is the PAT copied from github. If you are on a Mac, make sure you have an openMP supported compiler -- see [here](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) for instructions on how to get this setup. For Windows, download and install RTools. For R (>=4.0.0) follow the instructions [here](https://cran.r-project.org/bin/windows/Rtools/) -- older versions of R follow the instructions [here](https://cran.r-project.org/bin/windows/Rtools/history.html)
