% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_kappa}
\alias{calc_kappa}
\title{Calculate the intermediate Polya-gamma value kappa}
\usage{
calc_kappa(Y, Mi)
}
\arguments{
\item{Y}{is a \eqn{n \times J}{n x J} matrix of compositional count data.}

\item{Mi}{is a \eqn{n \times J-1}{n x J-1} matrix of transformed observations for Polya-gamma that is the output of the \code{calc_Mi()} function}
}
\description{
Calculate the intermediate Polya-gamma value kappa
}
