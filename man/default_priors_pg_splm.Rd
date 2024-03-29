% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default-priors.R
\name{default_priors_pg_splm}
\alias{default_priors_pg_splm}
\title{Initialized default values of priors for \code{pg_splm()}}
\usage{
default_priors_pg_splm(Y, X, corr_fun = "exponential")
}
\arguments{
\item{Y}{is a \eqn{n \times d}{n x d} matrix of compositional count data.}

\item{X}{is a \eqn{n \times p}{n x p} matrix of climate variables.}

\item{corr_fun}{is a character that denotes the correlation function form. Current options include "matern" and "exponential".}
}
\description{
A function for setting up the default \code{pg_splm()} priors
}
