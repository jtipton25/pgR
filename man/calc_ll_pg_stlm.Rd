% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_ll.R
\name{calc_ll_pg_stlm}
\alias{calc_ll_pg_stlm}
\title{Log-likelihood for pg_stlm() model for use in model selection}
\usage{
calc_ll_pg_stlm(Y, X, out)
}
\arguments{
\item{Y}{is a \eqn{N \times J \times T}{N x J x T} array of compositional count data.}

\item{X}{is a \eqn{N \times p}{n_sites x p} matrix of climate variables.}

\item{out}{is a list of MCMC outputs from \code{pg_stlm()}}
}
\description{
this function generates the log-linkelihood for data fit using the pg_splm()
function. The log-likelihood can be used for model selection and evaluation.
}
