% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{solve_ensemble}
\alias{solve_ensemble}
\title{The adaptive weighted ensemble-based learning method to integrate the multiple binary spots similarity matrix}
\usage{
solve_ensemble(
  Results.clustering,
  lambda = NULL,
  prob.quantile = 0.5,
  niter = 100,
  epsilon = 1e-05
)
}
\arguments{
\item{Results.clustering}{a list contains all the results of individual similarity matrix. The elements of list is a matrix, spots * spots.}

\item{lambda}{hyper-parameter constrain the weight of individual methods for ensemble. If the parameter is set to NULL, then, we will adopt the value in our algorithm.}

\item{prob.quantile}{numeric of probabilities with values in [0,1]. Default setting is 0.5.}

\item{niter}{a positive integer represents the maximum number of updating algorithm. Default setting is 100.}

\item{epsilon}{a parameter represents the stop criterion.}
}
\value{
a list contains a matrix of the ensemble similarity of spots and a vector of the weight assigned to base results.
}
\description{
The adaptive weighted ensemble-based learning method to integrate the multiple binary spots similarity matrix
}
