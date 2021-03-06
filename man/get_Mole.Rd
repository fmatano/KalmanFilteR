% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kalman_filter_functions.R
\name{get_Mole}
\alias{get_Mole}
\title{Function to compute the Mole matrix}
\usage{
get_Mole(H, Qm1, weights = NULL)
}
\arguments{
\item{H}{encoding matrix}

\item{Qm1}{inverse of the encoding covariance matrix}

\item{weights}{possible weights to weight the covariance matrix}
}
\value{
the Mole matrix, to be multiplied times the testing neural data for making
predictions
}
\description{
This function computes the Mole matrix given the encoring parameter H and the inverse of the covariance matrix Q.
}
