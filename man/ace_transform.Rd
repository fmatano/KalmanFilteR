% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ace_transformation_function.R
\name{ace_transform}
\alias{ace_transform}
\title{ACE transformation}
\usage{
ace_transform(sc_train, kine_train, sc_new, categ = 0)
}
\arguments{
\item{sc_train}{vector of spike count in the training set}

\item{kine_train}{design matrix of the kinematic linear predictor in the training set}

\item{sc_new}{vector of spike counts whose tranform you need}

\item{categ}{0 for spike counts, NULL for waveforms}
}
\value{
the transformed data
}
\description{
This function computes the ACE transformation for spike counts or any enter in SCtrain argument
}
