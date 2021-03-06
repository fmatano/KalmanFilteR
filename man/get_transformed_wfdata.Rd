% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setting_data.R
\name{get_transformed_wfdata}
\alias{get_transformed_wfdata}
\title{Transforming waveform data}
\usage{
get_transformed_wfdata(enc_wfmat1, enc_wfmat2, dec_wfmat1, dec_wfmat2,
  vel_train = velocity.Train, vel_test = velocity.Test,
  dead_channels = NULL)
}
\arguments{
\item{enc_wfmat1}{the encoding matrix of data, first moment}

\item{enc_wfmat2}{the encoding matrix of data, first moment}

\item{dec_wfmat1}{the decoding matrix of data, first moment}

\item{dec_wfmat2}{the decoding matrix of data, first moment}

\item{vel_train}{matrix of velocity from training set}

\item{vel_test}{matrix of velocity from testing set}

\item{dead_channels}{non active channels}
}
\value{
a list with WF matrix for encoding and decoding
}
\description{
This function transformes data using sqrt and ACE transformation
}
