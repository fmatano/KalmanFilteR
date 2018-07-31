#' ACE transformation
#'
#' This function computes the ACE transformation for spike counts or any enter in SCtrain argument
#' @param sc_train vector of spike count in the training set
#' @param kine_train design matrix of the kinematic linear predictor in the training set
#' @param sc_new vector of spike counts whose tranform you need
#' @param categ 0 for spike counts, NULL for waveforms
#' @return the transformed data
ace_transform <- function(sc_train, kine_train, sc_new, categ = 0){

  ace_out <- acepack::ace(kine_train, sc_train, cat = categ, lin = c(1))
  return( Hmisc::approxExtrap(sc_train, ace_out$ty, xout=sc_new)$y )
}

