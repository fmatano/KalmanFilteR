#' Get tuning curves r-square
#'
#' This function computes r squared for each tuning-curve model
#' @param z matrix of spike counts
#' @param cov matrix of covariates
get_r2 <-  function(z, cov){

  # Fit the model
  fit_enc <- lm(z ~ cov)

  # Get the R2
  r2 <- sapply(summary(fit_enc), function(x) x$r.squared)

  # Rename the coloumns
  names(r2) <- 1:ncol(z)

  return(r2)

}




#' Get tuning curves sigma-square
#'
#' This function computes sigma square for each tuning-curve model
#' @param z matrix of spike counts
#' @param cov matrix of covariates
get_sigma2 <- function(z, cov){

  # Fit the model
  fit_enc <- lm(z ~ cov)

  # Get the R2
  sigma2 <- (sapply(summary(fit_enc), function(x) x$sigma))^2

  # Rename the coloumns
  names(sigma2) <- 1:ncol(z)


  return(sigma2)

}


#' Get tuning curves bias
#'
#' This function computes the bias for each tuning-curve model
#' @param z matrix of spike counts
#' @param cov matrix of covariates
get_bias <- function(z, cov){

  require(gam)
  mse_lin <- c()
  mse_spl <- c()

  # Fit the model
  for(neur in 1:ncol(z)){
    #cat("\n~~~~~~~~~  Equation [", neur, "]  ~~~~~~~~~\n")
    lin_mod <- gam(z[,neur] ~ cov)
    spl_mod <- gam(z[,neur] ~ s(cov[,1], 2) + s(cov[,2], 2) + s(cov[,3], 2))

    mse_lin[neur] <- mean(residuals(lin_mod)^2)
    mse_spl[neur] <- mean(residuals(spl_mod)^2)
  }

  bias <- (mse_lin - mse_spl)^2

  # Rename the coloumns
  names(bias) <- 1:ncol(z)
  return(bias)

}



#' Get tuning curves optimal boxcox transformation
#'
#' This function computes the optimal boxcox transformation for each tuning-curve model
#' @param z matrix of spike counts
#' @param cov matrix of covariates
#' @param shift shift for computing boxcox transformation for avoiding dividing by zero
get_optimal_boxcox <- function(z, cov, shift=0.0001){

  library(MASS)
  lambda <- c()
  for(i in 1:ncol(z)){

    y <- (z[,i] + shift)
    lm_mod_shift <- lm(Y~., data = data.frame(Y=y, X=cov))
    bx_mod       <- boxcox(lm_mod_shift, plotit=FALSE,
                           data=data.frame(Y=y, X=cov))
    lambda[i]    <- bx_mod$x[which.max(bx_mod$y)]   # <--- optimal lambda

  } # <---- End by neuron

  return(lambda)
}



#' Build the equation info dataframe
#'
#' This function builds the equation info dataframe, with info for every equation about electrode, lag, transformation, waveforms characteristics and r-square score
#' @param r2 r-square score
#' @param z_colnames names of the spike count matrix columns
#' @param n_units number of units in the data set
#' @param n_lags number of lags in the data set
#' @param n_transf number of spike count transformation in the data set
#' @param n_wf number of waveforms in the data set
#' @param n_feat number of waveform features in the data set
build_df_equations <- function(r2, z_colnames, n_units, n_lags = 13,
                               n_transf = 3, n_wf = 0, n_feat = 4,
                               n_transf_wf = 1){

  if(is.null(z_colnames)) stop("z_colnames cannot be null")

  # producing the equation key, in terms of correspondent electrode
  var_key  <- rep(1:n_units, n_lags*n_transf)
  transf   <- rep(1:n_transf, each=n_lags*n_units)
  lag      <- rep(rep(1:n_lags, each=n_units), n_transf)
  wf <- wf_key <- wf_transf <- rep(0, n_lags*n_units*n_transf)

  # producing the equation key, if waveforms are included
  if(n_wf > 0){
    var_key <- c(var_key, rep(1:n_units, n_wf*n_feat*n_lags*n_transf_wf))
    transf  <- c(transf, rep(0, n_wf*n_feat*n_lags*n_units*n_transf_wf))
    lag     <- c(lag, rep(rep(1:n_lags, each= n_units*n_feat), n_wf*n_transf_wf))
    wf      <- c(wf, rep(1:n_wf, each = n_units*n_feat*n_lags*n_transf_wf))
    wf_key  <- c(wf_key, rep(rep(1:n_feat, each = n_units), n_lags*n_wf*n_transf_wf))
    wf_transf <- c(wf_transf, rep(rep(1:n_transf_wf, each = n_units*n_feat*n_lags),
                                  n_wf))
  }

  n_eqs <- length(var_key)

  eq_info <- cbind(1:n_eqs, var_key, transf, lag, wf, wf_key, wf_transf)
  colnames(eq_info) <- c('eq', 'key', 'transf', 'lag', 'wf', 'wf_key', 'wf_transf')

  return(eq_info)
}




#' Extracting equation info
#'
#' This function extracts equation info, and saves it in save_name
#' @param experiment_fname experiment file name to be loaded
#' @param z spike count training set
#' @param n_units number of units in the data set
#' @param n_lags number of lags in the data set
#' @param n_wf number of waveforms in the data set
#' @param n_feat number of waveform features in the data set
#' @param n_transf number of spike count transformation in the data set
#' @param save_name name of the file where saving the equations info
#' @param vel velocity training set
get_equations_info <- function(experiment_fname, z, n_units, n_lags, n_wf, n_feat,
                               n_transf, save_name, vel, n_transf_wf = 1,
                               compute_rsq = TRUE){

  load(experiment_fname)

  # check dimensions, not all columns have ACE if they are zero
  rsq  <- rep(NA, n_units*n_lags*n_transf + n_wf*n_units*n_feat*n_lags*n_transf_wf)
  # bias <- rep(NA, n_units*n_lags*n_transf + n_wf*n_units*n_feat*n_lags)

  if(ncol(z) != length(rsq))
    stop("There is an error in computing the spike count matrix")

  colnames(z) <- 1:ncol(z)
  spike_center <- apply(z, 2, mean)
  v_center <- apply(vel, 2, mean)

  # computing the R2
  if(compute_rsq)
    rsq  <- get_r2(cov = scale(vel, center = v_center, scale = FALSE),
                  z = scale(z, center = spike_center, scale = FALSE))

  cat("r squared", ifelse(!compute_rsq , "not", ""), "computed \n")

  # bias <- get_bias(cov = scale(velocity.Train, center = v.center, scale = FALSE),
  #                 z = scale(Z, center = spike.center, scale = FALSE))

  eqs_info <- build_df_equations(rsq, 1:ncol(z),
                                  n_units = n_units, n_lags = n_lags,
                                  n_transf = n_transf, n_wf = n_wf,
                                  n_feat = n_feat, n_transf_wf = n_transf_wf)

  cat("equation info extracted! \n")

  save(rsq, eqs_info, file=save_name)
}







