#' Transforming data
#'
#' This function transformes data using sqrt and ACE transformation
#' @param enc_mat the encoding matrix of data
#' @param dec_mat the decoding matrix of data
#' @param vel_train matrix of velocity from training set
#' @param vel_test matrix of velocity from testing set
#' @param dead_channels non active channels
get_transformed_scdata <- function(enc_mat, dec_mat, vel_train = velocity.Train,
																 vel_test = velocity.Test, dead_channels=NULL,
																 categorical_transf = TRUE){

	# Remove dead channels
	if(length(dead_channels) > 0){
		enc_mat <- sapply(enc_mat, function (x) x[,-c(dead_channels)], simplify = FALSE)
		dec_mat <- sapply(dec_mat, function (x) x[,-c(dead_channels)], simplify = FALSE)
	}

  cat('SC matrix transformation: \n')
	z <- do.call(cbind, enc_mat)
	ztest <- do.call(cbind, dec_mat)

	# Ace transform of each lag:
	# nonzero_col <- which(colSums(z)!=0)
	categ_selected <- ifelse(categorical_transf, 0, 1)
	ACE_train <- ACE_bycol(z, vel_train, categ = categ_selected)
	ACE_test  <- ACE_bycol(ztest, vel_test, categ = categ_selected)
	cat("- ace succesful! \n")

	# Transformed data
	z <- cbind(z, ACE_train, sqrt(z))
	ztest <- cbind(ztest, ACE_test, sqrt(ztest))
	cat("- sqrt succesful! \n")

	return(list(Z=z, Ztest=ztest))
}


#' Applying ACE transformation by column
#'
#' This function applyies ACE transformation by column for either SC or WF
#' @param mat sc or wf matrix
#' @param kine kinematic matrix, either training or testing, according to mat
#' @param categ indicates whether categorica variables are present,
#' 0 for spike counts, NULL for waveforms
#' @return ace transformed data
ACE_bycol <- function(mat, kine, categ){

  if(nrow(kine) != nrow(mat)) stop("You are using incompatible neural
                                   and kinematic data")

  if(!is.null(categ) && categ == 1) categ <- NULL
  ACE <-  NULL
  for(j in 1:ncol(mat))
    ACE <- cbind(ACE, ace_transform(mat[,j], kine, mat[,j], categ))

  return(ACE)
}

#' Building WF data eliminating dead channel
#'
#' This function builds the waveform matrix, that by default should be a list
#' with 13 elements (lag), each element is a list of 4 matrixes (features),
#' each matrix is T x n_units
#' @param wf_mat waveform matrix
#' @param dead_channels channels to eliminate
#' @return a list of a list, as in input, but without non-active channels
build_wf_mat <- function(wf_mat, dead_channels)
  wf_mat %>%
    lapply(function(x) lapply(x, function(y) y[,-c(dead_channels)]))



#' Transform WF matrix
#'
#' This function transform the waveform matrix, that is now a matrix T x
#' n_units n_lags n_feat in matrix where each column has been transformed
#' according to the FUN function passed
#' @param wf_mat waveform matrix
#' @param FUN function to use to modify the data
#' @param kine kinematic data, needed only for ace transformation
#' @return a matrix with the same dimension of the original one, but transformed
#' data
transform_wf_mat <- function(wf_mat, FUN, kine = NULL){

  if(FUN == "ace_transform"){
    wf_transformed <- wf_mat %>%
      lapply(function(x) lapply(x, ACE_bycol, kine, categ = NULL))
  } else{
    wf_transformed <- wf_mat %>%
      lapply(function(x) lapply(x, function(x) match.fun(FUN)(x)))
  }

  if(!all.equal(dim(wf_transformed[[1]][[1]]), dim(wf_mat[[1]][[1]])))
    stop("Output matrix doesn't match input matrix dimension")
  return(wf_transformed)

}


#' Transforming waveform data
#'
#' This function transformes data using sqrt and ACE transformation
#' @param enc_wfmat1 the encoding matrix of data, first moment
#' @param enc_wfmat2 the encoding matrix of data, first moment
#' @param dec_wfmat1 the decoding matrix of data, first moment
#' @param dec_wfmat2 the decoding matrix of data, first moment
#' @param vel_train matrix of velocity from training set
#' @param vel_test matrix of velocity from testing set
#' @param dead_channels non active channels
#' @return a list with WF matrix for encoding and decoding
get_transformed_wfdata <- function(enc_wfmat1, enc_wfmat2,
                                   dec_wfmat1, dec_wfmat2,
                                   vel_train = velocity.Train,
                                   vel_test = velocity.Test, dead_channels=NULL){

  if(length(dead_channels) > 0){
    WF_Enc <- build_wf_mat(enc_wfmat1, dead_channels)
    WF_Enc_sq <- build_wf_mat(enc_wfmat2, dead_channels)
    WF_Dec <- build_wf_mat(dec_wfmat1, dead_channels)
    WF_Dec_sq <- build_wf_mat(dec_wfmat2, dead_channels)
  }

  # transform also WF data
  cat('WF matrix transformation: \n')

  WF_Enc_sqrt <- transform_wf_mat(WF_Enc, "sqrt")
  cat("-")
  WF_Enc_sq_sqrt <- transform_wf_mat(WF_Enc_sq, "sqrt")
  cat("-")
  WF_Dec_sqrt <- transform_wf_mat(WF_Dec, "sqrt")
  cat("-")
  WF_Dec_sq_sqrt <- transform_wf_mat(WF_Dec_sq, "sqrt")
  cat("- sqrt succesful! \n")

  WF_Enc_ace <- transform_wf_mat(WF_Enc, "ace_transform", vel_train)
  cat("-")
  WF_Enc_sq_ace <- transform_wf_mat(WF_Enc_sq, "ace_transform", vel_train)
  cat("-")
  WF_Dec_ace <- transform_wf_mat(WF_Dec, "ace_transform", vel_test)
  cat("-")
  WF_Dec_sq_ace <- transform_wf_mat(WF_Dec_sq, "ace_transform", vel_test)
  cat("- ace succesful! \n")

  WF_Enc_tot <- cbind(do.call(cbind, do.call(cbind, WF_Enc)),
                      do.call(cbind, do.call(cbind, WF_Enc_ace)),
                              do.call(cbind, do.call(cbind, WF_Enc_sqrt)),
                      do.call(cbind, do.call(cbind, WF_Enc_sq)),
                      do.call(cbind, do.call(cbind, WF_Enc_sq_ace)),
                      do.call(cbind, do.call(cbind, WF_Enc_sq_sqrt)))

  WF_Dec_tot <- cbind(do.call(cbind, do.call(cbind, WF_Dec)),
                      do.call(cbind, do.call(cbind, WF_Dec_ace)),
                      do.call(cbind, do.call(cbind, WF_Dec_sqrt)),
                      do.call(cbind, do.call(cbind, WF_Dec_sq)),
                      do.call(cbind, do.call(cbind, WF_Dec_sq_ace)),
                      do.call(cbind, do.call(cbind, WF_Dec_sq_sqrt)))

  return(list(WF_Enc_tot = WF_Enc_tot, WF_Dec_tot = WF_Dec_tot))
}





#' Select portion of movement
#'
#' This function selects a portion of movement
#' @param int vector of intervals
#' @param how_long how many steps to select, from 0 to how_long
select_movement <- function(int, how_long = 10){

	# test the length chosen
	if(how_long > min(diff(int))) stop("The length of the trial selected
																		 is longer than some of the trials")

	# select start and end of the reduced trial
	start <- int
	end   <- int + how_long

	# compute the final interval for all the trials
	interval_start_end <- cbind(start, end)
	interval_selected  <- c(apply(interval_start_end, 1, function(x)
														 	eval(parse(text=paste(x[1], ":", x[2])))))

	return(interval_selected)
}



