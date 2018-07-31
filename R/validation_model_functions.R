#' MSE on a single trial
#'
#' This function computes the mse between true velocity and predicted velocity on a single trial.
#'
#' @param predicted_vel centered predicted velocity
#' @param true_vel centered true velocity
#' @return mse between true velocity and predicted velocity for a given trial
get_mse_trial <- function(predicted_vel, true_vel){

  # Edge cases
  if(is.null(predicted_vel) | !is.matrix(predicted_vel))
    stop("predicted_vel should be a matrix, T x p")
  if(is.null(true_vel) | !is.matrix(true_vel))
    stop("true_vel should be a matrix, T x p")
  if((ncol(predicted_vel) > nrow(predicted_vel)) |
     (ncol(true_vel) > nrow(true_vel))) stop("Number of rows should be greater
                                            than number of column")
  testthat::expect_equal(dim(predicted_vel), dim(true_vel))

  d <- ncol(true_vel)
  if(is.null(d)) d <- 1

  if(d == 1) mse <- mean((predicted_vel - true_vel)^2)
  if(d > 1)  mse <- mean(apply((predicted_vel - true_vel)^2, 1, sum))

  return(mse)

}


#' MSE on the whole experiment
#'
#' This function computes the mse between true velocity and predicted velocity
#' on the whole experiment. Variables need to be previously centered
#'
#' @param int_trial interval of time, can be training or testing
#' @param time_trial Time of the experiment, can be training or testing
#' @param vtrue centered true velocity matrix
#' @param vhat centered decoded velocity matrix. Can be in a list form.
#' @return mse for all trials in the training set
get_tot_mse <- function(int_trial, time_trial, vtrue, vhat){

  # Edge Cases
  if(is.list(vhat)) vhat = Reduce("rbind", vhat)
  if(is.null(vtrue) | !is.matrix(vtrue))
    stop("true_vel should be a matrix, T x p")
  if(ncol(vtrue) > nrow(vtrue)) stop("Number of rows should be greater
                                            than number of column")
  if(max(int_trial) >= length(time_trial))
     stop("int_trial should contain integer correspondent to beginning of each
          trial. the max in int_trial can't be greater than the length of the
          time_trial vector")
  testthat::expect_equal(nrow(vtrue), length(time_trial))
  testthat::expect_equal(nrow(vhat), length(time_trial))


  # Inizialize vectors, and reshape velocity if needed
  err <- vector(mode = "numeric", length = length(int_trial))

  # Computes error trial by trial
  for(j in 1:length(int_trial)){
    ix = c(int_trial[j]:ifelse(j == length(int_trial), length(time_trial),
                               int_trial[j+1] - 1))
    err[j] = get_mse_trial(vhat[ix,], vtrue[ix,])

  }
  return(err)
}

#' Decoded velocity
#'
#' This function computes the decoded velocity. Variables needs to be centered
#' @param Mole matrix needed only for computing the MLE
#' @param int_trial interval of time, can be training or testing
#' @param time_trial Time of the experiment, can be training or testing
#' @param z centered spike count matrix
#' @param v centered velocity matrix
#' @param MLE boolean for whether or not decoding is done including a prior model,
#' MLE=TRUE/FALSE exclude/includes the prior model
#' @param A state matrix
#' @param W state covariance matrix
#' @param H encoding matrix
#' @param Q encoding covariance matrix
#' @return predicted velocity in either MLE or KF scenario
vhat_fun <- function(Mole = NULL, int_trial, time_trial, z, v, mle = TRUE,
                     A = NULL, W = NULL, H = NULL, Q = NULL){

  # Edge Cases
  if(mle & is.null(Mole)) stop("Mole must be passed under MLE assumption")
  if(is.null(H) & is.null(Q) & is.null(Mole)) stop("H Q and Mole can't all be null")

  if(!mle){
    if(is.null(A)) stop("A should be passed under KF assumptions")
    if(is.null(W)) stop("W should be passed under KF assumptions")
    if(nrow(A) != ncol(A)) stop("A should be squared")
    if(nrow(W) != ncol(W)) stop("W should be squared")
    if(nrow(A) != ncol(W)) stop("A and W should have same dimensions")
    if(nrow(A) != ncol(H)) stop("H should have same number of column as the square
                                matrices A and W")
    if(nrow(Q) != ncol(Q)) stop("Q should be a square matrix")
    if(nrow(H) != ncol(Q)) stop("H should have as many rows and Qm1")

  } else{
    if(ncol(Mole) != ncol(z) | nrow(Mole) != ncol(v))
      stop("If the Mole mtrix is passed it needs to be a p x n matrix")
  }
  if(max(int_trial) >= length(time_trial))
    stop("int_trial should contain integer correspondent to beginning of each
         trial. the max in int_trial can't be greater than the length of the
         time_trial vector")
  testthat::expect_equal(nrow(v), length(time_trial))
  testthat::expect_equal(nrow(z), length(time_trial))


  vhat <- vector(mode = "list", length = length(int_trial))

  for(j in 1:length(int_trial)){
    ix = c(int_trial[j]:ifelse(j == length(int_trial), length(time_trial),
                               int_trial[j+1] - 1))

    # vhat gets computed depending on whether I'm under mle/kf settings
    vhat[[j]] = kalman_decode(z = z[ix,], MLE = mle, Mole = Mole,
                              initial_state = v[ix[1],],
                              A = A, W = W, H = H, Q = Q)
  }

  return(vhat)
}




#' Computes the current error
#'
#' This function computes current error, given a certain model
#'
#' @param vel velocity in the training set
#' @param veltest velocity in the testing set
#' @param z neural data in the training set
#' @param ztest neaural data in the testing set
#' @param A the state matrix
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param int_train vector of integer on when the trial starts in the training set
#' @param time_train vector of time of the intervals in the training set
#' @param int_test vector of integer on when the trial starts in the testing set
#' @param time_test vector of time of the intervals in the testing set
#' @param encoding whether to run the encoding model or use given A,W,H,Q
#' @param decoding whether to run the decoding model or extract only encoding pars
#' @param mle whether compute prediction for MLE or KF
#' @param decoding_train whether decoding on training set
#' @param weights set of weights for mle
#' @return a list of values with encoding parameter and decoded mse
get_current_model_info <- function(vel = velocity.Train, veltest = velocity.Test,
                                   z, ztest, A=NULL, W=NULL, H=NULL, Q=NULL,
                                   int_train = int.Train, time_train = Time.Train,
                                   int_test = int.Test, time_test = Time.Test,
                                   encoding = TRUE, decoding = TRUE, mle = TRUE,
                                   decoding_train = TRUE, weights = NULL){

  # center the variables
  v_center     <- apply(vel, 2, mean)
  spike_center <- apply(z, 2, mean)

  # encoding stage
  if(encoding){
    encoding_model <- kalman_train(x = scale(vel, center = v_center, scale = FALSE),
                                   z = scale(z, center = spike_center, scale = FALSE),
                                   fit_A = TRUE,
                                   unit_interactions = TRUE)
    A <- encoding_model$A
    W <- encoding_model$W
    H <- encoding_model$H
    Q <- encoding_model$Q
  }

  # decoding
  if(decoding){

    # by defualt decodeon testing set
    int_dec  <- int_test
    time_dec <- time_test
    zdec     <- ztest
    vdec     <- veltest

    # decide on which set decoding
    if(decoding_train)	{
      int_dec  <- int_train
      time_dec <- time_train
      zdec     <- z
      vdec     <- vel
    }

    # decoding depends on mle assumption
    Mole <- NULL
    if(mle){
      Qm1  <- solve(Q)
      Mole <- get_Mole(H = H, Qm1 = Qm1, weights = weights)
    }
    # compute prediction and relative error
    vhat <- vhat_fun(Mole = Mole, int_trial = int_dec, time_trial = time_dec,
                     z = scale(zdec, center = spike_center, scale = FALSE),
                     v = scale(vdec, center = v_center, scale = FALSE),
                     mle = mle, A = A, W = W, H = H, Q = Q)
    mse <- get_tot_mse(int_trial = int_dec, time_trial = time_dec,
                       vtrue = scale(vdec, center = v_center, scale = FALSE),
                       vhat = vhat)
  } else{
    mse = NULL
  }

  return(list(mse = mse, A = A, W = W, H = H, Q = Q))
}
