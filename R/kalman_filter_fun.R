#' Function to encode parameters of a Kalman filter
#'
#' This function encodes parameters for a state space model with state
#' parameters A and W and observation parameters H and Q
#' @param x state (usually T x 3)
#' @param z spike counts (usually T x 96)
#' @param interruptions NA unless the data is segmented (i.e. time line
#' is interrupted) A printed string
#' @param  fit_A fit a transition matrix. fit_A = FALSE uses an identity state
#' transition matrix
#' @param noise_units forces flat tuning curves for a subset if 1
#' @param  unit_interactions if TRUE the state transition var-covar matrix is
#'  diagonal
#' @param  single_neuron if = 0 all are considered, otherwise the numb of
#' specific neuron has to be specified
#' @param	lambda parameter for ridge
#' @return a list of estimated parameters for the state space model
kalman_train <- function (x, z, interruptions = NA, fit_A = FALSE,
  noise_units = rep(0, ncol(z)), unit_interactions = TRUE, single_neuron = 0,
  lambda = 0){

  # Edge cases
  if(is.null(x) | !is.matrix(x)) stop("x has to be a matrix")
  if(is.null(z) | !is.matrix(z)) stop("z has to be a matrix")
  if(nrow(x) != nrow(z)) stop("x and z need to have the same number of rows")
  if(length(noise_units) != ncol(z)) stop("Noise units should have same
                                          length as the number of columns in z")

  TT <- nrow(x)

  # 1. Fit evolution parameters A, W ~~~~~~~~~~~~~~~~~~~~~~~~
  n_obs <- TT - length(interruptions)
  # number of observations to be used in estimating the state evolution parameters.
  # One observation is discarded for each segment of data.

  # the index of points for which the previous state is observed.

  if(all(is.na(interruptions))){
    pnts <- c(2:TT)
  }else{
    pnts <- c(1:TT)[-interruptions]
  }
  #	cat("length(pnts) = ", length(pnts),"\n")
  if(fit_A){
    A1 <- A2 <- w1 <- w2 <- 0
    for(k in pnts){
      A1 <- A1 + outer(x[k,], x[k-1,])
      A2 <- A2 + outer(x[k-1, ], x[k-1, ])
      w1 <- w1 + outer(x[k,], x[k,])
      w2 <- w2 + outer(x[k-1,], x[k,])
    }

    A <- A1%*%solve(A2)

    w <- w1 - A%*%w2

  }else{
    A <- diag(ncol(x))
    w <- 0
    for(k in pnts){
      w <- w + outer(x[k,]-x[k-1,],x[k,]-x[k-1,])
    }
  }

  W <- w/n_obs


  # 2. Fit encoding parameters H, Q  ~~~~~~~~~~~~~~~~~~~~~~~~

  Ha <- Hb <- Qa <- Qb <- 0

  for(k in 1:TT){
    Ha <- Ha + outer(z[k,], x[k,])
    Hb <- Hb + outer(x[k,], x[k,])
    Qa <- Qa + outer(z[k,], z[k,])
    Qb <- Qb + outer(x[k,], z[k,])
  }

  H <- Ha %*% solve(Hb + diag(lambda, nrow(Hb)))
  H[which(noise_units == 1),] <- 0


  Q <- (Qa - H %*% Qb)/TT + lambda * H %*% t(H)

  if(!unit_interactions){Q <- diag(diag(Q))}

  return(list(A = A, W = W, H = H, Q = Q))
}



#' Function to compute the Mole matrix
#'
#' This function computes the Mole matrix given the encoring parameter H and the inverse of the covariance matrix Q.
#' @param H encoding matrix
#' @param Qm1 inverse of the encoding covariance matrix
#' @param weights possible weights to weight the covariance matrix
#' @return the Mole matrix, to be multiplied times the testing neural data for making
#' predictions
get_Mole <- function(H, Qm1, weights = NULL){

  # Edge cases
  if(is.null(H) | is.null(Qm1)) stop("Neither H nor Qm1 can't be null")
  if(nrow(Qm1) != ncol(Qm1)) stop("Qm1 should be a square matrix")
  if(nrow(H) != ncol(Qm1)) stop("H should have as many rows and Qm1")
  if(!is.null(weights) & length(weights) != nrow(Qm1))
    stop("Vector of weigths is either null or should be same length as
         row and columns of Qm1")

  Lambda <- Qm1

  # if weights, build the right Lambda
  if(!is.null(weights)) {
    weight_mat <- diag(sqrt(weights))
    Lambda <- weight_mat %*% Qm1 %*% weight_mat
  }

  tH <- t(H)
  M <- solve(tH %*% Lambda %*% H) %*% tH %*% Lambda

  return(M)

}


#' Function to decode the Kalman Filter trajectory
#'
#' This function decodes trajectory for a state space model with state parameters A and W and observation parameters H and Q
#' @param A state matrix
#' @param W state covariance matrix
#' @param H encoding matrix
#' @param Q encoding covariance matrix
#' @param weights possible weights to weight the covariance matrix
#' @param z centered spike count matrix
#' @param initial_state centered initial velocity state
#' @param MLE boolean for whether or not decoding is done including a prior model, MLE=TRUE/FALSE exclude/includes the prior model
#' @param weights weights to rescale the covariance matrix
#' @param Mole matrix needed only for computing the MLE
#' @return  decoded trajectory for a given trial if not MLE, for the whole testing
#' set if MLE
kalman_decode <- function(A = NULL, W = NULL, H = NULL, Q = NULL, z,
                          initial_state = rep(0, ncol(H)), MLE = TRUE,
                          weights = NULL, Mole = NULL){
  # Edge cases
  if(MLE & is.null(Mole)) stop("Mole must be passed under MLE assumption")
  if(is.null(H) & is.null(Q) & is.null(Mole)) stop("H Q and Mole can't all be null")

  if(!MLE){
    if(is.null(A)) stop("A should be passed under KF assumptions")
    if(is.null(W)) stop("W should be passed under KF assumptions")
    if(nrow(A) != ncol(A)) stop("A should be squared")
    if(nrow(W) != ncol(W)) stop("W should be squared")
    if(nrow(A) != ncol(W)) stop("A and W should have same dimensions")
    if(nrow(A) != ncol(H)) stop("H should have same number of column as the square
                              matrices A and W")
    if(nrow(Q) != ncol(Q)) stop("Q should be a square matrix")
    if(nrow(H) != ncol(Q)) stop("H should have as many rows and Qm1")
    if(!is.null(weights) & length(weights) != nrow(Q))
      stop("Vector of weigths is either null or should be same length as
           row and columns of Q")
    if(length(initial_state) != ncol(H)) stop("Intial state should have same
                                              length as the number of columns in H")

  } else{
    if(!is.null(weights) & length(weights) != ncol(Mole))
      stop("Vector of weigths is either null or should be same length as
           row and columns of Q")
    if(length(initial_state) != nrow(Mole)) stop("Intial state should have same
                                            length as the number of columns in H")
    if(ncol(Mole) != ncol(z) | nrow(Mole) != length(initial_state))
      stop("If the Mole mtrix is passed it needs to be a p x n matrix")
  }


  if(!MLE){

    TT <- nrow(z)
    d <- ncol(A)

    x_hat   <- matrix(NA, nrow = TT, ncol = d)
    x_hat_p <- matrix(NA, nrow = TT, ncol = d)

    x_hat[1,] <- initial_state
    Kt <- Pt_m1 <- NULL
    Pt_m1 <- matrix(0, nrow = d, ncol = d)

    for(t in 2:TT){

      # 1. predict
      x_hat_p[t,] <- c(A %*% x_hat[(t-1),])
      P_p <- A %*% Pt_m1 %*% t(A) + W

      # 2. update
      Mt <- H %*% P_p %*% t(H) + Q
      # cat('symmetric matrix:', isSymmetric(Mt), '\n')
      Kt <- P_p %*% t(H) %*% solve(Mt)
      x_hat[t,] <- x_hat_p[t,] + c( Kt %*% (z[t,] - c(H %*% x_hat_p[t,])) )
      Pt_m1 <- (diag(1, ncol = d, nrow = d) - Kt %*% H) %*% P_p

    }
  }

  else {

    x_hat <- t(Mole %*% t(z))
    x_hat[1, ] <- initial_state

  }

  return(pred_state = x_hat)
}

#' Decoding trajectory with short cut formula for leaving one equation out
#'
#' Function to decode the trajectory using a short cut formula for leaving one equation out.
#'
#' @param A the state matrixe
#' @param W the state covariance matrix
#' @param H the encoding matrix
#' @param Q the encoding covariance matrix
#' @param weights possible weights to weight the covariance matrix
#' @param z the centered spike count matrix
#' @param length_current length of current model, needed only for add_one_in
#' @param initial_state is the centered initial velocity state
#' @return  a list of decoded trajectory for a given trial obtained by dropping
#' one equation at the time. Every element of the list corresponds to an equation
#' dropped
kalman_decode_loo <- function(A = NULL, W = NULL, H = NULL, Q = NULL, z,
                             length_current = NULL,
                             initial_state = rep(0, ncol(H))){


  # Edge cases
  if(is.null(W) | is.null(A) | is.null(H) | is.null(Q))
    stop("Neither A nor W nor H nor Q can ever be null")
  if(nrow(A) != ncol(A)) stop("A should be squared")
  if(nrow(W) != ncol(W)) stop("W should be squared")
  if(nrow(A) != ncol(W)) stop("A and W should have same dimensions")
  if(nrow(A) != ncol(H)) stop("H should have same number of column as the square
                              matrices A and W")
  if(is.null(H) | is.null(Q)) stop("Neither H nor Q can ever be null")
  if(nrow(Q) != ncol(Q)) stop("Q should be a square matrix")
  if(nrow(H) != ncol(Q)) stop("H should have as many rows and Qm1")
  if(is.null(z) | !is.matrix(z)) stop("z has to be a matrix")
  if(length(initial_state) != ncol(H)) stop("Intial state should have same
                                            length as the number of columns in H")


  TT <- nrow(z)
  d <- ncol(A)

  # Initialize objects for kf iteration
  Kt <- Pt_m1 <- NULL
  other_shortcut_obj <- vector(mode = 'list', length = nrow(H))
  Pt_m1 <- matrix(0, nrow = d, ncol = d)
  kf_update_obj <- kf_initialize(d, TT, initial_state, nrow(H))

  for(t in 2:TT){

    # 1. predict
    P_p <- A %*% Pt_m1 %*% t(A) + W

    # 2. update
    S      <- H %*% P_p %*% t(H) + Q
    Sm1    <- solve(S)
    Kt     <- P_p %*% t(H) %*% Sm1
    Pt_m1  <- (diag(1, ncol = d, nrow = d) - Kt %*% H) %*% P_p

    for(ii in 1:nrow(H)){
      other_shortcut_obj[[ii]] <- kf_loo_objects(H, Sm1, z, ii, t)
      kf_update_obj[[ii]] <- kf_update_kernel(A, W, kf_update_obj[[ii]],
                                              other_shortcut_obj[[ii]], t, d)
    }
  }
  x_hat <- lapply(kf_update_obj, function(x) x$x_hat)
  return(pred_state = x_hat)
}



#' Decoding trajectory with short cut formula for adding one equation in
#'
#' Function to decode the trajectory using a short cut formula for adding
#'  one equation in.
#'
#' @param A state matrix
#' @param W state covariance matrix
#' @param H encoding matrix for the full model
#' @param Q encoding covariance matrix for the full model
#' @param length_current length of current model
#' @param z centered spike count matrix
#' @param initial_state centered initial velocity state
#' @return  a list of decoded trajectory for a given trial obtained by adding
#' one equation at the time. Every element of the list corresponds to an equation
#' added
kalman_decode_aoo <- function(A = NULL, W = NULL, H = NULL, Q = NULL,
                             length_current, z, initial_state = rep(0, ncol(H))){


  # Edge cases
  if(is.null(W) | is.null(A) | is.null(H) | is.null(Q))
    stop("Neither A nor W nor H nor Q can ever be null")
  if(is.null(length_current)) stop("length of current model should be specified")
  if(nrow(A) != ncol(A)) stop("A should be squared")
  if(nrow(W) != ncol(W)) stop("W should be squared")
  if(nrow(A) != ncol(W)) stop("A and W should have same dimensions")
  if(nrow(A) != ncol(H)) stop("H should have same number of column as the square
                              matrices A and W")
  if(is.null(H) | is.null(Q)) stop("Neither H nor Q can ever be null")
  if(nrow(Q) != ncol(Q)) stop("Q should be a square matrix")
  if(nrow(H) != ncol(Q)) stop("H should have as many rows and Qm1")
  if(is.null(z) | !is.matrix(z)) stop("z has to be a matrix")
  if(length(initial_state) != ncol(H)) stop("Intial state should have same
                                            length as the number of columns in H")


  Hfull <- H
  Qfull <- Q

  TT <- nrow(z)
  d  <- ncol(A)
  incl <- length_current
  left <- nrow(Hfull) - incl

  # get the current model parameters
  H <- Hfull[1:incl,]
  Q <- Qfull[1:incl, 1:incl]

  # each element of the list is a matrix
  Kt <- Pt_m1 <- NULL
  other_shortcut_obj <- vector(mode = 'list', length = left)
  kf_update_obj <- kf_initialize(d, TT, initial_state, left)

  # know initial state for sure
  Pt_m1 <- matrix(0, nrow = d, ncol = d)

  for(t in 2:TT){

    # 1. predict
    P_p <- A %*% Pt_m1 %*% t(A) + W

    # 2. update
    S      <- H %*% P_p %*% t(H) + Q
    Sm1    <- solve(S)
    Kt <- P_p %*% t(H) %*% Sm1
    Pt_m1 <- (diag(1, ncol = d, nrow = d) - Kt %*% H) %*% P_p

    # compute add-one in predicton at time t
    for(ii in 1:left){
      other_shortcut_obj[[ii]] <- kf_aoo_objects(Hfull, Qfull, Sm1, incl, z, ii,
                                                 t, P_p)
      kf_update_obj[[ii]] <- kf_update_kernel(A, W, kf_update_obj[[ii]],
                                              other_shortcut_obj[[ii]], t, d)
    }
  }
  x_hat <- lapply(kf_update_obj, function(x) x$x_hat)
  return(pred_state = x_hat)
}

#' Initialize kf objects
#'
#' This function initialize the kf objects
#' @param d dimension of the prediction object, should be 3
#' @param TT length of the trial time
#' @param initial_state initial state for the kf
#' @param n_row how many element to initialize, should match number of equations
#' to leave out/add in
#' @return a list of initialized
kf_initialize <- function(d, TT, initial_state, n_row){

  if(is.null(d) | is.null(TT) | is.null(initial_state) | is.null(n_row))
    stop("None of the variables can be NULL")
  if(d != length(initial_state)) stop("d and length of initial state should match")

  # list for final predictions
  x_hat   <- vector("list", n_row)
  x_hat_p <- vector("list", n_row)
  kf_update_obj <- list()

  for(i in 1:n_row){
    x_hat     <- matrix(NA, nrow = TT, ncol = d)
    x_hat_p   <- matrix(NA, nrow = TT, ncol = d)
    x_hat[1,] <- initial_state
    Pt_m1_i   <- matrix(0, nrow = d, ncol = d)
    kf_update_obj[[i]] <- list(x_hat_p = x_hat_p, x_hat = x_hat,
                               Pt_m1_i = Pt_m1_i)
  }
  return(kf_update_obj)
}


kf_objects <- function(Hfull, Qfull, Sm1, incl, z, index, t, P_p){

  if(is.null(Hfull) | is.null(Qfull) | is.null(Sm1) | is.null(incl) |
     is.null(z) | is.null(index) | is.null(t))
    stop("None of the variables can be NULL")
  if((incl + index) > nrow(Hfull)) stop("incl + index can't be greater than
                                        the number of rows in H")

  index_sel <- c(1:incl, incl + index)
  Hi  <- Hfull[index_sel,]
  Qi  <- Qfull[index_sel, index_sel]
  Sp1 <- Hi %*% P_p %*% t(Hi) + Qi
  Si  <- get_inverse_matrix_p1(Sp1, Sm1)
  zt  <- z[t,index_sel]

  return(list(Hi = Hi, Si = Si, zt = zt))
}


#' Prepare the add-one-in objects
#'
#' This function prepare the add-one-in objects for the kf update
#' @param Hfull H parameter-matrix with all equations to add
#' @param Qfull Q covarance-matrix with all equations to add
#' @param incl length of included equation
#' @param z spike count matrix
#' @param index index to include
#' @param t time step
#' @param P_p prior covariance matrix from KF formula
#' @return a list with parameters rescaled according to the add-one-in formula
kf_aoo_objects <- function(Hfull, Qfull, Sm1, incl, z, index, t, P_p){

  if(is.null(Hfull) | is.null(Qfull) | is.null(Sm1) | is.null(incl) |
     is.null(z) | is.null(index) | is.null(t))
    stop("None of the variables can be NULL")
  if((incl + index) > nrow(Hfull)) stop("incl + index can't be greater than
                                        the number of rows in H")

  index_sel <- c(1:incl, incl + index)
  Hi  <- Hfull[index_sel,]
  Qi  <- Qfull[index_sel, index_sel]
  Sp1 <- Hi %*% P_p %*% t(Hi) + Qi
  Si  <- get_inverse_matrix_p1(Sp1, Sm1)
  zt  <- z[t,index_sel]

  return(list(Hi = Hi, Si = Si, zt = zt))
}


#' Prepare the leave-one-out objects
#'
#' This function prepare the add-one-in objects for the kf update
#' @param H H parameter-matrix with all the equations
#' @param Q Q covariance-matrix with all the equations
#' @param z spike count matrix
#' @param index index to exclude
#' @param t time step
#' @return a list with parameters rescaled according to the leave-one-out formula
kf_loo_objects <- function(H, Sm1, z, index, t){

  if(is.null(H) | is.null(Sm1) | is.null(z) | is.null(index) |
     is.null(t)) stop("None of the variables can be NULL")
  if((index) > nrow(H)) stop("index can't be greater than the number of rows in H")

  Hi <- H[-index,]
  Si <- get_inverse_matrix_m1(Sm1, ind = index)
  zt <- z[t,-index]

  return(list(Hi = Hi, Si = Si, zt = zt))
}


#' Updates kalman filter solution regardless leave-one-out/add-one-in implementation
#'
#' This function updates kalman filter solution regardless
#' leave-one-out/add-one-in implementation
#' @param A prior model parameter matrix
#' @param W prior model covariance matrix
#' @param kf_obj_i kf objects at time t-1 for all the equations
#' @param other_shortcut_obj objects containing specific leave-one/add-one stuff
#' @param t time step
#' @param d dimension of the state
#' @return a list with parameters rescaled according to the leave-one-out formula
kf_update_kernel <- function(A, W, kf_obj_i, other_shortcut_obj, t, d){

  # Edge cases
  if(!is.list(kf_obj_i)) stop('kf_obj_i is expected to be a list')
  if(!("x_hat_p" %in% names(kf_obj_i)) |
     !("x_hat" %in% names(kf_obj_i)) |
     !("Pt_m1_i" %in% names(kf_obj_i))) stop('kf_obj_i is expected to be
                                             a list with 3 elements:
                                             x_hat_p, x_hat, Pt_m1_i')
  if(d != nrow(A)) stop("d and dimension of A should match")


  # Extract objects from update kf object
  x_hat_p <- kf_obj_i$x_hat_p
  x_hat   <- kf_obj_i$x_hat
  Pt_m1_i <- kf_obj_i$Pt_m1_i

  # get new parameters
  Hi <- other_shortcut_obj$Hi
  Si <- other_shortcut_obj$Si
  zt <- other_shortcut_obj$zt

  # 1. predict
  P_p_i <- A %*% Pt_m1_i %*% t(A) + W
  x_hat_p[t,] <- c(A %*% x_hat[(t-1),])

  Kt_i <- P_p_i %*% t(Hi) %*% Si
  x_hat[t,] <- x_hat_p[t,] + c( Kt_i %*% ( zt - c(Hi %*% x_hat_p[t,])) )
  Pt_m1_i <- (diag(1, ncol = d, nrow = d) - Kt_i %*% Hi) %*% P_p_i

  return(list(x_hat_p = x_hat_p, x_hat = x_hat, Pt_m1_i = Pt_m1_i))
}

