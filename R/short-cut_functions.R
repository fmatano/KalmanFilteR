#' Compute inverse of a symmetrix matrix, leaving out one row and column
#'
#' This function computes the inverse of a symmetrix matrix, leaving out one row
#' and column using the short-cut formula provived in my thesis
#'
#' @param Qm1 inverse of the current(full) symmetric matrix
#' @param ind index to leave out, or set of indeces to leave out
get_inverse_matrix_m1 <- function(Qm1, ind){

  # Edge cases
  # if(is.null(ind) | !is.numeric(ind))
  #   stop("ind should be a numeric object: either a single index or a sequence")
  # if(is.null(Qm1) | !is.matrix(Qm1)) stop("Qm1 should be a square matrix")
  # if(nrow(Qm1) != ncol(Qm1)) stop("Qm1 should be a square matrix")
  # if((length(ind) >= ncol(Qm1)) | (max(ind) > ncol(Qm1)))
  #   stop("the set of indexis should be contained in the number of rows in Qm1")


  # Schur-complement
  SA_m1 <- Qm1[ind,ind]
  M     <- Qm1[-ind, -ind]
  L     <- matrix(Qm1[-ind,ind], nrow=nrow(Qm1)-length(ind), ncol=length(ind))

  # leave-one-out precision matrix
  Qi_m1 <- M - (L%*%solve(SA_m1)%*%t(L))
  return(Qi_m1)
}


#' Compute inverse of a symmetrix matrix, adding one row and column
#'
#' This function computes the inverse of a symmetrix matrix, adding one row and column
#' using the short-cut formula provived in my thesis
#'
#' @param Qi full matrix composed of n+1 x n+1 element, the last one is
#' the one to add
#' @param Qm1 inverse of the current used symmetric matrix
get_inverse_matrix_p1 <- function(Qi, Qm1){

  # Edge cases
  if(is.null(Qm1) | !is.matrix(Qm1) | (nrow(Qm1) != ncol(Qm1)))
    stop("Qm1 should be a square matrix")
  if(is.null(Qi) | !is.matrix(Qi) | (nrow(Qi) != ncol(Qi)))
    stop("Qi should be a square matrix")
  if(ncol(Qi) < (ncol(Qm1))) stop('Qi shoudl be a n+k x n+k matrix,
                                       where n is the dimension of Qm1')


  # select the element to add
  l <- (nrow(Qm1) + 1) : nrow(Qi)

  # compute block matrices
  D <- Qi[l, l]
  C <- matrix(Qi[l,-l], nrow=length(l), ncol=nrow(Qi) - length(l))
  SA_m1 <- Qm1[-l,-l]
  Dm1 <- solve(D)
  L   <- -SA_m1 %*% t(C) %*% Dm1
  L_L <- C%*%SA_m1%*%t(C)
  M   <- Dm1 + Dm1 %*% L_L %*% Dm1

  # build the final matrix
  Qi_p1 <- cbind(rbind(SA_m1, t(L)), c(L, M))

  return(Qi_p1)

}
