#' Davis-Yin Splitting for PVCM
#'
#' @export dys_pvcm
#' @param y_vec A numeric vector of response variable
#' @param M_mat A matrix whose row is the kroneker product between the covariates X_i and the evaluated basis functions at U_i
#' @param tau A vector of the two tuning parameters
#' @param Theta_initial An initial value for the Theta matrix
#' @param weights_ann A vector of weights for the adaptive nuclear norm penalty
#' @param weights_agl A vector of weights for the adaptive group lasso penalty
#' @param gamma_val A numeric value for one over Lipschitz constant
#' @param max.iter A maximum iteration number
#' @param tol A numeric value for tolerance of relative deviation in order to stop DYS iteration
#' @return A list that contains estimated Theta matrices with different characteristics

dys_pvcm <- function(y_vec, M_mat, tau, Theta_initial, weights_ann, weights_agl, gamma_val, max.iter=100, tol=1e-5){

  n <- length(y_vec)
  p <- dim(Theta_initial)[1]
  d <- dim(Theta_initial)[2]


  n.iter <- 1

  converge <- 0

  rel_diff_vec <- rep(NA, max.iter)
  rank_trace <- c()

  X_k <- Theta_initial

  while(n.iter <= max.iter){    ## iteration
    Temp <- X_k

    # Update Y
    Y_k <- X_k
    Y_kplus1 <- X_k

    for(j in 1:p){
      y_j <- Y_k[j,]
      y_j_norm <- sqrt( t(y_j) %*% y_j )

      if(tau[2]==0 ){
        y_j_update <-  y_j
      } else {
        if( y_j_norm <= (tau[2] * weights_agl[j] * gamma_val) ){
          y_j_update <- rep(0, d )
        } else {
          y_j_update <- c(1 - ( tau[2] * weights_agl[j] * gamma_val ) / y_j_norm ) * y_j
        }
      }

      Y_kplus1[j,] <- y_j_update
    }
    ## Y update complete

    ## Update Z
    delta_H_y_t <-  1/n * crossprod( M_mat, M_mat %*% c(Y_kplus1) - y_vec )
    X_mat <- 2 * Y_kplus1 - gamma_val * matrix(delta_H_y_t, p, d) - X_k


    # Solve the proximal operator
    svd_X_mat <- svd(X_mat)

    temp_d <- (svd_X_mat$d - gamma_val * tau[1] * weights_ann)
    temp_d[ temp_d <= 0 ] <- 0
    rank <- sum( temp_d > 0 )

    Z_kplus1 <- svd_X_mat$u %*% diag(temp_d) %*% t(svd_X_mat$v)
    if( rank == 0){
      U <-NULL
      V <- NULL
      D <- NULL
    } else {
      U <- svd_X_mat$u[,1:rank]
      V <- svd_X_mat$v[,1:rank]
      D <- temp_d[1:rank]
    }


    X_kplus1 <- X_k +  Z_kplus1 - Y_kplus1

    ## New X_k for the next iteration
    X_k <- X_kplus1

    rel_diff <- norm(X_k - Temp, "F") / max(1, norm(Temp, "F"))
    rel_diff_vec[n.iter] <- rel_diff
    error.control <- (rel_diff <= tol)


    if(error.control){
      converge<-1
      break
    }

    n.iter <- n.iter + 1
    rank_trace <- c(rank_trace, rank)
  }


  return(
    list(
      "Theta_sparse" = Y_kplus1, "Theta_low_rank" = Z_kplus1, "Theta" = X_kplus1,
      "rank" = rank, "U" = U, "V"=V, "D"= D, "converge"=converge, "rank_trace"=rank_trace,
      "iteration" = n.iter, "rel_diff"= rel_diff_vec
    )
  )
}
