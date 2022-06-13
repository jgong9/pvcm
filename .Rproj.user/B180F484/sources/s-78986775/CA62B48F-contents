#' Principal Varying Coefficient Model with Adaptive weights
#'
#' @export pvcm
#' @param y_sample A numeric vector of response variable
#' @param X_sample A design matrix of covariates including an intercept if it is included in your model
#' @param U_sample A numeric vector of index variable
#' @param num_cores The number of CPU cores to implement parallel computing (Default NULL will use the total number of cores - 1)
#' @param knots_val A numeric value for the knots argument of select.knots function in the face R package
#' @param control A list that contains two numbers of tuning parameters except 0 and the stopping criteria for DYS algorithm
#' @return A list that contains all estimation results
#' @importFrom face select.knots
#' @importFrom splines splineDesign
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel detectCores makeCluster parRapply

pvcm <- function(y_sample, X_sample, U_sample, num_cores = NULL, knots_val=7,
                 control = list(num_tau_1=20, num_tau_2=20, max.iter=100, tol=1e-5)
){


  n <- length(y_sample)
  p <- dim(X_sample)[2]


  ## Create orthogonal basis functions with cubic B-spline
  ## with equally spaced points over the domain
  ## ranging from min and max of U_sample
  u_vec <- seq( min(U_sample), max(U_sample),length.out = 100)
  knots_u_vec <- select.knots(u_vec, knots = knots_val ,p=3)

  B_mat <- (splineDesign(knots_u_vec, u_vec, ord = 4))
  d <- dim(B_mat)[2]

  G <- (t(B_mat) %*% B_mat ) / length(u_vec)
  G_sqrt <- eigen(G)$vectors %*% diag( sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  G_sqrt_inverse <- eigen(G)$vectors %*% diag( 1/sqrt(eigen(G)$values)  ) %*% t( eigen(G)$vectors )
  B_U_mat <- t( (splineDesign(knots = knots_u_vec, x = U_sample, ord = 4)) %*% G_sqrt_inverse)


  ## Create the design matrix M with kronecker product
  M_mat <- matrix(NA, n, d* p )
  for(i in 1:n){
    M_mat[i,] <- kronecker((t(B_U_mat))[i,], X_sample[i,] )
  }


  ## Evaluate the initial value for Theta with a ridge regression
  temp <- cv.glmnet(x=M_mat, y=y_sample, alpha=0, intercept = F)
  theta_initial_vec <- coef(temp, s=temp$lambda.min)[-1]
  Theta_initial <- matrix(theta_initial_vec, p, d)
  svd_initial <- svd(Theta_initial)



  ## Adaptive weights
  # 1. Adaptive nuclear norm
  eigen_decom <- eigen( 1/n * t(M_mat) %*%M_mat)
  ANN_weights <- svd_initial$d^(-2)
  gamma_val <- 1/ (eigen_decom$values[1])




  # 2. Adaptive group lasso
  row_norms <- apply( (Theta_initial), 1, FUN = function(x){
    sqrt( sum( x^{2}  ) )
  } )

  AGL_weights <- row_norms^(-1)


  ## Tuning parameters
  tau_1_max <- svd_initial$d[1]^3
  tau_1_min <- min(svd_initial$d)^3
  tau_1_all <- c( 0, exp(seq(from = log(tau_1_min), to = log(tau_1_max), length.out = control[[1]])) )

  tau_2_max <- max(row_norms)^2
  tau_2_min <- min(row_norms)^2
  tau_2_all <- c( 0, exp(seq( log( tau_2_min ), log(tau_2_max), length = control[[2]]) ) )


  tuning = list() # all combination of tuning parameters
  par_mat_prox <- matrix(NA, length(tau_1_all) * length(tau_2_all), 2)
  # for parallel computing
  len_tuning = 0

  for(l in tau_1_all){
    for(k in tau_2_all){
      len_tuning = len_tuning + 1
      par_mat_prox[len_tuning,] <- c(l, k)
    }
  }


  ## Use parallel computing for DYS algorithm with a large number of tuning parameter combinations
  if(is.null(num_cores) == TRUE){
    num_cores <- detectCores() - 1
  }
  my_cluster <- makeCluster(num_cores)



  ## Application of DYS algorithm
  fit_list <- list()
  BIC_fit_list <- list()


  fit_list <- parRapply(my_cluster, par_mat_prox, FUN = dys_pvcm,
                        y_vec = y_sample,
                        M_mat = M_mat,
                        Theta_initial = Theta_initial,
                        weights_ann = ANN_weights,
                        weights_agl = AGL_weights,
                        gamma_val = gamma_val,
                        max.iter = control[[3]],
                        tol = control[[4]]
  )
  for(i in 1:len_tuning ){
    BIC_fit_list[[i]] <- Theta_BIC(y=y_sample, M = M_mat, Theta_fit=fit_list[[i]], rank_est=fit_list[[i]]$rank)
  }

  bic_vec <- vector("numeric", len_tuning)

  for( i in 1:len_tuning){
    bic_vec[i] <-  BIC_fit_list[[i]][2]
  }

  ## The best Theta with respect to BIC
  best_loc <- which( min(bic_vec) == bic_vec )
  best_fit <- fit_list[[best_loc]]

  ## Variable selection
  selected_covariates_loc <- apply(best_fit$Theta_sparse,1,FUN=function(x){
    sum(x) !=0
  })

  if( is.null(colnames(X_sample)) != TRUE){
    selected_covariates <- colnames(X_sample)[selected_covariates_loc]
  } else {
    selected_covariates <- "Cannot find the covariate names in X_sample"
  }

  result <- list( "Theta_low_rank" = best_fit$Theta_low_rank,
                  "Theta_sparse" = best_fit$Theta_sparse,
                  "Theta" = best_fit$Theta,
                  "rank" = best_fit$rank,
                  "U" = best_fit$U,
                  "V"=best_fit$V,
                  "D"= best_fit$D,
                  "B_U_mat" = B_U_mat,
                  "fitted_values" =  diag( X_sample %*% best_fit$Theta_low_rank %*% B_U_mat),
                  "selected_covariates" = selected_covariates)
  return(result)

}




