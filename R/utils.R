##########################################
#####       Utility functions        #####
##########################################


Theta_df_nai <- function(Theta, rank_est){
  p <- dim(Theta)[1]
  d <- dim(Theta)[2]


  index <- apply(Theta==0, 1, sum)
  index <- (index != d)

  p0 <- min(sum(index), p)
  df <- d * rank_est  + p0 * rank_est - rank_est * rank_est
  return(c(df, p - sum(index) ))
}

Theta_BIC <- function(y, M, Theta_fit, rank_est){

  n <- length(y)

  p <- dim(Theta_fit$Theta_sparse)[1]
  d <- dim(Theta_fit$Theta_sparse)[2]

  index <- apply(Theta_fit$Theta_sparse==0, 1, sum)
  index_nonzero_row <- (index != d)
  p0 <- min(sum(index_nonzero_row), p)
  if(rank_est == 0){

    SSE <- sum( (y)^2 )

    df_nai <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[1]
    zero_rows <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[2]
    logMSE <- log(SSE / (n))
    BIC_nai <- SSE / (n)
  } else {
    if(p0 == 0 ){

      SSE <- sum( (y)^2 )

      df_nai <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[1]
      zero_rows <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[2]
      logMSE <- log(SSE / (n))
      BIC_nai <- SSE / (n)
    } else {

      SSE <- sum( (y - M %*% c(Theta_fit$Theta_low_rank)  )^2 )

      df_nai <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[1]
      zero_rows <- Theta_df_nai(Theta_fit$Theta_sparse, rank_est)[2]
      logMSE <- log(SSE / (n))
      BIC_nai <- logMSE + df_nai * log(n) / (n)
    }
  }
  return(c(logMSE = logMSE, BIC_nai = BIC_nai, df_nai = df_nai, zero_rows = zero_rows, SSE=SSE))
}




