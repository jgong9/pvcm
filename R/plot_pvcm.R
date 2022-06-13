#' Trajectory plots for vayring coefficient functions and principal functions
#'
#' @export plot_pvcm
#' @param pvcm_fit A list output of pvcm function
#' @param pallet_size A vector of two positive integer elements for the number of rows and columns of pallet
#' @param U_sample A numeric vector of index variable
#' @param beta_plot A logical value to indicate whether the coefficient functions (TRUE) or the principal functions (FALSE) are displayed
#' @param true_covariates A vector of names for the known true covariates. Default NULL will display the plots for all selected covariates when beta_plot is TRUE

plot_pvcm <- function(pvcm_fit, pallet_size = c(1,1), X_sample = NULL, U_sample, beta_plot = TRUE, true_covariates = NULL){

  par(mfrow=pallet_size)


  if(beta_plot){
    p <- dim(pvcm_fit$Theta)[1]
    selected_covariates_loc <- apply(pvcm_fit$Theta_sparse,1,FUN=function(x){
      sum(x) !=0
    })

    if(is.null(true_covariates)){
      for(i in (1:p)[selected_covariates_loc] ){

        beta_all_mat <- pvcm_fit$Theta_sparse%*% pvcm_fit$B_U_mat

        plot(sort(U_sample), ( pvcm_fit$Theta_sparse[i,] %*% pvcm_fit$B_U_mat )[order(U_sample)], col = "red",

             ylim= c(
               min(beta_all_mat),
               max(beta_all_mat)

             ),

             ylab=expression(hat(beta)[j]), xlab="U",
             main=(pvcm_fit$selected_covariates)[i],
             type="l",
             lwd= 1.2,
             cex=1.2,
             cex.axis=1.2,
             cex.lab=1.1
        )

      }
    } else {



      for(i in (1:p)[selected_covariates_loc & (colnames(X_sample) %in% true_covariates) ] ){

        beta_all_mat <- pvcm_fit$Theta_sparse%*% pvcm_fit$B_U_mat

        plot(sort(U_sample), ( pvcm_fit$Theta_sparse[i,] %*% pvcm_fit$B_U_mat )[order(U_sample)], col = "red",

             ylim= c(
               min(beta_all_mat),
               max(beta_all_mat)

             ),

             ylab=expression(hat(beta)[j]),
             xlab="U",
             main=colnames(X_sample)[i],
             type="l",
             lwd= 1.2,
             cex=1.2,
             cex.axis=1.2,
             cex.lab=1.1
        )

      }
    }

  } else {
    for(i in 1:pvcm_fit$rank){
      plot(sort(U_sample), ( pvcm_fit$D[i] * t(pvcm_fit$V[,i]) %*% pvcm_fit$B_U_mat )[order(U_sample)], col = "red",
           ylim = c(min(( pvcm_fit$D %*% t(pvcm_fit$V) %*% pvcm_fit$B_U_mat )),max(( pvcm_fit$D %*% t(pvcm_fit$V) %*% pvcm_fit$B_U_mat ))),
           type = "l", main = paste("Principal function ",i),
           ylab= "", xlab="U",
           lwd= 1.2,
           cex=1.2,
           cex.axis=1.2,
           cex.lab=1.2
      )
    }
  }



}
