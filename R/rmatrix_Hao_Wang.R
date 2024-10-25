#' rmatrix_Hao_Wang
#'
#' This is a private function that calls C++ code to run the MCMC sampler and
#' compute the average of equation 9 in associated paper. This function uses an
#' unrestricted Hao Wang sampler to update omega one row/column at a time,
#' iterating from 1 to p - 1 for each of the burnin + nmc iterations of the MCMC
#' sampler.
#' 
#' @param S The sample covariance matrix for the currently considered columns
#' of xx
#' @param n the number of samples of xx
#' @param burnin the number of iterations to run in the MCMC sampler before 
#' saving values for the calculation of MC average of equation 9
#' @param nmc the number of iterations that the MCMC sampler uses to calculate
#' equation 11
#' @param prior the specified prior distribution, one of 'GHS', 'BGL', 'Wishart'
#' @param dof the dof parameter for the Wishart prior
#' @param matrix_acc_gibbs the accumulated changes of previous calls
#' to updating omega created by storing a modified outer product of the last
#' column that is used to update indices which have 0's in the adjacency matrix
#' during the MCMC sampler
#' @param lambda the lambda parameter for either BGL or GHS prior
#' 
#' @returns A list that stores the updated omega matrix, tau matrix, as well as
#' the computed MC average of equation 9
#' @keywords internal
#' @noRd
rmatrix_Hao_Wang <- function(
  S,
  n,
  burnin,
  nmc,
  prior,
  dof = 0,
  matrix_acc_gibbs = 0,
  lambda = 0
) {

  coded_prior <- switch(
    prior,
    'Wishart' = 0,
    'BGL' = 1,
    'GHS'= 2
  )

  p <- nrow(S)
  
  ##################################
  # RcppArmadillo implementation  ##
  ##################################

  hw_results <- mcmc_hw_rmatrix(
    n, burnin, nmc, p, coded_prior, dof, lambda, S, matrix_acc_gibbs
  )
  
  ##################################

  return(list(
    post_mean_omega=hw_results[[1]],
    post_mean_tau=hw_results[[2]],
    MC_avg_eq_9=hw_results[[3]]
  ))
}