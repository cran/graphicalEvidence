#' rmatrix_last_col_fixed
#'
#' This is a private function that calls C++ code to run the MCMC sampler and
#' compute the average of equation 11 in associated paper. This function uses a
#' restricted Hao Wang sampler to update omega one row/column at a time,
#' iterating from 1 to p - 1 for each of the burnin + nmc iterations of the MCMC
#' sampler.
#' 
#' @param S The sample covariance matrix for the currently considered columns
#' of xx
#' @param n the number of samples of xx
#' @param burnin the number of iterations to run in the MCMC sampler before 
#' saving values for the calculation of MC average of equation 11
#' @param nmc the number of iterations that the MCMC sampler uses to calculate
#' equation 11
#' @param alpha the alpha parameter for G-Wishart prior
#' @param fixed_last_col the ith column excluding the ith element of the
#' currently considered precision matrix
#' @param scale_matrix the scale matrix parameter for the G-Wishart prior
#' @param prior the specified prior distribution, one of 'GHS', 'BGL', 'Wishart'
#' @param dof the dof parameter for the Wishart prior
#' @param matrix_acc_gibbs the accumulated changes of previous calls
#' to updating omega created by storing a modified outer product of the last
#' column that is used to update indices which have 0's in the adjacency matrix
#' during the MCMC sampler
#' @param post_mean_omega the posterior mean of the sampled precision matrix
#' in the unrestricted sampler
#' @param post_mean_tau the posterior mean of the sampled tau matrix
#' in the unrestricted sampler
#' @param lambda the lambda parameter for either BGL or GHS prior
#' 
#' @returns A list that stores the updated omega matrix, the initial starting
#' matrix for the next iteration of the unrestricted sampler, as well as the
#' computed MC average of equation 11 
#' @keywords internal
#' @noRd
rmatrix_last_col_fixed <- function(
  S,
  n,
  burnin,
  nmc,
  fixed_last_col,
  prior,
  dof = 0,
  matrix_acc_gibbs = 0,
  post_mean_omega = 0,
  post_mean_tau = 0,
  lambda = 0
) {
  
  # Get arguments for C++ call
  p <- nrow(S)
  
  coded_prior <- switch(
    prior,
    'Wishart' = 0,
    'BGL' = 1,
    'GHS'= 2
  )
  
  ##################################
  # RcppArmadillo implementation  ##
  ##################################

  # RcppArmadillo implementation
  ans_last_col <- mcmc_last_col_rmatrix(
    n, burnin, nmc, p, dof, lambda, coded_prior, fixed_last_col, S,
    post_mean_tau, matrix_acc_gibbs, post_mean_omega
  )

  ##################################

  return(list(
    MC_avg_eq_11=ans_last_col[[1]],
    start_point_first_gibbs=ans_last_col[[2]],
    post_mean_omega_22=ans_last_col[[3]]
  ))
  
}