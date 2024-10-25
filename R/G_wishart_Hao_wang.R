#' G_wishart_Hao_wang
#'
#' This is a private function that calls C++ code to run the MCMC sampler and
#' compute the average of equation 9 in associated paper. This function uses an
#' unrestricted Hao Wang sampler to update omega on row/column at a time,
#' iterating from 1 to p for each of the burnin + nmc iterations of the MCMC
#' sampler.
#' 
#' @param S The sample covariance matrix for the currently considered columns
#' of xx
#' @param n the number of samples of xx
#' @param burnin the number of iterations to run in the MCMC sampler before 
#' saving values for the calculation of MC average of equation 9
#' @param nmc the number of iterations that the MCMC sampler uses to calculate
#' equation 9
#' @param alpha the alpha parameter for G-Wishart prior
#' @param scale_matrix the scale matrix parameter for the G-Wishart prior
#' @param G_mat_adj the adjacency matrix parameter for the G-Wishart prior
#' @param matrix_accumulator_gibbs the accumulated changes of previous calls
#' to updating omega created by storing a modified outer product of the last
#' column that is used to update indices which have 0's in the adjacency matrix
#' during the MCMC sampler
#' @param start_point_first_gibbs the current omega that will be used as a 
#' starting point for the MCMC sampling
#' 
#' @returns A list that stores the updated omega matrix as well as the computed
#' MC average of equation 9
#' @keywords internal
#' @noRd
G_wishart_Hao_wang <- function(
  S,
  n,
  burnin,
  nmc,
  alpha,
  scale_matrix,
  G_mat_adj,
  matrix_accumulator_gibbs,
  start_point_first_gibbs
) {
  
  # Dimension of this iteration
  p <- nrow(S)
  
  # ##################################
  # # RcppArmadillo implementation ###
  # ##################################
  
  ans_hw <- mcmc_hw(
    n, burnin, nmc, alpha, p, S, scale_matrix, G_mat_adj,
    matrix_accumulator_gibbs, start_point_first_gibbs
  )

  return(list(post_mean_omega=ans_hw[[1]], MC_average_Equation_9=ans_hw[[2]]))
}