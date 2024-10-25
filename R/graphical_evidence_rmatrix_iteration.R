#' graphical_evidence_rmatrix_iteration
#' 
#' Compute the log ratio of likelihood for the  of input data xx, labeled as 
#' term III in associated paper for one iteration.
#' 
#' @param xx The input data specified by a user for which the marginal 
#' likelihood is to be calculated. This should be input as a matrix like object
#' with each individual sample of xx representing one row
#' @param S the sample covariance matrix of the data xx
#' @param n the number of samples, in this case rows, of the data xx
#' @param p the dimension, in this case the number of columns, of the data xx
#' @param burnin The number of iterations the MCMC sampler should iterate 
#' through and discard before beginning to save results
#' @param nmc The number of samples that the MCMC sampler should use to estimate
#' marginal likelihood
#' @param prior The name of the prior for which the marginal should be 
#' calculated, this is one of 'Wishart', 'BGL', 'GHS'
#' @param lambda A number specifying lambda for the priors of 'BGL' and 'GHS'
#' prior
#' @param alpha A number specifying alpha for the priors of 'Wishart'
#' @param V The scale matrix when specifying 'Wishart'
#' @param print_progress A boolean which indicates whether progress should be 
#' displayed on the console as each row of the telescoping sum is computed
#' @param matrix_acc the accumulated changes of previous calls
#' to updating omega created by storing a modified outer product of the last
#' column
#' @param num_rmat the current iteration of the algorithm from 1 to p
#' 
#' @returns A list that stores III result, the fixed last column in the
#' restricted sampler, and the update matrix accumulator 
#' @keywords internal
#' @noRd
graphical_evidence_rmatrix_iteration <- function(
  xx,
  S,
  n,
  p,
  burnin,
  nmc,
  prior,
  lambda,
  alpha,
  V,
  print_progress,
  matrix_acc,
  num_rmat
) {
  
  # Print each iteration to track progress if requested
  if (print_progress) {
    cat(paste0("Working on ", num_rmat, "th row of the telescoping sum\n"))
  }
  
  # Reduce data by one column
  reduced_data_xx <- xx[, 1:(p - num_rmat + 1)]
  p_reduced <- ncol(reduced_data_xx)
  S_reduced <- as.matrix(t(reduced_data_xx) %*% reduced_data_xx)
  
  # General case

  # Matrix accumulator is not needed in Wishart, will never be used
  matrix_acc_gibbs <- matrix_acc[1:p_reduced, 1:p_reduced]
  
  # Run unrestricted Hao Wang sampler which will
  # be used to approximate the Normal density in the
  # evaluation of the term IV_{p-j+1} - Equation (9)
  hw_results <- rmatrix_Hao_Wang(
    S_reduced, n, burnin, nmc, prior, dof=(alpha + 1 - num_rmat),
    matrix_acc_gibbs=matrix_acc_gibbs, lambda=lambda
  )
  
  fixed_last_col <- (
    hw_results$post_mean_omega[1:(p_reduced - 1), p_reduced]
  )
  
  last_col_results <- rmatrix_last_col_fixed(
    S_reduced, n, burnin, nmc, fixed_last_col, prior,
    dof=(alpha + 1 - num_rmat), matrix_acc_gibbs=matrix_acc_gibbs,
    post_mean_omega=hw_results$post_mean_omega,
    post_mean_tau=hw_results$post_mean_tau, lambda=lambda
  )
  
  # Computing IV_{p-j+1}
  log_posterior_density <- (
    hw_results$MC_avg_eq_9 + last_col_results$MC_avg_eq_11
  )
  
  # Computing I_{p-j+1}
  st_dev <- sqrt(1 / last_col_results$post_mean_omega_22)
  mean_vec <- -(st_dev^2) * (
    reduced_data_xx[, 1:(p_reduced - 1)] %*% 
      as.matrix(hw_results$post_mean_omega[p_reduced, 1:(p_reduced - 1)])
  )
  log_data_likelihood <- sum(log(
    dnorm(reduced_data_xx[, p_reduced], mean_vec, st_dev)
  ))
  
  # Computing III_{p-j+1}
  log_ratio_of_likelihood <- (
    log_data_likelihood - log_posterior_density
  )
  
  # Not needed for Wishart
  if (prior != 'Wishart') {
    matrix_acc[1:(p - num_rmat), 1:(p - num_rmat)] <- (
      matrix_acc[1:(p - num_rmat), 1:(p - num_rmat)] +
        (1 / last_col_results$post_mean_omega_22) * (
          fixed_last_col %*% t(fixed_last_col)
        )
    )
  }
  
  last_col_store_item <- c(
    fixed_last_col, last_col_results$post_mean_omega_22
  )
  
  return(list(
    log_ratio_of_likelihood=log_ratio_of_likelihood,
    last_col_store_item=last_col_store_item,
    matrix_acc=matrix_acc
  ))
}