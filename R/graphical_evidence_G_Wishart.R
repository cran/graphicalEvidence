#' @title Compute Marginal Likelihood using Graphical Evidence under G Wishart
#' 
#' @description
#' Computes the marginal likelihood of input data xx under G-Wishart prior using
#' graphical evidence.
#' 
#' @param xx The input data specified by a user for which the marginal 
#' likelihood is to be calculated. This should be input as a matrix like object
#' with each individual sample of xx representing one row
#' @param burnin The number of iterations the MCMC sampler should iterate 
#' through and discard before beginning to save results
#' @param nmc The number of samples that the MCMC sampler should use to estimate
#' marginal likelihood
#' @param alpha A number specifying alpha for G-Wishart prior
#' @param V The scale matrix of G-Wishart prior
#' @param G The adjacency matrix of G-Wishart prior
#' @param print_progress A boolean which indicates whether progress should be 
#' displayed on the console as each row of the telescoping sum is computed
#' 
#' 
#' @returns An estimate for the marginal likelihood under G-Wishart prior with
#' the specified parameters
#' 
#' @examples
#' # Compute the marginal likelihood of xx for G-Wishart prior using 
#' # 2,000 burnin and 10,000 sampled values at each call to the MCMC sampler
#' g_params <- gen_params_evidence('G_Wishart')
#' marginal_results <- graphical_evidence_G_Wishart(
#'   g_params$x_mat, 2e3, 1e4, 2, g_params$scale_mat, g_params$g_mat
#' )
#' @export
graphical_evidence_G_Wishart <- function(
  xx,
  burnin,
  nmc,
  alpha,
  V,
  G,
  print_progress = FALSE
) {
  
  # Ensure arguments are matrices
  xx <- as.matrix(xx)
  V <- as.matrix(V)
  G <- as.matrix(G)
  
  # Calculate sample covariance matrix
  S <- as.matrix(t(xx) %*% xx)
  
  # Extract n and p
  n <- nrow(xx)
  p <- ncol(xx)
  
  # Initialize storage locations
  log_ratio_of_liklelihoods <- rep(0, p)
  
  log_data_likelihood <- rep(0, p)
  log_posterior_density <- rep(0, p)
  log_normal_posterior_density <- rep(0, p)
  log_gamma_posterior_density <- rep(0, p)
  
  # Terms III_{p-j+1}
  log_prior_density <- rep(0, p)
  
  # Gamma density in Equation (20)
  log_prior_density_scalar_gamma <- rep(0, p)
  
  # Normal density in Equation (20)
  log_prior_density_vector_normal <- rep(0, p)
  
  # Cumulative linear shifts storage
  matrix_accumulator <- matrix(0, nrow=p, ncol=p)
  
  # Get initial starting point
  start_point_first_gibbs <- as.matrix(prior_sampler_G_Wishart(
    p, 9, 1, as.matrix(G), as.matrix(V), alpha
  ))

  # Main graphical evidence loop
  for (num_G_Wishart in 1:p) {
    
    if (num_G_Wishart <= (p - 1)) {
      if (print_progress) { 
        cat(paste0(num_G_Wishart, 'th num_G_Wishart\n'))
      }
      
      # for every iteration we need smaller and smaller blocks of
      # the data matrix xx. When num_G_Wishart = 1, we need the entire
      # matrix xx. When num_G_Wishart = 2, we need the (p-1) x (p-1)
      # block of the matrix xx
      reduced_data_xx <- as.matrix(xx[, 1:(p - num_G_Wishart + 1)])
      p_reduced <- ncol(reduced_data_xx)
      S_reduced <- as.matrix(t(reduced_data_xx) %*% reduced_data_xx)
      
      matrix_accumulator_gibbs <- as.matrix(
        matrix_accumulator[1:p_reduced, 1:p_reduced]
      )
      G_mat_adj_reduced <- as.matrix(
        G[1:p_reduced, 1:p_reduced]
      )
      scale_matrix_reduced <- as.matrix(
        V[1:p_reduced, 1:p_reduced]
      )
      
      # Run the unrestricted sampler to get samples, which will
      # be used to approximate the Normal density in the
      # evaluation of the term IV_{p-j+1} - Equation (21)
      Hao_wang_results <-  G_wishart_Hao_wang(
        S_reduced, n, burnin, nmc, alpha, scale_matrix_reduced,
        G_mat_adj_reduced, matrix_accumulator_gibbs, start_point_first_gibbs
      )

      post_mean_omega <- Hao_wang_results$post_mean_omega
      MC_average_Equation_9 <- Hao_wang_results$MC_average_Equation_9
      
      fixed_last_col <- post_mean_omega[1:(p_reduced - 1), p_reduced]
      
      last_col_results <- G_wishart_last_col_fixed(
        S_reduced, n, burnin, nmc, alpha, fixed_last_col, scale_matrix_reduced,
        G_mat_adj_reduced, matrix_accumulator_gibbs, post_mean_omega
      )
      
      start_point_first_gibbs <- last_col_results$start_point_first_gibbs
      post_mean_omega_22 <- last_col_results$post_mean_omega_22
      MC_average_Equation_11 <- last_col_results$MC_average_Equation_11
      
      log_normal_posterior_density[num_G_Wishart] <-  MC_average_Equation_9
      
      log_gamma_posterior_density[num_G_Wishart] <-  MC_average_Equation_11
      
      log_posterior_density[num_G_Wishart] <- (
        MC_average_Equation_9 + MC_average_Equation_11
      )
      
      st_dev <- sqrt(1 / post_mean_omega_22)
      
      ind_to_use <- 1:(p_reduced - 1)
      
      mean_vec <- -1 * st_dev^2 * (
        as.matrix(reduced_data_xx[, ind_to_use]) %*% 
        as.matrix(post_mean_omega[p_reduced, ind_to_use])
      )
          

      log_data_likelihood[num_G_Wishart] <- (
          sum(log(dnorm(reduced_data_xx[, p_reduced], mean_vec, st_dev)))
      )
      
      # Computing III_{p-j+1}
      # Prior density at posterior mean
      
      V_mat_22 <- scale_matrix_reduced[p_reduced, p_reduced]
      V_mat_12 <- scale_matrix_reduced[1:(p_reduced - 1), p_reduced]
      vec_accumulator_21 <- -1 * (
        matrix_accumulator_gibbs[1:(p_reduced - 1), p_reduced]
      )
      s_21 <- S[1:(p_reduced - 1), p_reduced]
      
      # Note the -1 above. This is done to make sure to respect the
      # crucial indicator function.
      
      G_mat_current_col <- G_mat_adj_reduced[1:(p_reduced - 1), p_reduced]
      which_ones <- as.integer(G_mat_current_col == 1)
      logi_which_ones <- which(which_ones == 1)
      which_zeros <- as.integer(G_mat_current_col == 0)
      logi_which_zeros <- which(which_zeros == 1)
      
      if (sum(which_ones) >= 1) {
        V_mat_12_required <- V_mat_12[logi_which_ones]
        s_21_required <- s_21[logi_which_ones]
        
        if (sum(which_zeros) >= 1) {
          vec_accumulator_21_required <- vec_accumulator_21[logi_which_zeros]
        }
      }
      
      inv_C <- (1 / post_mean_omega_22) * as.matrix(
        scale_matrix_reduced[1:(p_reduced - 1), 1:(p_reduced - 1)]
      )
        
      if (sum(which_ones) >= 1) {
        inv_C_required <- inv_C[logi_which_ones, logi_which_ones]
        
        if (sum(which_zeros) >= 1) {
          inv_C_not_required <- inv_C[logi_which_zeros, logi_which_ones]
          
          vec_accumulator_21_required_modified <- as.vector(
            t(vec_accumulator_21_required) %*% inv_C_not_required
          )
          
          mu_i_reduced <- -solve(
            inv_C_required,
            (V_mat_12_required + vec_accumulator_21_required_modified)
          )
        }
        else {
          mu_i_reduced <- -solve(
            inv_C_required, V_mat_12_required + s_21_required
          )
        }
        
        mean_vec <- mu_i_reduced
        
        # Normal density in Equation S.6
        log_prior_density_vector_normal[num_G_Wishart] <- (
          log(mvtnorm::dmvnorm(
            post_mean_omega[p_reduced, logi_which_ones], t(mean_vec), 
            solve(inv_C_required)
          ))
        )
      }
      else {
        log_prior_density_vector_normal[num_G_Wishart] <- 0
      }
      
      # Computing the GIG
      if (!sum(which_zeros)) {
        
        # Gamma density in Equation (20)
        log_prior_density_scalar_gamma[num_G_Wishart] <- (
          log(dgamma(
            post_mean_omega_22, alpha + (sum(which_ones) / 2) + 1,
            scale=(2 / V_mat_22)
          ))
        )
      }
      else {
        GIG_a <- V_mat_22
        
        diag_V_mat <- diag(scale_matrix_reduced)
        diag_V_mat_required <- diag_V_mat[logi_which_zeros]
        vec_accumulator_21 <- -1 * matrix_accumulator_gibbs[1:(p_reduced - 1),
                                                            p_reduced]
        vec_accumulator_21_required <- vec_accumulator_21[logi_which_zeros]
        
        GIG_b <- sum(diag_V_mat_required * (vec_accumulator_21_required^2))
        
        # Treat any value below machine precision as 0
        if (abs(GIG_b)^2 < .Machine$double.eps) {
          log_prior_density_scalar_gamma[num_G_Wishart] <- (
            log(dgamma(
              post_mean_omega_22, alpha + (sum(which_ones) / 2) + 1, 
              scale=(2 / V_mat_22)
            ))
          )
        }
        else {
          GIG_p <- alpha + (sum(which_ones) / 2) + 1
          log_prior_density_scalar_gamma[num_G_Wishart] <- (
            (GIG_p /2 ) * log(GIG_a / GIG_b) - log(2) -
            log(besselK(sqrt(GIG_a * GIG_b), GIG_p)) + 
            (GIG_p - 1) * log(post_mean_omega_22) - 
            0.5 * (GIG_a * post_mean_omega_22 + (GIG_b / post_mean_omega_22))
          )
        }
      }
      
      log_prior_density[num_G_Wishart] <- (
        log_prior_density_vector_normal[num_G_Wishart] + 
        log_prior_density_scalar_gamma[num_G_Wishart]
      )
      
      log_ratio_of_liklelihoods[num_G_Wishart] <- (
        log_data_likelihood[num_G_Wishart] + 
        log_prior_density[num_G_Wishart] -
        log_posterior_density[num_G_Wishart]
      )
      
      matrix_accumulator[1:(p - num_G_Wishart), 1:(p - num_G_Wishart)] <- (
        matrix_accumulator[1:(p - num_G_Wishart), 1:(p - num_G_Wishart)] + 
        (1 / post_mean_omega_22) * (fixed_last_col %*% t(fixed_last_col))
      )
    }
    else {
      # To evaluate log f(y_1) or the marginal likelihood of the
      # remaining first column, it's like
      # imposing a univariate G-wishart prior on the precision of a
      # univariate normal random variable

      xx_reduced <- xx[, 1]
      S_reduced <- t(xx_reduced) %*% xx_reduced
      p_reduced <- 1
      
      V_mat_22 <- scale_matrix_reduced[p_reduced, p_reduced]
      
      LOG_marginal_first_col_only <- -(n * p_reduced / 2) * log(pi) + (
        logmvgamma(alpha + (n / 2) + 1, p_reduced) -
        (alpha + (n / 2) + 1) * log(det(V_mat_22 + S_reduced)) - 
        logmvgamma(alpha + 1, p_reduced) + (alpha + 1) * log(V_mat_22)
      )
      
      log_ratio_of_liklelihoods[p] <- LOG_marginal_first_col_only
    }
  }
  return(sum(log_ratio_of_liklelihoods))
}