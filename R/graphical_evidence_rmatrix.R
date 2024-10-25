#' @title Compute Marginal Likelihood using Graphical Evidence for Wishart,
#' BGL, and GHS
#' 
#' @description
#' Computes the marginal likelihood of input data xx under one of the following
#' priors: Wishart, Bayesian Graphical Lasso (BGL), and Graphical Horseshoe 
#' (GHS), specified under prior_name. 
#' 
#' @param xx The input data specified by a user for which the marginal 
#' likelihood is to be calculated. This should be input as a matrix like object
#' with each individual sample of xx representing one row
#' @param burnin The number of iterations the MCMC sampler should iterate 
#' through and discard before beginning to save results
#' @param nmc The number of samples that the MCMC sampler should use to estimate
#' marginal likelihood
#' @param prior_name The name of the prior for which the marginal should be 
#' calculated, this is one of 'Wishart', 'BGL', 'GHS'
#' @param lambda A number specifying lambda for the priors of 'BGL' and 'GHS'
#' prior
#' @param alpha A number specifying alpha for the priors of 'Wishart'
#' @param V The scale matrix when specifying 'Wishart'
#' @param print_progress A boolean which indicates whether progress should be 
#' displayed on the console as each row of the telescoping sum is computed
#' 
#' @returns An estimate for the marginal likelihood under specified prior with
#' the specified parameters
#' 
#' @examples
#' # Compute the marginal likelihood of xx for GHS prior using 1,000 
#' # burnin and 5,000 sampled values at each call to the MCMC sampler
#' g_params <- gen_params_evidence('GHS')
#' marginal_results <- graphical_evidence_rmatrix(
#'   g_params$x_mat, 1e3, 5e3, 'GHS', lambda=1
#' )
#' @export
graphical_evidence_rmatrix <- function(
  xx,
  burnin,
  nmc,
  prior_name = c('Wishart', 'BGL', 'GHS'),
  lambda = 0,
  alpha = 0,
  V = 0,
  print_progress = FALSE
) {
  
  # Ensure arguments are matrices
  xx <- as.matrix(xx)
  V <- as.matrix(V)
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Calculate sample covariance matrix
  S <- as.matrix(t(xx) %*% xx)
  
  # Extract n and p
  n <- nrow(xx)
  p <- ncol(xx)
  
  # Ensure scale matrix is matrix type
  V <- as.matrix(V)
  
  # Initialize storage for I_{p-j+1} -IV_{p-j+1} for every j from 1 to p
  log_ratio_of_likelihood <- numeric(p - 1)
  
  # Errors possible when calculating log-marginal -> log f(y_{1:p}),
  # return NA in that case
  result <- tryCatch({
    
    # storage of posterior means
    last_col_store <- vector("list", length = p)
    
    # Accumulator for shifts
    matrix_acc <- matrix(numeric(p * p), nrow=p)
    
    # Main loop from 1 to p, for Wishart only this can be run in parallel
    if ((prior_name != 'Wishart') | (p < 15)) { 
      for (num_rmat in 1:(p - 1)) {
        
        updated_density <- graphical_evidence_rmatrix_iteration(
          xx, S, n, p, burnin, nmc, prior_name, lambda, alpha, V, print_progress,
          matrix_acc, num_rmat
        )
        
        log_ratio_of_likelihood[num_rmat] <- updated_density$log_ratio_of_likelihood
        last_col_store[[num_rmat]] <- updated_density$last_col_store_item
        matrix_acc <- updated_density$matrix_acc
      }

      # Calculate final iteration
      direct_eval_log_prior_density <- graphical_evidence_rmatrix_last_iter(
        xx, S, n, p, nmc, prior_name, lambda, alpha, V, last_col_store, matrix_acc
      )
      
      return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    }
    
    # Use doParallel to execute top level loop for Wishart case
    else {
      par_res <- foreach::foreach(num_rmat = 1:(p - 1)) %dopar% {
        
        graphical_evidence_rmatrix_iteration(
          xx, S, n, p, burnin, nmc, prior_name, lambda, alpha, V, print_progress,
          matrix_acc, num_rmat
        )
      }
      
      # Aggregate parallel computation results
      for (i in 1:(p - 1)) {
        last_col_store[[i]] <- par_res[[i]]$last_col_store_item
        log_ratio_of_likelihood[i] <- par_res[[i]]$log_ratio_of_likelihood
      }
      
      # Calculate final iteration
      direct_eval_log_prior_density <- graphical_evidence_rmatrix_last_iter(
        xx, S, n, p, nmc, prior_name, lambda, alpha, V, last_col_store, 
        matrix_acc
      )

      return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    }
    
  }, error = function(e) {
    print(e)
    return(NaN)
  })
  
  return(result)
}
