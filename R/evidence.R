#' @title Compute Marginal Likelihood using Graphical Evidence
#' 
#' @description
#' Computes the marginal likelihood of input data xx under one of the following
#' priors: Wishart, Bayesian Graphical Lasso (BGL), 
#' Graphical Horseshoe (GHS), and G-Wishart, specified under prior_name. 
#' The number of runs is specified by num_runs, where each run is by default
#' using a random permutation of the columns of xx, as marginal likelihood 
#' should be independent of column permutation.
#' 
#' @param xx The input data specified by a user for which the marginal 
#' likelihood is to be calculated. This should be input as a matrix like object
#' with each individual sample of xx representing one row
#' @param burnin The number of iterations the MCMC sampler should iterate 
#' through and discard before beginning to save results
#' @param nmc The number of samples that the MCMC sampler should use to estimate
#' marginal likelihood
#' @param prior_name The name of the prior for which the marginal should be 
#' calculated, this is one of 'Wishart', 'BGL', 'GHS', 'G_Wishart'
#' @param runs The number of complete runs of the graphical evidence method that
#' will be executed. Specifying multiple runs allows estimation of the variance
#' of the estimator and by default will permute the columns of xx such that 
#' each run uses a random column ordering, as marginal likelihood should be 
#' independent of column permutations
#' @param print_progress A boolean which indicates whether progress should be 
#' displayed on the console as each row of the telescoping sum is computed and
#' each run is completed
#' @param permute_columns A boolean which indicates whether columns of xx for
#' runs beyond the first should be randomly permuted to ensure that 
#' marginal calculation is consistent across different column permutations
#' @param alpha A number specifying alpha for the priors of 'Wishart' and 
#' 'G_Wishart'
#' @param lambda A number specifying lambda for the priors of 'BGL' and 'GHS'
#' prior
#' @param V The scale matrix when specifying 'Wishart' or 'G_Wishart' prior
#' @param G The adjacency matrix when specifying 'G_Wishart' prior
#' 
#' 
#' @returns A list of results which contains the mean marginal likelihood, the
#' standard deviation of the estimator, and the raw results in a vector
#' 
#' @examples
#' # Compute the marginal 10 times with random column permutations of xx at each
#' # individual run for G-Wishart prior using 2,000 burnin and 10,000 sampled
#' # values at each call to the MCMC sampler
#' g_params <- gen_params_evidence('G_Wishart')
#' marginal_results <- evidence(
#'   g_params$x_mat, 2e3, 1e4, 'G_Wishart', 3, alpha=2, 
#'   V=g_params$scale_mat, G=g_params$g_mat
#' )
#' @export
evidence <- function(
  xx,
  burnin,
  nmc,
  prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart'),
  runs = 1,
  print_progress = FALSE,
  permute_columns = TRUE,
  alpha = NULL,
  lambda = NULL,
  V = NULL,
  G = NULL
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Detect available cores and assign to compiled program
  num_cores <- parallel::detectCores() - 1
  set_cores(num_cores)
  
  # If dimension is sufficient and prior is Wishart, create cluster for 
  # parallel computation in top level R code
  if (ncol(xx) >= 15 & (prior_name == 'Wishart')) {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  }
  
  # Extract sample size and dimension of xx
  xx <- as.matrix(xx)
  p <- ncol(xx)
  n <- nrow(xx)
  
  # Generate permutations for columns of xx
  permutation_matrix <- t(replicate(runs, sample(p)))
  
  # Store results of repeated runs
  results <- numeric(runs)
  
  for (i in 1:runs) {
    
    # Use a random permutation of columns of xx
    test_perm <- permutation_matrix[i, ]
    
    # Keep xx column order for first run, or if permute_columns is not requested
    if ((i == 1) | !permute_columns) {
      
      xx_perm <- xx
      if (!is.null(G)) {
        G_perm <- G 
      }
      if (!is.null(V)) {
        V_perm <- V 
      }
    }
    
    # Permute xx columns and then calculate S for multiple runs
    else {
      
      xx_perm <- as.matrix(xx[, test_perm])
      if (!is.null(G)) {
        G_perm <- as.matrix(G[test_perm, test_perm])
      }
      if (!is.null(V)) {
        V_perm <- as.matrix(V[test_perm, test_perm])
      }
    }

    # Run graphical evidence method on identified prior
    results[i] <- switch(
      prior_name,
      
      # Largely implemented in C++
      'Wishart' = graphical_evidence_rmatrix(
        xx_perm, burnin, nmc, prior_name, alpha=alpha, V=V_perm,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'BGL' = graphical_evidence_rmatrix(
        xx_perm, burnin, nmc, prior_name, lambda=lambda,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'GHS' = graphical_evidence_rmatrix(
        xx_perm, burnin, nmc, prior_name, lambda=lambda,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'G_Wishart' = graphical_evidence_G_Wishart(
        xx_perm, burnin, nmc, alpha, V_perm, G_perm,
        print_progress=print_progress
      )
    )
    if (print_progress) {
      cat("done with run", i, '\n')
    }
  }
  
  # End parallel cluster
  if (ncol(xx) >= 15 & prior_name == 'Wishart') {
    parallel::stopCluster(cl)
  }
  
  # Filter out infinite and NA values from results
  init_len <- length(results)
  results <- results[!is.na(results) & !is.infinite(results)]
  trunc_len <- length(results)
  message(init_len - trunc_len, " runs excluded for NaN\n")
  
  return(
    list(mean=mean(results), sd=sqrt(var(results)), results=results)
  )
}
