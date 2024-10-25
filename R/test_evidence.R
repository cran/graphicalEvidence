#' @title Test Graphical Evidence
#' 
#' @description
#' Tests any of the allowed priors with preexisting test inputs which
#' should yield reproducible results, as the random seed is reset during
#' this function call
#' 
#' @param num_runs An integer number that specifies how many runs of 
#' graphical evidence will be executed on the test parameters, as multiple
#' runs allows us to quantify uncertainty on the estimator.
#' 
#' @param prior_name The name of the prior for being tested with preexisting
#' test parameters, this is one of 'Wishart', 'BGL', 'GHS', 'G_Wishart'
#' 
#' @param input_seed An integer value that will be used as a random seed to
#' make outputs repeatable.
#' 
#' @returns A list of results which contains the mean marginal likelihood, the
#' standard deviation of the estimator, and the raw results in a vector
#' 
#' @examples
#' # Compute the marginal 10 times with random column permutations of the 
#' # preexisting test parameters for G-Wishart prior 
#' test_evidence(num_runs=3, 'G_Wishart')
#' @export
test_evidence <- function(
    num_runs, 
    prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart'),
    input_seed = NULL
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Get predetermined test inputs
  params <- gen_params_evidence(prior_name)
  
  # Make results repeatable if input seed
  if (!is.null(input_seed)) {
    set_seed_evidence(input_seed)
  }
  
  # Save start in R program time
  our_time <- proc.time()

  # Call top level function
  estimated_marginal_store <- switch(
    prior_name,
    
    'Wishart' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, alpha=7, V=params$scale_mat,
    ),
    
    'BGL' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, lambda=1
    ),
    
    'GHS' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, lambda=1
    ),
    
    'G_Wishart' = evidence(
      xx=params$x_mat, burnin=2e3, nmc=1e4, prior_name=prior_name, 
      runs=num_runs, alpha=2, V=params$scale_mat, G=params$g_mat
    )
  )
  
  # Cumulative execution time in R program
  message('Execution time in R program per run (seconds) added to results...\n')
  estimated_marginal_store$time <- (proc.time() - our_time) / num_runs

  message("Test params added to results object... \n")
  estimated_marginal_store$test_params <- params
  
  # Visually plot multiple runs in a histogram
  if (num_runs > 1) {
    hist(
      estimated_marginal_store$results, xlab='estimated marginal', 
      main='Histogram of Results', breaks=40
    )
  }
  return(estimated_marginal_store)
}