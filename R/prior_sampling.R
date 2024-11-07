#' @title Sample The Precision Matrix
#' 
#' @description
#' Takes specified prior_name and relevant parameters to sample the precision
#' matrix nmc times after discarding the first number of runs specified by 
#' burnin. 
#' 
#' @param p The dimension of the precision matrix that will be sampled
#' @param burnin The number of iterations the MCMC sampler should iterate 
#' through and discard before beginning to save results
#' @param nmc The number of samples that will be drawn
#' @param prior_name The name of the prior for which the marginal should be 
#' calculated, this is one of 'Wishart', 'BGL', 'GHS', 'G_Wishart'
#' @param alpha A number specifying alpha for the priors of 'Wishart' and 
#' 'G_Wishart'
#' @param lambda A number specifying lambda for the priors of 'BGL' and 'GHS'
#' prior
#' @param V The scale matrix when specifying 'Wishart' or 'G_Wishart' prior
#' @param G The adjacency matrix when specifying 'G_Wishart' prior
#' 
#' @returns An array of dim nmc x p x p where each p x p slice is one sample of 
#' the precision matrix
#' 
#' @examples
#' # Draw 5000 samples of the precision matrix for GHS prior distribution with
#' # parameter lambda set to 1
#' ghs_samples <- prior_sampling(5, 1e3, 5e3, 'GHS', lambda=1)
#' @export
prior_sampling <- function(
    p, 
    burnin, 
    nmc,
    prior_name = c('BGL', 'GHS', 'G_Wishart'),
    G = NULL, 
    V = NULL, 
    alpha = NULL,
    lambda = NULL
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)

  coded_prior <- switch(
    prior_name,
    'BGL' = 1,
    'GHS'= 2
  )

  # Draw samples using compiled sampler dependent on prior
  sampled_omegas <- switch(
    prior_name,
    
    # Implemented in C++
    'BGL' = prior_sampler_rmatrix(
      p, burnin, nmc, coded_prior, lambda
    ),
    
    # Implemented in C++
    'GHS' = prior_sampler_rmatrix(
      p, burnin, nmc, coded_prior, lambda
    ),
    
    # Implemented in C++
    'G_Wishart' = prior_sampler_G_Wishart(
      p, burnin, nmc, as.matrix(G), 
      as.matrix(V), alpha
    )
  )
  
  return(sampled_omegas)
}