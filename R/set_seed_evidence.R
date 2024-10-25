#' @title Set the Random Seed
#'
#' @description
#' Sets the random seed of both the R session (using set.seed) and the compiled
#' sampler, as both samplers are used during any calls to evidence(...) or 
#' prior_sampling(...)
#' 
#' @param seed a random seed that will be passed to the interpreted random 
#' number generator using set.seed, and will be passed to the compiled random
#' number generator using private Rcpp package function set_seed
#' 
#' @returns No return value, this function is called to set the random seed.
#' @examples
#' set_seed_evidence(42)
set_seed_evidence <- function(
    seed
) {
  
  # Set interpreted seed
  set.seed(seed)
  
  # Set compiled seed
  set_seed(seed)
}