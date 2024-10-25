#' logmvgamma
#'
#' Compute logarithm multivariate Gamma function.
#' Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
#' log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)
#' Written by Michael Chen (sth4nth@gmail.com).
#' 
#' @param x data passed to the multivariate gamma function
#' @param d the dimension of the data 
#' 
#' @returns The log multivariate gamma function results
#' @keywords internal
#' @noRd
logmvgamma <- function(x, d) {
  s <- dim(as.matrix(x))
  x <- array(x, c(1, s[1] * s[2]))
  repeated_x <- matrix(rep(x, each=d), ncol=length(x), byrow=TRUE)
  vector_to_add <- (1 - (1:d)) / 2
  x <- repeated_x + vector_to_add
  y <- d * (d - 1) / 4 * log(pi) + sum(lgamma(x))
  return(array(y, s))
}