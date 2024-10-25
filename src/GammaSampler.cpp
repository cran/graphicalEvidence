#include "graphical_evidence.h"

/* Initialize the global gamma sampler  */
GammaSampler g_rgamma;


/*
 * Returns the density of the gamma distribution under the specified 
 * shape and scale parameters for input x
 */

double gamma_density(const double x, const double shape, const double scale) {
  double log_density = (
    (shape - 1) * log(x) - (x / scale) - shape * log(scale) - lgamma(shape)
  );
  return exp(log_density);
}


/*
 * Returns n samples from Gamma distribution with input rates of vec_rates 
 * and input shapes of vec_shapes
 */
// [[Rcpp::export]]
NumericVector rgamma_compiled(int n, NumericVector vec_shapes, NumericVector vec_rates) {
  
  /* Sampled values to be returned */
  NumericVector rval(n);

  for (int i = 0; i < n; i++) {
    rval(i) = g_rgamma.GetSample(vec_shapes[i % vec_shapes.size()], 1 / vec_rates[i % vec_rates.size()]);
  }

  return(rval);
}
