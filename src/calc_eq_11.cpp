#include "graphical_evidence.h"

/*
 * Calculates equation 11 in paper using sampled results for omega_reduced
 * and the elements of omega at column and row p throughout the MCMC sampling.
 * This can be parallelized with OpenMP as the calculation is not inherently
 * sequential.
 */

double calc_eq_11(
  const double omega_22_mean,
  const double shape_param,
  const double scale_param,
  const unsigned int nmc,
  arma::vec const& gamma_subtractors
) {

  /* Iterate through accumulated gamma subtractors and calculate gamma density  */
  double gamma_acc = 0.0;
  for (arma::uword i = 0; i < nmc; i++) {

    /* Gamma value to check is mean pth row/col of omega minus saved subtractors  */
    double temp_gamma = omega_22_mean - gamma_subtractors[i];

    /* Check if temp_gamma is not positive  */
    if (temp_gamma > 0) {

      /* Accumulate the gamma density */
      gamma_acc += gamma_density(temp_gamma, shape_param, scale_param);
    }
  }
  gamma_acc /= nmc;

  return log(gamma_acc);
}
