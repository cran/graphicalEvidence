#include "graphical_evidence.h"

/*
 * Top level function to calculate MC average of equation 9 in graphical evidence paper
 */

double calc_eq_9(
  arma::uvec const& find_which_ones,
  arma::cube const& arma_inv_sigma_stores,
  arma::mat const& arma_mean_vecs,
  arma::mat const& arma_mean_omega,
  arma::uword const xdim,
  const unsigned int p,
  const unsigned int nmc
 ) {

  /* If no nonzero entries of adj matrix last col, we can immediately return 0  */
  double mc_avg_eq_9 = 0;
  if (xdim) {

    /* Create x data to be evaluated with multivariate normal pdf */
    arma::rowvec x = arma_mean_omega.submat(arma::uvec({p - 1}), find_which_ones);

    /* Calculate log mean of multivariate normal density  */
    mc_avg_eq_9 = looping_dmvnrm_arma(x, arma_mean_vecs, arma_inv_sigma_stores);
  }

  return mc_avg_eq_9;
}
