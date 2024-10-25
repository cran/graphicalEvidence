#include "graphical_evidence.h"

/*
 * Calculates log multivariate normal density for some constant row vector
 * x and varying means and sigmas, stored in mean_vecs and inv_sigma_stores,
 * with the inverse of each index in inv_sigma_stores value taken as sigma.
 * Adapted from RcppArmadillo docs
 */

double looping_dmvnrm_arma(
  arma::rowvec const& x,
  arma::mat const& mean_vecs,
  arma::cube const& inv_sigma_stores
) {

  /* Extract number of MC iterations n  */
  arma::uword const n = mean_vecs.n_cols;

  /* Dimension of x data*/
  arma::uword const xdim = x.n_cols;

  /* Armadillo vector to return */
  arma::vec out(n);
  
  /* Armadillo vars extracted each iteration  */
  arma::rowvec z;
  arma::mat sigma;
  arma::mat rooti;

  /* PDF constants  */
  double const constants = -((double) xdim) / 2.0 * LOG2PI;
    
  /* Parallelization of non sequentail calcluations if high dimension */
  if ((xdim >= PARALLEL_XDIM_THRESHHOLD) && (n >= MCMC_LEN_THRESHHOLD)){
    #pragma omp parallel for schedule(static) private(z, sigma, rooti)
    for (arma::uword i = 0; i < n; i++) {
      looping_mvpdf_process_iteration(
        out, x, inv_sigma_stores, mean_vecs, rooti, z, sigma, constants, i
      );
    }
  }
  /* If xdim is below threshold for parallelizaiton, use sequential */
  else {
    for (arma::uword i = 0; i < n; i++) {
      looping_mvpdf_process_iteration(
        out, x, inv_sigma_stores, mean_vecs, rooti, z, sigma, constants, i
      );
    }
  }

  return log(arma::mean(exp(out)));
}


/*
 * One loop iteration of extracting inverse sigmas and mean vectors
 * to calculate multivariate normal PDF for some constant x value
 */

void looping_mvpdf_process_iteration(
  arma::vec& out,
  arma::rowvec const& x,
  arma::cube const& inv_sigma_stores,
  arma::mat const& mean_vecs,
  arma::mat& rooti,
  arma::rowvec& z,
  arma::mat& sigma,
  double const constants,
  arma::uword i
) {

  /* Extract relevant inv sigma and derive sigma  */
  sigma = arma::inv_sympd(inv_sigma_stores.slice(i), arma::inv_opts::allow_approx);

  /* Multivariate normal PDF  */

  /* Calculation of the inverse of the square root of a positive definite matrix */
  rooti = arma::inv(arma::trimatu(arma::chol(sigma)));

  /* Constants needed for multivariate normal PDF */
  double const rootisum = arma::sum(log(rooti.diag()));
  double const other_terms = rootisum + constants;

  /* PDF given input x for currently considered mean  */
  z = (x - mean_vecs.col(i).t());
  inplace_tri_mat_mult(z, rooti);
  out[i] = other_terms - 0.5 * arma::dot(z, z);

}


/*
 * Calculates matrix product of upper triangular matrix and row vector x
 * in place for input row vec x.
 * Sourced from RcppArmadillo docs
 */

void inplace_tri_mat_mult(
  arma::rowvec& x, 
  arma::mat const& trimat
) {
  arma::uword const n = trimat.n_cols;

  for (unsigned j = n; j-- > 0;) {
    double tmp(0.);
    for (unsigned i = 0; i <= j; ++i) {
      tmp += trimat.at(i, j) * x[i];
    }
    x[j] = tmp;
  }
}
