#include "graphical_evidence.h"


/*
 * Use graphical evidence to run an MCMC sampler nmc + burnin times
 * and return the last nmc iterations for G Wishart prior sampling
 */
 // [[Rcpp::export]]
NumericVector prior_sampler_rmatrix(
  int p,
  int burnin,
  int nmc,
  int prior,
  double lambda
) {

  /* Create accumulating storage for mcmc sampling  */
  arma::cube omega_save = arma::cube(p, p, nmc);
  arma::mat omega = arma::mat(p, p, arma::fill::eye);
  arma::mat sigma = arma::inv_sympd(omega);
  arma::mat tau = arma::ones<arma::mat>(p, p);
  arma::mat nu = arma::ones<arma::mat>(p, p);

  /* Initialize indices with ith row excluded for each col  */
  arma::umat ind_noi_mat(p - 1, p);
  initialize_indices(ind_noi_mat);

  /* Allocate working memory in top level scope */
  arma::mat inv_omega_11 = arma::zeros<arma::mat>(p - 1, p - 1);
  arma::mat inv_c = arma::zeros<arma::mat>(p - 1, p - 1);
  arma::vec beta = arma::zeros<arma::vec>(p - 1);

  /* Iterate burnin + nmc times and save results past burnin  */
  arma::uword total_iters = static_cast<arma::uword>(burnin + nmc);

  /* Generate gamma paramters for the MCMC sampling */
  double scale_param = 0;
  if (prior == BGL) {
    scale_param = 2.0 / lambda;
  }
  else if (prior == GHS) {
    scale_param = 2.0 * lambda;
  }

  /* Run MCMC */
  for (arma::uword i = 0; i < total_iters; i++) {
    prior_sample_omega_rmatirx(
      i, burnin, lambda, prior, scale_param, beta, omega,
      inv_omega_11, inv_c, omega_save, ind_noi_mat,
      sigma, tau, nu
    );
  }

  /* Return array of sampled values past burnin */
  return(Rcpp::wrap(omega_save));
}
