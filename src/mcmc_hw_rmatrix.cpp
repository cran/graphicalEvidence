#include "graphical_evidence.h"


/*
 * Use Hao Wang decomposition to run an MCMC sampler nmc + burnin times
 * and accumulate omega, mean_vec, and inv_c_required to calculate
 * mc average equation 9 after mcmc sampling is complete for
 * the G-Wishart prior.
 */
// [[Rcpp::export]]
List mcmc_hw_rmatrix(
  int n,
  int burnin,
  int nmc,
  int p,
  int prior,
  double dof,
  double lambda,
  NumericVector s_mat_nvec,
  NumericVector gibbs_mat_nvec
) {

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::mat s_mat(s_mat_nvec.begin(), p, p);

  /* Create accumulating storage for mcmc sampling  */
  arma::mat mean_vec_store = arma::mat(p - 1, nmc);
  arma::cube inv_c_required_store = arma::cube(p - 1, p - 1, nmc);
  arma::mat omega_save = arma::zeros<arma::mat>(p, p);
  arma::mat tau_save = arma::zeros<arma::mat>(p, p);

  /* Only BGL and GHS require gibbs accumulator matrix  */
  arma::mat gibbs_mat;
  if (prior == BGL || prior == GHS) {
    gibbs_mat = arma::mat(gibbs_mat_nvec.begin(), p, p);
  }

  /* Initialize indices with ith row excluded for each col  */
  arma::umat ind_noi_mat(p - 1, p);
  initialize_indices(ind_noi_mat);

  /* Create variables updated through mcmc interation */
  arma::mat omega = arma::eye(p, p);
  arma::mat cur_sigma = arma::inv_sympd(omega);
  arma::mat tau;
  arma::mat nu;

  /* Allocate working memory in top level scope */
  arma::mat inv_omega_11 = arma::zeros<arma::mat>(p - 1, p - 1);
  arma::mat inv_c = arma::zeros<arma::mat>(p - 1, p - 1);
  arma::vec beta = arma::zeros<arma::vec>(p - 1);

  /* Allocate additional iterating variables dependent on prior */
  if (prior == BGL || prior == GHS) {

    /* Tau required for BGL and GHS */
    tau = arma::ones(p, p);

    if (prior == GHS) {

      /* Nu required only for GHS */
      nu = arma::ones(p, p);
    }
  }

  /* Iterate burnin + nmc times and save results past burnin  */
  arma::uword total_iters = static_cast<arma::uword>(burnin + nmc);

  /* Generate gamma paramters for the MCMC sampling */
  double shape = 0;
  double scale_vec[MAX_DIM] = {0};
  get_gamma_params_hw_rmatrix(&shape, scale_vec, prior, dof, n, lambda, s_mat);

  /* Run MCMC */
  for (arma::uword i = 0; i < total_iters; i++) {
    sample_omega_hw_rmatrix(
      i, burnin, prior, dof, lambda, beta, omega, cur_sigma, tau, nu,
      inv_omega_11, inv_c, omega_save, tau_save, mean_vec_store,
      inv_c_required_store, gibbs_mat, ind_noi_mat, s_mat, shape, scale_vec
    );
  }

  /* Get posterior mean of sampled omega  */
  omega_save /= nmc;
  tau_save /= nmc;

  /* Calcluate mc average equation 9 and return post mean omega as well */
  arma::uvec which_ones = ind_noi_mat.unsafe_col(p - 1);
  double mc_avg_eq_9 = calc_eq_9(
    which_ones, inv_c_required_store, mean_vec_store,
    omega_save, p - 1, p, nmc
  );

  List z = List::create(Rcpp::wrap(omega_save), Rcpp::wrap(tau_save), mc_avg_eq_9);

  return z;
}


/*
 * Generate the gamma parameters for Hao Wang unrestricted sampler dependent on
 * the prior specified, these are unchanging for the entire MCMC process
 */

void get_gamma_params_hw_rmatrix(
  double* shape,
  double* scale_vec,
  const int prior,
  const int dof,
  const int n,
  const double lambda,
  arma::mat const& s_mat
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  if (prior == WISHART) {
    *shape = ((double) dof + n - p + 1) / 2;
    for (unsigned int i = 0; i < p; i++) {
      scale_vec[i] = 2 / (s_mat.at(i, i) + 1);
    }
  }
  else if (prior == BGL) {
    *shape = ((double) n / 2) + 1;
    for (unsigned int i = 0; i < p; i++) {
      scale_vec[i] = 2 / (s_mat.at(i, i) + lambda);
    }
  }
  else if (prior == GHS) {
    *shape = ((double) n / 2) + 1;
    for (unsigned int i = 0; i < p; i++) {
      scale_vec[i] = 2 / (s_mat.at(i, i) + (1 / lambda));
    }
  }
}
