#include "graphical_evidence.h"


/*
 * Use Hao Wang decomposition to run an MCMC sampler nmc + burnin times
 * and accumulate omega, mean_vec, and inv_c_required to calculate
 * mc average equation 9 after mcmc sampling is complete for
 * the G-Wishart prior.
 */
// [[Rcpp::export]]
List mcmc_hw(
  int n,
  int burnin,
  int nmc,
  double alpha,
  int p,
  NumericVector s_mat_nvec,
  NumericVector scale_mat_nvec,
  NumericVector g_mat_adj_nvec,
  NumericVector gibbs_mat_nvec,
  NumericVector init_gibbs_nvec
) {

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::mat s_mat(s_mat_nvec.begin(), p, p);
  arma::mat scale_mat(scale_mat_nvec.begin(), p, p);
  arma::mat omega(init_gibbs_nvec.begin(), p, p);
  arma::mat gibbs_mat(gibbs_mat_nvec.begin(), p, p);
  arma::mat g_mat_adj(g_mat_adj_nvec.begin(), p, p);

  arma::mat sigma = arma::inv_sympd(omega);

  /* Identify indices where ones and zeros exist for adjacency matrix */
  arma::umat ind_noi_mat(p - 1, p);
  std::vector<arma::uvec> find_which_ones(p);
  std::vector<arma::uvec> find_which_zeros(p);
  initialize_indices(
    g_mat_adj, ind_noi_mat, find_which_ones, find_which_zeros
  );
  
  /* Create accumulating storage for mcmc sampling  */
  arma::uword xdim = find_which_ones[p - 1].n_elem;
  arma::mat mean_vec_store = arma::zeros<arma::mat>(xdim, nmc);
  arma::cube inv_c_required_store = arma::zeros<arma::cube>(xdim, xdim, nmc);
  arma::mat omega_save = arma::zeros<arma::mat>(p, p);

  /* Allocate working memory in top level scope */
  arma::mat inv_omega_11 = arma::zeros(p - 1, p - 1);
  arma::mat inv_c = arma::zeros(p - 1, p - 1);
  arma::vec beta = arma::zeros(p - 1);

  /* Iterate burnin + nmc times and save results past burnin  */
  arma::uword total_iters = static_cast<arma::uword>(burnin + nmc);

  /* Generate gamma sampling parameters */
  double scale_params[MAX_DIM];
  const double shape_param = alpha + ((double)n / 2) + 1;
  for (unsigned int i = 0; i < (unsigned int) p; i++) {
    scale_params[i] = 2 / (s_mat.at(i, i) + scale_mat.at(i, i));
  }

  /* Iterate and update omega */
  for (arma::uword i = 0; i < total_iters; i++) {
    sample_omega_hw(
      i, burnin, alpha, beta, omega, inv_omega_11, inv_c, omega_save,
      mean_vec_store, inv_c_required_store, gibbs_mat, g_mat_adj, ind_noi_mat,
      find_which_ones, find_which_zeros, scale_mat, s_mat, sigma, shape_param,
      scale_params
    );
  }

  /* Get posterior mean of sampled omega  */
  omega_save /= nmc;

  /* Calcluate mc average equation 9 in associated paper  */
  double mc_avg_eq_9 = calc_eq_9(
    find_which_ones[p - 1], inv_c_required_store, mean_vec_store,
    omega_save, xdim, p, nmc
  );

  /* Return values to R code  */
  List z = List::create(Rcpp::wrap(omega_save), mc_avg_eq_9);

  return z;
}
