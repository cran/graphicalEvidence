#include "graphical_evidence.h"


/*
 * Use graphicalEvidence to run an MCMC sampler nmc + burnin times
 * and return the last nmc iterations for G Wishart prior sampling
 */
 // [[Rcpp::export]]
NumericVector prior_sampler_G_Wishart(
  int p,
  int burnin,
  int nmc,
  NumericVector g_mat_adj_nvec,
  NumericVector scale_mat_nvec,
  int alpha
) {

  /* Reset global memory  */
  memset(g_vec1, 0, sizeof(double) * MAX_DIM);
  memset(g_vec2, 0, sizeof(double) * MAX_DIM);
  memset(g_mat1, 0, sizeof(double) * MAX_DIM * MAX_DIM);

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::mat scale_mat(scale_mat_nvec.begin(), p, p);
  arma::mat g_mat_adj(g_mat_adj_nvec.begin(), p, p);

  /* Initialize omega */
  arma::mat omega = arma::diagmat(arma::ones(p));
  for (arma::uword i = 0; i < static_cast<arma::uword>(p); i++) {
    for (arma::uword j = 0; j < static_cast<arma::uword>(p); j++) {
      if (g_mat_adj.at(i, j) && (i != j)) {
        omega.at(i, j) = 0.01;
      }
    }
  }

  /* Initialize sigma */
  arma::mat sigma = arma::inv_sympd(omega);

  /* Identify indices where ones and zeros exist for adjacency matrix */
  arma::umat ind_noi_mat(p - 1, p);
  std::vector<arma::uvec> find_which_ones(p);
  std::vector<arma::uvec> find_which_zeros(p);
  initialize_indices(
    g_mat_adj, ind_noi_mat, find_which_ones, find_which_zeros
  );

  /* Create storage for mcmc sampling  */
  arma::cube omega_save = arma::zeros<arma::cube>(p, p, nmc);

  /* Allocate working memory in top level scope */
  arma::mat inv_omega_11 = arma::zeros(p - 1, p - 1);
  arma::mat inv_c = arma::zeros(p - 1, p - 1);
  arma::vec beta = arma::zeros(p - 1);

  /* Iterate burnin + nmc times and save results past burnin  */
  arma::uword total_iters = static_cast<arma::uword>(burnin + nmc);

  /* Iterate and update omega */
  for (arma::uword i = 0; i < total_iters; i++) {
    prior_sample_omega(
      i, burnin, alpha, beta, omega, inv_omega_11, inv_c,
      omega_save, g_mat_adj, ind_noi_mat, find_which_ones,
      find_which_zeros, scale_mat, sigma
    );
  }

  /* Return array of sampled values past burnin */
  return(Rcpp::wrap(omega_save));
}
