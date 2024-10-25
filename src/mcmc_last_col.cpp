#include "graphical_evidence.h"


/*
 * Use Hao Wang decomposition to run an MCMC sampler nmc + burnin times
 * and accumulate omega_reduced, gamma_subtractors, and posterior mean of
 * omega_22 and then use those values to calculate eq. 11 in paper,
 * this function is called for G_Wishart prior only
 */
 // [[Rcpp::export]]
List mcmc_last_col(
  const unsigned int n,
  const unsigned int burnin,
  const unsigned int nmc,
  const unsigned int p,
  const double alpha,
  NumericVector fixed_last_col_nvec,
  NumericVector s_mat_nvec,
  NumericVector scale_mat_nvec,
  NumericVector g_mat_adj_nvec,
  NumericVector gibbs_mat_nvec,
  NumericVector post_mean_omega_nvec
) {

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::vec fixed_last_col(fixed_last_col_nvec.begin(), p - 1);
  arma::mat s_mat(s_mat_nvec.begin(), p, p);
  arma::mat scale_mat(scale_mat_nvec.begin(), p, p);
  arma::mat omega(post_mean_omega_nvec.begin(), p, p);
  arma::mat gibbs_mat(gibbs_mat_nvec.begin(), p, p);
  arma::mat g_mat_adj(g_mat_adj_nvec.begin(), p, p);

  /* Sampler with last column fixed relies on p_reduced sampling  */
  const int p_reduced = p - 1;

  /* Set omega_22 gamma parameters  */
  double scale_params[MAX_DIM];
  const double shape_param = alpha + ((double) n / 2) + 1;
  for (unsigned int i = 0; i < p; i++) {
    scale_params[i] = 2.0 / (s_mat.at(i, i) + scale_mat.at(i, i));
  }

  /* Calculate outer product of last col  */
  arma::mat last_col_outer = fixed_last_col * fixed_last_col.t();

  /* Calculate inital inverse of omega  */
  arma::mat sigma = arma::inv_sympd(omega);

  /* Use global memory for variables used to update sigma/omega */
  arma::mat inv_omega_11_full = arma::mat(g_mat1, p - 1, p - 1, false, true);

  /* Reduce matrices */
  omega = omega.submat(0, 0, p_reduced - 1, p_reduced - 1);
  gibbs_mat = gibbs_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
  g_mat_adj = g_mat_adj.submat(0, 0, p_reduced - 1, p_reduced - 1);
  scale_mat = scale_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
  s_mat = s_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);

  /* Identify indices where ones and zeros exist for adjacency matrix */
  arma::umat ind_noi_mat(p_reduced - 1, p_reduced);
  std::vector<arma::uvec> find_which_ones(p_reduced);
  std::vector<arma::uvec> find_which_zeros(p_reduced);
  initialize_indices(
    g_mat_adj, ind_noi_mat, find_which_ones, find_which_zeros
  );

  /* Initialize accumulators  */
  arma::mat omega_reduced_acc = arma::zeros(p_reduced, p_reduced);
  arma::vec gamma_subtractors = arma::zeros(nmc);
  double omega_22_acc = 0;

  /* Initialize calculation memory  */
  arma::mat inv_c = arma::zeros(p_reduced - 1, p_reduced - 1);
  arma::mat inv_omega_11 = arma::zeros(p_reduced - 1, p_reduced - 1);
  arma::vec beta = arma::zeros(p_reduced - 1);

  /* Perform main loop of mcmc sampling */
  for (arma::uword i = 0; i < (burnin + nmc); i++) {

    /* Get random gamma value and calculate omega_22  */
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[p - 1]);

    /* Efficient calculation of inv_omega_11_full */
    last_col_calc_inv_omega_11_full(inv_omega_11_full, sigma);

    /* Calculate gamma_subtractor and omega_22  */
    double gamma_subtractor = calc_gamma_subtractor(
      fixed_last_col, inv_omega_11_full
    );
    double omega_22 = gamma_subtractor + gamma_sample;

    /* Prepare sigma_reduced for iteration of 1 to p_reduced sampling */
    /* of omega_reduced, each requiring calculation of inv_omega_11   */
    last_col_prepare_sigma_reduced(
      sigma, inv_omega_11_full, fixed_last_col, gamma_sample
    );

    /* Save gamma_subtractor now if burnin is past, this avoids recalculation   */
    /* of gamma_subtractor in calculation of eq. 11, although it shifts the     */
    /* the comparison of gamma_subtractors back by one in the iteration process */
    /* when comparing against posterior mean of omega_22, but this is fine      */
    if (i >= burnin) {
      gamma_subtractors[i - burnin] = gamma_subtractor;
    }

    /* In general case call omega reduced sampler */
    if (p_reduced > 1) {
      sample_omega_last_col(
        p_reduced, shape_param, scale_params, omega_22, beta, omega, inv_c,
        inv_omega_11, sigma, find_which_ones, find_which_zeros, ind_noi_mat,
        s_mat, scale_mat, gibbs_mat, last_col_outer
      );
    }

    /* In special case where the matrices considered are single elements  */
    /* the sampling process is highly simplified                          */
    else {
      gamma_sample = g_rgamma.GetSample(shape_param, scale_params[0]);
      omega.at(0, 0) = gamma_sample + (
        fixed_last_col[0] * fixed_last_col[0] / omega_22
      );
    }

    /* Update sigma's last column and row for next iteration  */
    update_sigma_last_col(sigma, fixed_last_col, omega_22);

    /* Save results if burnin is complete */
    if (i >= burnin) {
      omega_22_acc += omega_22;
      omega_reduced_acc += omega;
    }
  }

  /* Get posterior means of sampled values */
  omega_22_acc /= nmc;
  omega_reduced_acc /= nmc;

  /* Calculate equation 11 in paper */
  double mc_avg_eq_11 = calc_eq_11(
    omega_22_acc, shape_param, scale_params[p - 1], nmc, gamma_subtractors
  );

  /* Return values to R code  */
  List z = List::create(
    mc_avg_eq_11, Rcpp::wrap(omega_reduced_acc), omega_22_acc
  );

  return z;

}


/* 
 * Calculates the full p - 1 by p - 1 inv_omega_11_full in place using
 * O(n^2) operations, given the return memory and sigma (inverse of omega).
 * This function is called once prior to iterating through 1 to p_reduced
 * to sample each column of omega
 */

void last_col_calc_inv_omega_11_full(
  arma::mat& inv_omega_11_full,
  arma::mat const& sigma
) {

  unsigned int p = sigma.n_rows;

  // arma::cout << "initial sigma: \n" << sigma << arma::endl;

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      inv_omega_11_full.at(j, i) = sigma.at(j, i) - (
        sigma.at(j, p - 1) * sigma.at(i, p - 1) /
        sigma.at(p - 1, p - 1)
      );
    }
  }
}


/*
 * Calculate gamma_subtractor using global memory with formula
 * gamma_subtractor = t(fixed_last_col) %*% inv_omega_11 %*% fixed_last_col
 */

double calc_gamma_subtractor(
  arma::vec const& fixed_last_col,
  arma::mat const& inv_omega_11_full
) {

  const unsigned int p_reduced = fixed_last_col.n_elem;
  double gamma_subtractor = 0;
  for (unsigned int j = 0; j < p_reduced; j++) {

    /* Store fixed_last_col.t() * inv_omega_11_full in g_vec1  */
    g_vec1[j] = 0;
    for (unsigned int k = 0; k < p_reduced; k++) {

      /* First fixed_last_col.t() * inv_omega_11_full[, k] */
      g_vec1[j] += (fixed_last_col[k] * inv_omega_11_full.at(k, j));
    }

    /* Accumulate fixed_last_col_omega[j] * fixed_last_col[j] */
    gamma_subtractor += (fixed_last_col[j] * g_vec1[j]);
  }

  return gamma_subtractor;
}


/*
 * Prepares reduced sigma such that at each iteration of the last col
 * fixed sampler, inv_omega_11 can be calculated from sigma_reduced
 * (where sigma_reduced is just sigma from 0 to p_reduced) according to
 * x = inv_omega_11_full %*% fixed_last col
 * sigma_reduced = inv_omega_11_full + x %*% t(x) / gamma_param
 */

void last_col_prepare_sigma_reduced(
  arma::mat& sigma,
  arma::mat const& inv_omega_11_full,
  arma::vec const& fixed_last_col,
  const double gamma_param
) {

  unsigned int p_reduced = fixed_last_col.n_elem;

  /* Calculate inv_omega_11_full %*% fixed_last_col */
  for (unsigned int i = 0; i < p_reduced; i++) {
    double sum = 0.0;
    for (unsigned int j = 0; j < p_reduced; j++) {
      sum += (inv_omega_11_full.at(i, j) * fixed_last_col[j]);
    }
    g_vec1[i] = sum;
  }

  /* Update sigma from indices 0 to p_reduced */

  /* Iterate over cols  */
  for (unsigned int i = 0; i < p_reduced; i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < p_reduced; j++) {

      sigma.at(j, i) = inv_omega_11_full.at(j, i) + (
        g_vec1[j] * g_vec1[i] / gamma_param
      );
    }
  }
}


/* 
 * Updates the last column and row of sigma after completing sampling 
 * iterations from 1 to p_reduced in the fixed last col sampler, 
 * note that sigma_reduced is already updated in the sampling process
 */

void update_sigma_last_col(
  arma::mat& sigma,
  arma::vec const& fixed_last_col,
  const double omega_pp
) {

  const unsigned int p = sigma.n_rows;

  /* Initialize sigma_22 =                       */
  /* 1 / omega_pp + (1 / omega_pp)^2 *           */
  /* last_col.t() %*% sigma_reduced %*% last_col */
  sigma.at(p - 1, p - 1) = 1 / omega_pp;

  /* Update sigma_*/
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Update sigma_12 = -sigma_reduced %*% last_col / omega_pp in place  */
    sigma.at(i, p - 1) = 0;
    for (unsigned int j = 0; j < (p - 1); j++) {
      sigma.at(i, p - 1) += (sigma.at(i, j) * fixed_last_col[j] / omega_pp);
    }

    /* Accumulate sigma_22  */
    sigma.at(p - 1, p - 1) += (sigma.at(i, p - 1) * fixed_last_col[i] / omega_pp);

    /* Update sigma_12 and sigma_21 */
    sigma.at(i, p - 1) *= -1;
    sigma.at(p - 1, i) = sigma.at(i, p - 1);
  }
}
