#include "graphical_evidence.h"


/*
 * Use Hao Wang decomposition to run an MCMC sampler nmc + burnin times
 * and accumulate omega_reduced, gamma_subtractors, and posterior mean of
 * omega_22 and then use those values to calculate eq. 11 in paper, this
 * function is called for BGL, GHS, and Wishart priors
 */
 // [[Rcpp::export]]
List mcmc_last_col_rmatrix(
  const unsigned int n,
  const unsigned int burnin,
  const unsigned int nmc,
  const unsigned int p,
  const double dof,
  const double lambda,
  const int prior,
  NumericVector fixed_last_col_nvec,
  NumericVector s_mat_nvec,
  NumericVector tau_mat_nvec,
  NumericVector gibbs_mat_nvec,
  NumericVector post_mean_omega_nvec
) {

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::vec fixed_last_col(fixed_last_col_nvec.begin(), p - 1);
  arma::mat s_mat(s_mat_nvec.begin(), p, p);
  arma::mat omega(post_mean_omega_nvec.begin(), p, p);

  /* Initialize prior specific parameters */
  arma::mat nu;
  arma::mat tau;
  arma::mat gibbs_mat;

  /* Sampler with last column fixed relies on p_reduced sampling  */
  const int p_reduced = p - 1;

  /* Calculate inital inverse of omega  */
  arma::mat sigma = arma::inv_sympd(omega);

  /* Use global memory for variables used to update sigma/omega */
  arma::mat inv_omega_11_full = arma::mat(g_mat1, p - 1, p - 1, false, true);

  /* Set gamma parameters and calculation memory dependent on prior,  */
  /* and reduce matrices to p_reduced dim where needed                */  
  double scale_params[MAX_DIM];
  double shape_param = 0;
  arma::mat inv_c = arma::zeros(p_reduced - 1, p_reduced - 1);
  arma::mat inv_omega_11 = arma::zeros(p_reduced - 1, p_reduced - 1);
  arma::vec beta = arma::zeros(p_reduced - 1);
  omega = omega.submat(0, 0, p_reduced - 1, p_reduced - 1);

  /* Wishart case */
  if (prior == WISHART) {

    /* Gamma parameters */
    shape_param = (dof + (double)n - p + 1) / 2;
    for (unsigned int i = 0; i < p; i++) {
      scale_params[i] = 2.0 / (s_mat.at(i, i) + 1);
    }
  }

  /* BGL case */
  else if (prior == BGL) {

    /* Gamma parameters */
    shape_param = ((double)n / 2) + 1;
    for (unsigned int i = 0; i < p; i++) {
      scale_params[i] = 2.0 / (s_mat.at(i, i) + lambda);
    }

    /* Deep copy matrices */
    gibbs_mat = arma::mat(gibbs_mat_nvec.begin(), p, p);
    tau = arma::mat(tau_mat_nvec.begin(), p, p);
    
    /* Reduce matrices  */
    gibbs_mat = gibbs_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
    tau = tau.submat(0, 0, p_reduced - 1, p_reduced - 1);
  }

  /* GHS case */
  else if (prior == GHS) {

    /* Gamma parameters */
    shape_param = ((double)n / 2) + 1;
    for (unsigned int i = 0; i < p; i++) {
      scale_params[i] = 2.0 / (s_mat.at(i, i) + (1.0 / lambda));
    }

    /* Deep copy matrices */
    gibbs_mat = arma::mat(gibbs_mat_nvec.begin(), p, p);
    tau = arma::mat(tau_mat_nvec.begin(), p, p);
    nu = arma::ones(p, p);

    /* Reduce matrices  */
    gibbs_mat = gibbs_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
    tau = tau.submat(0, 0, p_reduced - 1, p_reduced - 1);
    nu = nu.submat(0, 0, p_reduced - 1, p_reduced - 1);
  }

  /* Reduced s_mat after extracting relevant gamma parameters */
  s_mat = s_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);

  /* Calculate outer product of last col  */
  arma::mat last_col_outer = fixed_last_col * fixed_last_col.t();

  /* Initialize indices which exclude i */
  arma::umat ind_noi_mat(p_reduced - 1, p_reduced);
  initialize_indices(ind_noi_mat);

  /* Initialize accumulator memory  */
  arma::mat omega_reduced_acc = arma::zeros(p_reduced, p_reduced);
  arma::vec gamma_subtractors = arma::zeros(nmc);
  double omega_22_acc = 0;

  /* Perform main loop of mcmc sampling */
  for (arma::uword i = 0; i < (burnin + nmc); i++) {

    /* Get random gamma value and calculate omega_22  */
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[p - 1]);

    /* Efficient in place calculation of inv_omega_11_full */
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
      sample_omega_last_col_rmatrix(
        p_reduced, shape_param, scale_params, omega_22, lambda, dof, prior,
        beta, omega, inv_c, inv_omega_11, sigma, tau, nu, ind_noi_mat, s_mat,
        gibbs_mat, last_col_outer
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

  double mc_avg_eq_11 = calc_eq_11(
    omega_22_acc, shape_param, scale_params[p - 1], nmc, gamma_subtractors
  );

  List z = List::create(
    mc_avg_eq_11, Rcpp::wrap(omega_reduced_acc), omega_22_acc
  );

  return z;
}
