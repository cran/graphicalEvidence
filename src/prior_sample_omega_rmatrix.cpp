#include "graphical_evidence.h"


/*
 * Sample Omega for prior sampling of BGL/GHS using 
 * graphical evidence
 */
void prior_sample_omega_rmatirx(
  const int iter,
  const int burnin,
  const double lambda,
  const int prior,
  const double scale_param,
  arma::vec& beta,
  arma::mat& omega,
  arma::mat& inv_omega_11,
  arma::mat& inv_c,
  arma::cube& omega_save,
  arma::umat const& ind_noi_mat,
  arma::mat& sigma,
  arma::mat& tau,
  arma::mat& nu
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = omega.n_rows;

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Calculate this constant that is used in GHS  */
  const double lambda_sq = pow(lambda, 2);

  /* Iterate over each column */
  for (arma::uword i = 0; i < p; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Sample from prior specific gamma density */
    const double gamma_param = g_rgamma.GetSample(1, scale_param);

    /* Inverse of omega excluding row i and col i can be solved in O(n^2) */
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, sigma, p, i);

    /* BGL case */
    if (prior == BGL) {

      inv_c = inv_omega_11 * lambda;
      for (unsigned int j = 0; j < (p - 1); j++) {

        /* Update inv_c diagonal  */
        inv_c.at(j, j) += (1.0 / tau.at(ind_noi[j], i));
      }
    }

    /* GHS case */
    else if (prior == GHS) {

      inv_c = inv_omega_11 * (1 / lambda);
      for (unsigned int j = 0; j < (p - 1); j++) {

        /* Update inv_c diagonal  */
        inv_c.at(j, j) += (1.0 / (tau.at(ind_noi[j], i) * lambda_sq));
      }
    }

    /* Cholesky decomp of inv_c with LAPACK, store result in inv_c  */
    inv_c = arma::chol(inv_c);

    /* Solve chol(inv_c) x = randn(), store result in beta  */
    beta.randn();
    beta = arma::solve(arma::trimatu(inv_c), beta);

    /* Update flex_mem to store inv_omega_11 %*% beta */
    update_omega_inplace(omega, inv_omega_11, beta, ind_noi, gamma_param, i, p);

    /* Update sigma */
    update_sigma_inplace(
      sigma, inv_omega_11, g_vec1, ind_noi, gamma_param, p, i
    );

    /* BGL case */
    if (prior == BGL) {

      /* Calculate a_gig_tau and resuse variable to sample tau_12 */
      for (unsigned int j = 0; j < (p - 1); j++) {

        double gig_tau = pow(beta[j], 2);
        gig_tau = 1 / gigrnd(-1.0 / 2.0, gig_tau, lambda_sq);
        tau.at(ind_noi[j], i) = gig_tau;
        tau.at(i, ind_noi[j]) = gig_tau;
      }
    }

    /* GHS case */
    if (prior == GHS) {

      /* Sample tau_12 and nu_12 */
      for (unsigned int j = 0; j < (p - 1); j++) {

        const double cur_rate = (
          pow(beta[j], 2) / (2.0 * lambda_sq) + (1.0 / nu.at(ind_noi[j], i))
        );

        const double cur_tau = 1.0 / g_rgamma.GetSample(1.0, 1.0 / cur_rate);
        const double cur_nu = 1.0 / g_rgamma.GetSample(1.0, 1.0 / (1.0 + (1.0 / cur_tau)));

        tau.at(ind_noi[j], i) = cur_tau;
        tau.at(i, ind_noi[j]) = cur_tau;
        nu.at(ind_noi[j], i) = cur_nu;
        nu.at(i, ind_noi[j]) = cur_nu;
      }
    }

  }

  /* If iteration is past burnin period, save omega */
  if ((iter - burnin) >= 0) {
    omega_save.slice(iter - burnin) = omega;
  }
}
