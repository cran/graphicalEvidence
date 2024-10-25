#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates omega accumulator, inv_c accumulator, and mean_vec accumulator as
 * determined by iterations exceeding burnin period for Wishart, BGL, and
 * GHS prior, with conditional branching determined by the input prior.
 */

void sample_omega_hw_rmatrix(
  const int iter,
  const int burnin,
  const int prior,
  const int dof,
  const double lambda,
  arma::vec& beta,
  arma::mat& omega,
  arma::mat& cur_sigma,
  arma::mat& tau,
  arma::mat& nu,
  arma::mat& inv_omega_11,
  arma::mat& inv_c,
  arma::mat& omega_save,
  arma::mat& tau_save,
  arma::mat& mean_vec_store,
  arma::cube& inv_c_required_store,
  arma::mat const& gibbs_mat,
  arma::umat const& ind_noi_mat,
  arma::mat const& s_mat,
  const double shape,
  const double* scale_vec
) {
  
  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Calculate this constant that is used in GHS  */
  const double lambda_sq = pow(lambda, 2);

  /* Use existing global memory to avoid constant reallocaiton  */
  arma::vec flex_mem = arma::vec(g_vec1, p - 1, false, true);
  arma::vec solve_for = arma::vec(g_vec2, p - 1, false, true);

  /* Iterate p times and update col p and row p at each iteration */
  for (unsigned int i = 0; i < p; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Update preparation step is dependent on prior  */

    /* Sample from prior specific gamma density */
    const double gamma_param = g_rgamma.GetSample(shape, scale_vec[i]);

    /* Inverse of omega excluding row i and col i can be solved in O(n^2) */
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, cur_sigma, p, i);

    /* Wishart case */
    if (prior == WISHART) {

      inv_c = (s_mat.at(i, i) + 1) * inv_omega_11;
      
      /* Solve for vector is just S_21 in Wishart  */
      for (unsigned int j = 0; j < (p - 1); j++) {
        solve_for[j] = s_mat.at(i, ind_noi[j]);
      }
      
    }

    /* BGL case */
    else if (prior == BGL) {

      inv_c = inv_omega_11 * (s_mat.at(i, i) + lambda);
      for (unsigned int j = 0; j < (p - 1); j++) {

        /* Update inv_c diagonal  */
        inv_c.at(j, j) += (1 / tau.at(ind_noi[j], i));

        /* update solve for */
        solve_for[j] = s_mat.at(ind_noi[j], i) + (
          gibbs_mat.at(ind_noi[j], i) / tau.at(ind_noi[j], i)
        );
      }
    }

    /* GHS case */
    else if (prior == GHS) {
      
      inv_c = inv_omega_11 * (s_mat.at(i, i) + (1 / lambda));
      for (unsigned int j = 0; j < (p - 1); j++) {

        /* Update inv_c diagonal  */
        inv_c.at(j, j) += (1 / (tau.at(ind_noi[j], i) * lambda_sq));

        /* update solve for */
        solve_for[j] = s_mat.at(ind_noi[j], i) + (
          gibbs_mat.at(ind_noi[j], i) / (tau.at(ind_noi[j], i) * lambda_sq)
        );
      }
    }

    /* Save inv_c before solving for chol(inv_c)  */
    if (((iter - burnin) >= 0) && (i == (p - 1))) {
      inv_c_required_store.slice(iter - burnin) = inv_c;
    }

    /* -mu_i = solve(inv_c, solve_for), store chol(inv_c) in the pointer of inv_c */
    solve_for = arma::solve(inv_c, solve_for, arma::solve_opts::likely_sympd);
    inv_c = arma::chol(inv_c);

    /* Save mu_i before memory is used to calculate beta  */
    if (((iter - burnin) >= 0) && (i == (p - 1))) {
      mean_vec_store.col(iter - burnin) = -solve_for;
    }

    /* Generate random normals needed to solve for beta */
    flex_mem.randn();

    /* Solve chol(inv_c) x = randn(), store result in flex_mem  */
    flex_mem = arma::solve(arma::trimatu(inv_c), flex_mem);

    beta = -solve_for + flex_mem;

    /* Update flex_mem to store inv_omega_11 %*% beta */
    update_omega_inplace(omega, inv_omega_11, beta, ind_noi, gamma_param, i, p);

    /* Update sigma */
    update_sigma_inplace(
      cur_sigma, inv_omega_11, flex_mem.memptr(), ind_noi, gamma_param, p, i
    );
    
    /* Update sampling parameters for next iteration in prior specific method */
    
    /* BGL case */
    if (prior == BGL) {

      /* Calculate a_gig_tau and resuse variable to sample tau_12 */
      for (unsigned int j = 0; j < (p - 1); j++) {

        double gig_tau = pow(beta[j] + gibbs_mat.at(i, ind_noi[j]), 2);
        gig_tau = 1 / gigrnd(-1.0 / 2.0, gig_tau, pow(lambda, 2));
        tau.at(ind_noi[j], i) = gig_tau;
        tau.at(i, ind_noi[j]) = gig_tau;
      }

    }

    /* GHS case */ 
    else if (prior == GHS) {

      /* Sample tau_12 and nu_12 */
      for (unsigned int j = 0; j < (p - 1); j++) {

        const double cur_rate = (
          pow(beta[j], 2) / (2 * lambda_sq) + (1 / nu.at(ind_noi[j], i))
        );

        const double cur_tau = 1 / g_rgamma.GetSample(1, 1 / cur_rate);
        const double cur_nu = 1 / g_rgamma.GetSample(1, 1 / (1 + (1 / cur_tau)));

        tau.at(ind_noi[j], i) = cur_tau;
        tau.at(i, ind_noi[j]) = cur_tau;
        nu.at(ind_noi[j], i) = cur_nu;
        nu.at(i, ind_noi[j]) = cur_nu;
      }

    }
  }

  /* If iteration is past burnin period, accumulate Omega */
  if ((iter - burnin) >= 0) {
    omega_save += omega;
    if (prior != WISHART) {
      tau_save += tau;
    }
  }

}
