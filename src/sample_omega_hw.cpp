#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates omega accumulator, inv_c accumulator, and mean_vec accumulator as 
 * determined by iterations exceeding burnin period, G-Wishart prior
 */

void sample_omega_hw(
  const int iter,
  const int burnin,
  const int alpha,
  arma::vec& beta,
  arma::mat& omega,
  arma::mat& inv_omega_11,
  arma::mat& inv_c,
  arma::mat& omega_save,
  arma::mat& mean_vec_store,
  arma::cube& inv_c_required_store,
  arma::mat const& gibbs_mat,
  arma::mat const& g_mat_adj,
  arma::umat const& ind_noi_mat,
  std::vector<arma::uvec> const& find_which_ones,
  std::vector<arma::uvec> const& find_which_zeros,
  arma::mat const& scale_mat,
  arma::mat const& s_mat,
  arma::mat& sigma,
  const double shape_param,
  const double* scale_params
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  for (arma::uword i = 0; i < p; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Number of ones in current col */
    const arma::uword reduced_dim = find_which_ones[i].n_elem;

    /* Generate random gamma sample based */
    const double gamma_sample = g_rgamma.GetSample(
      shape_param, scale_params[i]
    );

    /* Inverse of omega excluding row i and col i can be solved in O(n^2) */
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, sigma, p, i);

    /* Update ith row and col of omega by calculating beta, if i == (p - 1),  */
    /* and burnin is completed, save results in accumulator variables         */

    /* Fill in any zero indices (may be empty) to beta  */
    for (unsigned int j = 0; j < find_which_zeros[i].n_elem; j++) {

      /* Current index of zero  */
      const unsigned int which_zero = find_which_zeros[i][j];

      /* Store vec_acc_21 in g_vec1 */
      const double extracted_gibbs = -gibbs_mat.at(ind_noi[which_zero], i);
      g_vec1[j] = extracted_gibbs;

      /* Set vec to relevant elements of beta */
      beta[which_zero] = extracted_gibbs;
    }

    /* Calculate mean mu and assign to one indices if they exist */
    if (reduced_dim) {

      /* Case where some ones are found in current col  */
      inv_c = inv_omega_11 * (scale_mat.at(i, i) + s_mat.at(i, i));

      /* Manual memory management, initialize solve for vector to find mu reduced */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        int row_index = ind_noi_mat.at(find_which_ones[i][j], i);
        g_vec2[j] = s_mat.at(row_index, i) + scale_mat.at(row_index, i);
      }

      /* Solve for mu_i in place  */
      arma::vec mu_reduced = solve_mu_reduced_hw(
        i, find_which_ones[i], find_which_zeros[i], inv_c
      );

      /* Store results if the column considered (i) is the last (p) */
      if (((iter - burnin) >= 0) && (i == (p - 1))) {

        inv_c_required_store.slice((iter - burnin)) = inv_c.submat(
          find_which_ones[i], find_which_ones[i]
        );

        double* cur_col = mean_vec_store.colptr(iter - burnin);
        memcpy(cur_col, mu_reduced.memptr(), reduced_dim * sizeof(double));
      }

      /* Get chol(inv_c[which_ones, which_ones])  */
      arma::mat chol_inv_c(g_mat1, reduced_dim, reduced_dim, false);
      chol_inv_c = inv_c.submat(find_which_ones[i], find_which_ones[i]);
      chol_inv_c = arma::chol(chol_inv_c);

      /* Generate random normals in g_vec1  */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        g_vec1[j] = arma::randn();
      }
      arma::vec flex_mem(g_vec1, reduced_dim, false);

      /* Solve chol(inv_c) x = randn(), store result in g_vec1  */
      flex_mem = arma::solve(arma::trimatu(chol_inv_c), flex_mem);

      /* Update one indices of beta with mu_i + solve(chol(inv_c_ones), randn())  */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        beta[find_which_ones[i][j]] = mu_reduced[j] + g_vec1[j];
      }
    }

    /* Update omega in place using newly calculated beta  */
    update_omega_inplace(omega, inv_omega_11, beta, ind_noi, gamma_sample, i, p);

    /* After omega is updated, update sigma where beta.t() %*% inv_omega_11       */
    /* is still stored in global memory g_vec1 from our previous update of omega  */
    update_sigma_inplace(
      sigma, inv_omega_11, g_vec1, ind_noi, gamma_sample, p, i
    );

  }

  /* If iteration is past burnin period, accumulate Omega */
  if ((iter - burnin) >= 0) {
    omega_save += omega;
  }
  
}
