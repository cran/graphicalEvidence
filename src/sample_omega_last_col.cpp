#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates current reduced omega using last col restricted sampler, this
 * function operates only on the G_Wishart prior
 */

void sample_omega_last_col(
  const unsigned int p_reduced,
  const double shape_param,
  const double* scale_params,
  const double omega_pp,
  arma::vec& beta,
  arma::mat& omega_reduced,
  arma::mat& inv_c,
  arma::mat& inv_omega_11,
  arma::mat& sigma,
  std::vector<arma::uvec> const& find_which_ones,
  std::vector<arma::uvec> const& find_which_zeros,
  arma::umat const& ind_noi_mat,
  arma::mat const& s_mat,
  arma::mat const& scale_mat,
  arma::mat const& gibbs_mat,
  arma::mat const& last_col_outer
) {

  /* Tilda parameterization sampling  */
  omega_reduced -= ((1 / omega_pp) * last_col_outer);

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Iterate through 1 to p_reduced for restricted sampler  */
  for (arma::uword i = 0; i < p_reduced; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Get sampled gamma value  */
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[i]);

    /* Update inv_omega_11  */
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, sigma, p_reduced, i);

    /* Initialize beta indices where zeros occur  */
    for (unsigned int j = 0; j < find_which_zeros[i].n_elem; j++) {

      /* Reduced zero index */
      const unsigned int which_zero = find_which_zeros[i][j];

      /* Update beta */
      beta[which_zero] = -(
        gibbs_mat.at(ind_noi[which_zero], i) +
        (last_col_outer.at(ind_noi[which_zero], i) / omega_pp)
      );
    }

    /* Number of ones in the adjacency matrix column  */
    const arma::uword reduced_dim = find_which_ones[i].n_elem;
    if (reduced_dim) {

      /* Calculate inv_c  */
      inv_c = inv_omega_11 * (scale_mat.at(i, i) + s_mat.at(i, i));

      /* Fill global memory with inv_c[which_ones, which_ones]  */
      int assign_index = 0;
      for (unsigned int j = 0; j < reduced_dim; j++) {
        for (unsigned int k = 0; k < reduced_dim; k++) {
          g_mat1[assign_index++] = inv_c.at(find_which_ones[i][k], find_which_ones[i][j]);
        }
      }

      if (find_which_zeros[i].n_elem) {

        /* Update g_vec2 to store V[ind_noi, i] + S[ind_noi, i] +                 */
        /* + Gibbs[reduced_zeros, i].t() * inv_c[reduced_zeros, reduced_ones] +   */
        /* + col_outer[reduced_zeros, i].t() * inv_c[reduced_zeros, reduced_ones] */
        for (unsigned int j = 0; j < find_which_ones[i].n_elem; j++) {

          /* Reduced one index  */
          const unsigned int which_one = ind_noi[find_which_ones[i][j]];

          /* initialize memory with V and S */
          g_vec2[j] = scale_mat.at(which_one, i) + s_mat.at(which_one, i);

          /* Loop through current row of inv_c[zeros, ones] */
          double dot1 = 0.0;
          double dot2 = 0.0;
          for (unsigned int k = 0; k < find_which_zeros[i].n_elem; k++) {
            
            /* Reduced zero index */
            const unsigned int which_zero = ind_noi[find_which_zeros[i][k]];

            /* Accumulate dot of inv_c_not_required and gibbs/last_col_outer  */
            dot1 += (-gibbs_mat.at(which_zero, i) * inv_c.at(which_zero, j));
            dot2 += (-last_col_outer.at(which_zero, i) * inv_c.at(which_zero, j));
          }
          dot2 /= omega_pp;
          g_vec2[j] += (dot1 + dot2);
        }
      }
      else {

        for (unsigned int j = 0; j < reduced_dim; j++) {
          g_vec2[j] = arma::randn();
        }
      }

      /* -mu_i = solve(inv_c, g_vec2), store chol(inv_c) in the pointer of inv_c */

      /* Create a vec that will be the target to solve for based on g_vec2  */
      arma::vec solve_for(g_vec2, reduced_dim, false);

      arma::mat chol_inv_c(g_mat1, reduced_dim, reduced_dim, false);
      arma::vec mu_reduced = -arma::solve(
        chol_inv_c, solve_for, arma::solve_opts::likely_sympd
      );
      chol_inv_c = arma::chol(chol_inv_c);

      /* Assign random normals to g_vec1 to solve for beta ones */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        g_vec1[j] = arma::randn();
      }
      arma::vec flex_mem(g_vec1, reduced_dim, false);

      /* Solve chol(inv_c) x = randn(), store result in g_vec1  */
      flex_mem = arma::solve(arma::trimatu(chol_inv_c), flex_mem);

      /* Update beta[which_ones] = difference of g_vec1 and g_vec2  */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        beta[find_which_ones[i][j]] = g_vec1[j] + flex_mem[j];
      }
    }

    /* Update omega in place using newly calculated beta  */
    update_omega_inplace(
      omega_reduced, inv_omega_11, beta, ind_noi, gamma_sample, i, p_reduced
    );

    /* After omega_reduced is updated, update sigma where beta.t() %*% inv_omega_11 */
    /* is still stored in global memory g_vec1 from our previous update of omega    */
    update_sigma_inplace(
      sigma, inv_omega_11, g_vec1, ind_noi, gamma_sample, p_reduced, i
    );

  }

  /* Update omega reduced with last col outer product */
  omega_reduced += ((1 / omega_pp) * last_col_outer);
}
