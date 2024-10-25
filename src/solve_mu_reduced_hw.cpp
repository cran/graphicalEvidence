#include "graphical_evidence.h"


/*
 * Calculate the mean mu of the multivariate density for eq9 in Hao
 * Wang decomposition MCMC sampling.
 * Native Armadillo implementation
 */

arma::vec solve_mu_reduced_hw(
  const unsigned int i,
  arma::uvec const& reduced_one_ind,
  arma::uvec const& reduced_zero_ind,
  arma::mat const& inv_c
) {

  /* Create a vec that will be the target to solve for based on g_vec2  */
  arma::vec solve_for(g_vec2, reduced_one_ind.n_elem, false);

  /* If there are any zero values in currently considered column  */
  if (reduced_zero_ind.n_elem) {

    /* Create arma::rowvec from gvec1 */
    arma::rowvec vec_acc_21(g_vec1, reduced_zero_ind.n_elem, false);

    /* Multiply in place with relevant inv_C submat */
    vec_acc_21 *= inv_c.submat(reduced_zero_ind, reduced_one_ind);

    solve_for += vec_acc_21.t();
  }

  /* Solve for mu i reduced */
  arma::vec mu_reduced = -arma::solve(
    inv_c.submat(reduced_one_ind, reduced_one_ind), solve_for, 
    arma::solve_opts::likely_sympd
  );

  return mu_reduced;
}
