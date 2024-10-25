#include "graphical_evidence.h"


/*
 * Update the ith column and row of Omega while also accumulating
 * beta.t() %*% inv_omega_11 in global memory g_vec1
 * Used if no SIMD capabilities are detected
 */

void update_omega_inplace_no_simd(
  arma::mat& omega,
  arma::mat const& inv_omega_11,
  arma::vec const& beta,
  arma::uvec const& ind_noi,
  const double gamma_sample,
  const unsigned int ith,
  const unsigned int p
) {

  /* Update ith col and row of omega and */
  /* calculate omega_22 = gamma_sample + (beta.t() * inv_omega_11 * beta) in g_vec1 */
  double omega_22 = gamma_sample;
  for (unsigned int j = 0; j < (p - 1); j++) {

    /* Update the col and row indices excluding the diagonal  */
    omega.at(ind_noi[j], ith) = beta[j];
    omega.at(ith, ind_noi[j]) = beta[j];

    /* Store beta.t() * inv_omega_11 in g_vec1  */
    g_vec1[j] = 0;
    for (unsigned int k = 0; k < (p - 1); k++) {

      /* First beta.t() * inv_omega_11[, k] */
      g_vec1[j] += (beta[k] * inv_omega_11.at(k, j));
    }

    /* Accumulate beta_omega[j] * beta[j] */
    omega_22 += (beta[j] * g_vec1[j]);
  }
  omega.at(ith, ith) = omega_22;
}


/*
 * Update the ith column and row of Omega while also accumulating
 * beta.t() %*% inv_omega_11 in global memory g_vec1
 * Uses explicit SIMD with support for SSE, AVX, or AVX512 if available
 */

void update_omega_inplace(
  arma::mat& omega,
  arma::mat const& inv_omega_11,
  arma::vec const& beta,
  arma::uvec const& ind_noi,
  const double gamma_sample,
  const unsigned int ith,
  const unsigned int p
) {

  /* Check for supported SIMD */

#if defined(__AVX512F__) || defined(__AVX__) || defined(__SSE2__)

  unsigned int u1 = ((p - 1) / SIMD_WIDTH) * SIMD_WIDTH;

  /* First calculate beta.t() %*% inv_omega_11 in g_vec1  */
  double double_acc[SIMD_WIDTH] = { 0 };
  for (unsigned int j = 0; j < (p - 1); j++) {

    const double* cur_inv_11_col = inv_omega_11.colptr(j);

    /* Create accumulator simd storage  */
    _simd_type acc = _simd_setzero_pd();

    /* Use SIMD to get dot product of inv_omega_11[, j] with beta */
    for (unsigned int k = 0; k < u1; k += SIMD_WIDTH) {

      /* Load inv_omega_11[k, j] */
      _simd_type m1 = _simd_loadu_pd(&cur_inv_11_col[k]);

      /* Load beta[k] */
      _simd_type m2 = _simd_loadu_pd(&beta[k]);

      /* Multiply beta[k] * inv_omega_11[k, j]  */
      m2 = _simd_mul_pd(m1, m2);

      /* Accumulate results in acc  */
      acc = _simd_add_pd(acc, m2);
    }

    /* Sum each component of SIMD accumulator */
    _simd_storeu_pd(double_acc, acc);
    g_vec1[j] = 0;
    for (unsigned int iter = 0; iter < SIMD_WIDTH; iter++) {
      g_vec1[j] += double_acc[iter];
    }

    /* Add rows that were not evenly divisible by SIMD_WIDTH */
    for (unsigned int k = u1; k < (p - 1); k++) {
      g_vec1[j] += (beta[k] * inv_omega_11.at(k, j));
    }
  }

  /* Calculate omega_22 = gamma_sample + g_vec1 %*% beta and  */
  /* also update omega ith row and ith colum                  */
  double omega_22 = gamma_sample;
  for (unsigned int j = 0; j < (p - 1); j++) {

    /* Update the col and row indices excluding the diagonal  */
    omega.at(ind_noi[j], ith) = beta[j];
    omega.at(ith, ind_noi[j]) = beta[j];

    omega_22 += (beta[j] * g_vec1[j]);
  }
  omega.at(ith, ith) = omega_22;

#else

  /* Same calculation with no explicit SIMD, i.e., one iteration at a time  */
  update_omega_inplace_no_simd(
    omega, inv_omega_11, beta, ind_noi, gamma_sample, ith, p
  );

#endif
}
