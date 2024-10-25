#include "graphical_evidence.h"

/*
 * Computationally efficient calculation of inverse of omega_11:
 * For this implementation, inverse of omega_11 can be calculated as
 * sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
 * This function is called when no supported SIMD capabilities are detected
 */

void efficient_inv_omega_11_calc_no_simd(
  arma::mat& inv_omega_11,
  arma::uvec const& ind_noi,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith
) {

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      inv_omega_11.at(j, i) = sigma.at(ind_noi[j], ind_noi[i]) - (
        sigma.at(ind_noi[j], ith) * sigma.at(ind_noi[i], ith) /
        sigma.at(ith, ith)
      );
    }
  }
}


#if defined(__AVX512F__) || defined(__AVX__) || defined(__SSE2__)

/* 
 * Explicit SIMD code to update column indices of inv_omega_11 starting
 * at col_start and ending before col_stop
 */

void efficient_inv_omega_11_calc_columns_simd(
  arma::mat& inv_omega_11,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith,
  const unsigned int col_start,
  const unsigned int col_stop,
  const unsigned int u1,
  const unsigned int u2
) {

  /* Sigma and inv_omega_11 cols will be unaligned after ith col  */
  unsigned int col_offset = col_start > ith ? 1 : 0;

  /* The ith column pointer of sigma */
  const double* sigma_col_i = sigma.colptr(ith);

  /* Load sigma[ith, ith] to SIMD register  */
  _simd_type sigma_ii = _simd_set1_pd(sigma.at(ith, ith));

  /* Iterate over cols from [col_start, col_stop): */
  for (unsigned int cur_col = col_start; cur_col < col_stop; cur_col++) {

    /* The current column pointers of sigma and inv_omega_11 */
    const double* sigma_cur_col = sigma.colptr(cur_col);
    double* inv11_cur_col = inv_omega_11.colptr(cur_col - col_offset);

    /* Load sigma[cur_col, ith] to SIMD register  */
    _simd_type m1 = _simd_set1_pd(sigma.at(cur_col, ith));

    /* Iterate over rows [0, u1) using SIMD */
    for (unsigned int cur_row = 0; cur_row < u1; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11 */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row], res);
    }

    /* Finish remaining rows using standard looping from u1 to ith */
    for (unsigned int cur_row = u1; cur_row < ith; cur_row++) {
      inv_omega_11.at(cur_row, cur_col - col_offset) = (
        sigma.at(cur_row, cur_col) - (
          sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
          sigma.at(ith, ith)
        )
      );
    }

    /* Iterate over rows from [ith + 1, p): */
    /* First SIMD iteration up to u2    */
    for (unsigned int cur_row = ith + 1; cur_row < u2; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11, note     */
      /* that storing this value is now at cur_row - 1 and cur_col  */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row - 1], res);
    }

    /* Finish remaining rows using standard looping from u2 to p */
    for (unsigned int cur_row = u2; cur_row < p; cur_row++) {
      inv_omega_11.at(cur_row - 1, cur_col - col_offset) = (
        sigma.at(cur_row, cur_col) - (
          sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
          sigma.at(ith, ith)
        )
      );
    }
  }
}

#endif


/*
 * Computationally efficient calculation of inverse of omega_11:
 * For this implementation, inverse of omega_11 can be calculated as
 * sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
 * Uses explicit SIMD with support for SSE, AVX, or AVX512 if available
 */

void efficient_inv_omega_11_calc(
  arma::mat& inv_omega_11,
  arma::uvec const& ind_noi,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith
) {

  /* Check for supported SIMD */

#if defined(__AVX512F__) || defined(__AVX__) || defined(__SSE2__)

  /* Given the SIMD capabilities of the CPU, identify the largest group */
  /* of columns and rows that can computed using SIMD before and after  */
  /* the separation at index ith                                        */
  unsigned int u1 = (ith / SIMD_WIDTH) * SIMD_WIDTH;
  unsigned int u2 = ((p - ith - 1) / SIMD_WIDTH) * SIMD_WIDTH + ith + 1;

  /* Iterate over cols from [0, ith)  */
  efficient_inv_omega_11_calc_columns_simd(
    inv_omega_11, sigma, p, ith, 0, ith, u1, u2
  );

  /* Iterate over cols from [ith + 1, p)  */
  efficient_inv_omega_11_calc_columns_simd(
    inv_omega_11, sigma, p, ith, ith + 1, p, u1, u2
  );

#else

  /* Same calculation with no explicit SIMD, i.e., one iteration at a time  */
  efficient_inv_omega_11_calc_no_simd(
    inv_omega_11, ind_noi, sigma, p, ith
  );

#endif

}
