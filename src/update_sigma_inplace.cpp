#include "graphical_evidence.h"


/*
 * Update sigma for either BGL or GHS prior in place using:
 * sigma_11 = inv_omega_11 + (omega_beta %*% t(omega_beta) / gamma_param)
 * sigma_12 = -omega_beta / gamma_param
 * sigma_22 = 1 / gamma_param
 */

void update_sigma_inplace_no_simd(
  arma::mat& sigma,
  arma::mat const& inv_omega_11,
  double* omega_beta,
  arma::uvec const& ind_noi,
  const double gamma_param,
  const unsigned int p,
  const unsigned int ith
) {

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      /* Update indices in sigma_11 */
      sigma.at(ind_noi[j], ind_noi[i]) = inv_omega_11.at(j, i) + (
        omega_beta[j] * omega_beta[i] / gamma_param
      );
    }

    /* Update indices in sigma_12 */
    sigma.at(ind_noi[i], ith) = -omega_beta[i] / gamma_param;
    sigma.at(ith, ind_noi[i]) = -omega_beta[i] / gamma_param;
  }

  /* Update sigma_22  */
  sigma.at(ith, ith) = 1 / gamma_param;
}


#if defined(__AVX512F__) || defined(__AVX__) || defined(__SSE2__)

/*
 * Explicit SIMD code to update sigma in place from column indices
 * col_start to col_stop
 */

void update_sigma_inplace_columns_simd(
  arma::mat& sigma,
  arma::mat const& inv_omega_11,
  double* omega_beta,
  const double gamma_param,
  const unsigned int p,
  const unsigned int ith,
  const unsigned int col_start,
  const unsigned int col_stop,
  const unsigned int u1,
  const unsigned int u2
) {

  /* Sigma and inv_omega_11 cols will be unaligned after ith col  */
  unsigned int col_offset = col_start > ith ? 1 : 0;

  /* Load gamma param to SIMD register  */
  _simd_type gamma_val = _simd_set1_pd(gamma_param);

  /* Iterate over cols from [col_start, col_stop): */
  for (unsigned int cur_col = col_start; cur_col < col_stop; cur_col++) {

    /* The current column pointers of sigma and inv_omega_11  */
    double* sigma_cur_col = sigma.colptr(cur_col);
    const double* inv11_cur_col = inv_omega_11.colptr(cur_col - col_offset);

    /* Load omega_beta[cur_col - col_offset] to SIMD register */
    _simd_type m1 = _simd_set1_pd(omega_beta[cur_col - col_offset]);

    /* Iterate over rows [0, u1) using SIMD */
    for (unsigned int cur_row = 0; cur_row < u1; cur_row += SIMD_WIDTH) {

      /* Load omega_beta[cur_row] */
      _simd_type m2 = _simd_loadu_pd(&omega_beta[cur_row]);

      /* Multiply omega_beta[cur_row] * omega_beta[cur_col - col_offset] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by gamma_val  */
      m2 = _simd_div_pd(m2, gamma_val);

      /* Load in inv_omega_11[cur_row, cur_col - col_offset]  */
      _simd_type res = _simd_loadu_pd(&inv11_cur_col[cur_row]);

      /* Calculate final result to update in sigma  */
      res = _simd_add_pd(res, m2);
      _simd_storeu_pd(&sigma_cur_col[cur_row], res);
    }

    /* Finish remaining rows using standard looping from u1 to ith */
    for (unsigned int cur_row = u1; cur_row < ith; cur_row++) {
      sigma.at(cur_row, cur_col) = (
        inv_omega_11.at(cur_row, cur_col - col_offset) + (
          omega_beta[cur_row] * omega_beta[cur_col - col_offset] / 
          gamma_param
        )
      );
    }

    /* Iterate over rows from [ith + 1, p): */
    /* First SIMD iteration up to u2    */
    for (unsigned int cur_row = ith + 1; cur_row < u2; cur_row += SIMD_WIDTH) {

      /* Load omega_beta[cur_row - 1] */
      _simd_type m2 = _simd_loadu_pd(&omega_beta[cur_row - 1]);

      /* Multiply omega_beta[cur_col - col_offset] * omega_beta[cur_row - 1] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by gamma_val  */
      m2 = _simd_div_pd(m2, gamma_val);

      /* Load in inv_omega_11[cur_row - 1, cur_col - col_offset]  */
      _simd_type res = _simd_loadu_pd(&inv11_cur_col[cur_row - 1]);
      
      /* Calculate final result to update in sigma  */
      res = _simd_add_pd(res, m2);
      _simd_storeu_pd(&sigma_cur_col[cur_row], res);
    }

    /* Finish remaining rows using standard looping from u1 to ith */
    for (unsigned int cur_row = u2; cur_row < p; cur_row++) {
      sigma.at(cur_row, cur_col) = (
        inv_omega_11.at(cur_row - 1, cur_col - col_offset) + (
          omega_beta[cur_row - 1] * omega_beta[cur_col - col_offset] / 
          gamma_param
        )
      );
    }

    /* Update sigma row and col ith */
    sigma.at(cur_col, ith) = -omega_beta[cur_col - col_offset] / gamma_param;
    sigma.at(ith, cur_col) = sigma.at(cur_col, ith);
  }
}

#endif


/*
 * Update sigma for either BGL or GHS prior in place using:
 * sigma_11 = inv_omega_11 + (omega_beta %*% t(omega_beta) / gamma_param)
 * sigma_12 = -omega_beta / gamma_param
 * sigma_22 = 1 / gamma_param
 * Uses explicit SIMD with support for SSE, AVX, or AVX512 if available
 */

void update_sigma_inplace(
  arma::mat& sigma,
  arma::mat const& inv_omega_11,
  double* omega_beta,
  arma::uvec const& ind_noi,
  const double gamma_param,
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

  /* Update columns [0, ith) using SIMD */
  update_sigma_inplace_columns_simd(
    sigma, inv_omega_11, omega_beta, gamma_param, p, ith,
    0, ith, u1, u2
  );

  /* Update columns [ith + 1, p) using SIMD */
  update_sigma_inplace_columns_simd(
    sigma, inv_omega_11, omega_beta, gamma_param, p, ith,
    ith + 1, p, u1, u2
  );

  /* Update sigma_ii  */
  sigma.at(ith, ith) = 1 / gamma_param;

#else

  /* Same calculation with no explicit SIMD, i.e., one iteration at a time  */
  update_sigma_inplace_no_simd(
    sigma, inv_omega_11, omega_beta, ind_noi, gamma_param, p, ith
  );

#endif

}
