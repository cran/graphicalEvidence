#include "graphical_evidence.h"


/*
 * Initializes indices for a matrix which has indices for 
 * 1 to p, excluding ith index in each col i
 */

void initialize_indices(
  arma::umat& ind_noi_mat
) {

  /* For each column of adjacency matrix G, create relevant   */
  /* indices by excluding row i, then identifying which ones  */
  /* and zeros are present in the ith excluded linspace       */
  arma::uword n = ind_noi_mat.n_cols;
  arma::uvec noi_indices(n - 1);
  for (arma::uword i = 0; i < n; i++) {

    /* All indices except the ith row for iteration i */
    arma::uword k = 0;
    for (arma::uword j = 0; j < n; j++) {
      if (j != i) {
        ind_noi_mat.at(k++, i) = j;
      }
    }
  }

}


/*
 * Initializes indices for the following:
 * For each col i in the adjacency matrix, create and store a
 * vector of indices that excludes only the ith index.
 * When taking the indices descrribed above of the ith
 * col of the adjacency matrix, find the indices of both
 * all ones, and all zeros, storing the results in
 * the std::vectors passed as arguments
 */

void initialize_indices(
  arma::mat const& g_mat_adj,
  arma::umat& ind_noi_mat,
  std::vector<arma::uvec>& find_which_ones,
  std::vector<arma::uvec>& find_which_zeros
) {

  /* For each column of adjacency matrix G, create relevant   */
  /* indices by excluding row i, then identifying which ones  */
  /* and zeros are present in the ith excluded linspace       */
  arma::uword n = ind_noi_mat.n_cols;
  arma::uvec noi_indices(n - 1);
  for (arma::uword i = 0; i < n; i++) {
    
    /* All indices except the ith row for iteration i */
    arma::uword k = 0;
    for (arma::uword j = 0; j < n; j++) {
      if (j != i) {
        ind_noi_mat.at(k, i) = j;
        noi_indices[k++] = j;
      }
    }

    /* Check adjacency matrix */

    /* Identify which indices have 1's and store in vector  */
    find_which_ones[i] = arma::find(
      g_mat_adj.submat(noi_indices, arma::uvec({ i })) == 1
    );

    /* Identify which indices have 0's and store in vector  */
    find_which_zeros[i] = arma::find(
      g_mat_adj.submat(noi_indices, arma::uvec({ i })) == 0
    );
  }

}
