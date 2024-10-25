#include "graphical_evidence.h"

/*
 * Sets the global number of cores available for parallel programming
 * using OpenMP, this value is used for all parallel computation 
 */
// [[Rcpp::export]]
void set_cores(const int cores) {
  omp_set_num_threads(cores);
}
