#include "graphical_evidence.h"

/*
 * Sets the random seed for globabl samplers
 */
 // [[Rcpp::export]]
void set_seed(unsigned int seed) {
  g_rgamma.SetSeed(seed);
}
