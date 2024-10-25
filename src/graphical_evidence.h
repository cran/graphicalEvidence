/* graphical_evidence.h */

#pragma once

/* Parallel execution enabled through OpenMP if available */
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
  #include <omp.h>
#else
  // for machines with compilers void of openmp support
  #define omp_get_num_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_max_threads()  1
  #define omp_get_thread_limit() 1
  #define omp_get_num_procs()    1
  #define omp_set_nested(a)      // empty statement to remove the call
  #define omp_set_num_threads(a) // empty statement to remove the call
  #define omp_get_wtime()        0
#endif

/* RcppArmadillo used for some linear algebra structures  */
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

/* Package headers  */
#include <random>
#include <cmath>
#include "simd_intrinsics.h"
#include "global_storage.h"
#include "GammaSampler.h"
#include "prototypes.h"
#include "looping_dmvnrm_arma.h"

using namespace Rcpp;
