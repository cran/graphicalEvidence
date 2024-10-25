#include "graphical_evidence.h"

/* Global info storage/args for manual LAPACK/BLAS calls  */
int info_int;
char uplo = 'U';
double one = 1.0;
int nrhs = 1;

/* Define global static storage */
int g_ipv[MAX_DIM];
double g_vec1[MAX_DIM];
double g_vec2[MAX_DIM];
double g_mat1[MAX_DIM * MAX_DIM];
