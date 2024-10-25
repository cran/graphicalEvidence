#pragma once

/* global_storage.h */

/* Global constants for prior coding and size */
#define WISHART (0)
#define BGL (1)
#define GHS (2)
#define MAX_DIM (200)

/* Global static storage */
extern int g_ipv[MAX_DIM];
extern double g_vec1[MAX_DIM];
extern double g_vec2[MAX_DIM];
extern double g_mat1[MAX_DIM * MAX_DIM];
