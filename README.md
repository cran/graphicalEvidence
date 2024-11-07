# Graphical Evidence
RcppArmadillo package implementation for the proposed technique [1] to compute marginal likelihood in Gaussian graphical models under four different priors on the precision matrix:

1. Bayesian graphical lasso (BGL)

2. Graphical horseshoe (GHS)

3. Wishart

4. G-Wishart

## Installation
This package can be installed directly from github using the R package `devtools`.

```
> install.packages("devtools")
> library(devtools)
> install_github("dp-rho/graphicalEvidence")
```

Alternatively, a version of this package can be installed from CRAN using:

```
> install.packages('graphicalEvidence')
```

It should be noted that the version distributed on CRAN is more portable, but also slower. The version distributed directly from GitHub attempts to force compilation with AVX/SSE flags if possible using the `RcppXsimd` package, and also attempts to compile with both LAPACKE and OpenBLAS. These interfaces are not available on some systems, and so the version distributed on CRAN uses RcppArmadillo exclusively for all linear algebra routines. The version distributed on CRAN may be anywhere from 2-10x slower.

## Description:
This package allows estimation of marginal likelihood in Gaussian graphical models through a novel telescoping block decomposition of the precision matrix which allows estimation of model evidence via an application of Chib's [2] method. This package also provides an MCMC prior sampler for the priors of BGL, GHS, and G-Wishart.

## Examples

### Prebuilt tests
First, we can use entirely pre-determined test parameters and a fixed random seed to verify the package is working as expected. The function `test_evidence` can be called by specifying the number of complete runs for the marginal estimation and what prior should be tested (currently implemented options are: 'Wishart', 'BGL', 'GHS', 'G_Wishart').  The dimension of all data tested using this function is set to `p=5`. 
```
> test_evidence(num_runs=10, 'Wishart')
```
The console will display the parameters tested, as well as display the results using a histogram if more than 1 run is requested. The average computation time per run will be displayed, and the entire results vector, as well as the mean and standard deviation of this vector will be returned. The version of graphicalEvidence distributed on CRAN will attach these outputs as elements of the return list rather than print them explicitly.
```
Params are: 
$x_mat
          V1        V2        V3        V4       V5
1  -0.327510 -0.180860 -0.054457 -0.342300  0.99165
2  -0.092703 -0.010816 -0.094554 -0.321180 -0.53221
3   0.752180  0.885270  0.369450  0.685170 -0.56842
4   0.179810  0.018306 -0.022333 -0.019137 -0.39041
5   0.570070 -0.147280  0.391650  0.233730 -0.25239
6  -0.195330  0.087235 -0.292460  0.346490 -0.10594
7   1.036400  0.797540 -0.473170 -0.371080 -0.71219
8  -1.219000 -1.267100 -0.114760  0.342460  0.32281
9   0.126190  0.732330 -0.402690 -0.166400  0.55104
10 -0.184100 -0.288940  0.056100  0.345210 -0.29425

$scale_mat
        V1       V2       V3       V4       V5
1 0.142860 0.035714 0.000000 0.000000 0.000000
2 0.035714 0.142860 0.035714 0.000000 0.000000
3 0.000000 0.035714 0.142860 0.035714 0.000000
4 0.000000 0.000000 0.035714 0.142860 0.035714
5 0.000000 0.000000 0.000000 0.035714 0.142860

0 runs excluded for NaN
Execution time in R program per run (seconds):
   user  system elapsed 
  0.151   0.001   0.155 
$mean
[1] -84.13947

$sd
[1] 0.04027397

$results
 [1] -84.23698 -84.13235 -84.13171 -84.14765 -84.14498 -84.13035 -84.09755 -84.08554 -84.15069 -84.13690
```

Note, in the version distributed on CRAN, the results returned from `test_evidence` will be reproducible only if a random seed is passed explicitly as an argument. In the version distributed on GitHub, the random seed is reset by each call to `test_evidence` in all cases. The pre-existing parameters used for these tests can be generated using the function `gen_params_evidence`. 
```
> gen_params_evidence('Wishart')
$x_mat
          V1        V2        V3        V4       V5
1  -0.327510 -0.180860 -0.054457 -0.342300  0.99165
2  -0.092703 -0.010816 -0.094554 -0.321180 -0.53221
3   0.752180  0.885270  0.369450  0.685170 -0.56842
4   0.179810  0.018306 -0.022333 -0.019137 -0.39041
5   0.570070 -0.147280  0.391650  0.233730 -0.25239
6  -0.195330  0.087235 -0.292460  0.346490 -0.10594
7   1.036400  0.797540 -0.473170 -0.371080 -0.71219
8  -1.219000 -1.267100 -0.114760  0.342460  0.32281
9   0.126190  0.732330 -0.402690 -0.166400  0.55104
10 -0.184100 -0.288940  0.056100  0.345210 -0.29425

$scale_mat
        V1       V2       V3       V4       V5
1 0.142860 0.035714 0.000000 0.000000 0.000000
2 0.035714 0.142860 0.035714 0.000000 0.000000
3 0.000000 0.035714 0.142860 0.035714 0.000000
4 0.000000 0.000000 0.035714 0.142860 0.035714
5 0.000000 0.000000 0.000000 0.035714 0.142860
```

### Basic use case
We now consider the top level function used to estimate marginal likelihood. The function `evidence` requires an input data matrix, the number of iterations for both burnin and sampling, the specified prior, and any prior specific parameters.  The following code performs the same estimation task as the previously shown test, but note that random variation will lead to small differences unless the random seed is held constant. 
```
# Compute the marginal likelihood of x_mat for Wishart prior using 1,000 
# burnin and 5,000 sampled values at each call to the MCMC sampler with
# 10 total runs of the estimator
# Change prior name to 'Wishart', 'BGL', 'GHS', 'G_Wishart' as desired.
> g_params <- gen_params_evidence('Wishart')
> evidence(
	xx=g_params$x_mat, burnin=1e3, nmc=5e3, prior_name='Wishart', 
	runs=10, alpha=7, V=g_params$scale_mat,
)
```
As expected, we see small variation in the mean and standard deviation of the results.
```
0 runs excluded for NaN
$mean
[1] -84.15582

$sd
[1] 0.05245773

$results
 [1] -84.23467 -84.10023 -84.11618 -84.15848 -84.15428 -84.12452 -84.22029 -84.10393 -84.12049 -84.22512
```
To make results reproducible, it is not sufficient to use the R function `set.seed`, as the compiled library linked to graphicalEvidence calls a sampler independent of the R session. The function `set_seed_evidence` must be used.  Note that by setting the evidence seed to `123456789`, we can now perfectly replicate the results returned by `test_evidence`.
```
> g_params <- gen_params_evidence('Wishart')
> set_seed_evidence(123456789)
> evidence(
	xx=g_params$x_mat, burnin=1e3, nmc=5e3, prior_name='Wishart', 
	runs=10, alpha=7, V=g_params$scale_mat,
)
$mean
[1] -84.13947

$sd
[1] 0.04027397

$results
 [1] -84.23698 -84.13235 -84.13171 -84.14765 -84.14498 -84.13035 -84.09755 -84.08554 -84.15069 -84.13690
```
### Prior sampling
The function `prior_sampling` allows a user to specify one of BGL, GHS, or G_Wishart and any related parameters and sample `burnin + nmc` iterations of an MCMC sampler on the precision matrix modified slightly from a highly similar approach used in `evidence`. This code will execute 2,000 total iterations and discard the first half for the prior of GHS with `lambda=2` and dimension `p=5`.
```
> samples <- prior_sampling(p=5, burnin=1e3, nmc=1e3, prior_name='GHS', lambda=2)
> dim(samples)
[1]    5    5 1000
> samples[, , 1000]
          [,1]      [,2]      [,3]       [,4]       [,5]
[1,] 9.7316617 1.2097120 5.4897852  1.0853305  0.4333838
[2,] 1.2097120 0.4274469 0.2543827  0.9397151  0.2727555
[3,] 5.4897852 0.2543827 5.9593794  0.1234130  0.1205787
[4,] 1.0853305 0.9397151 0.1234130  4.2562709 -1.8937773
[5,] 0.4333838 0.2727555 0.1205787 -1.8937773  9.9557340
```

## References:
[1] Bhadra, A., Sagar, K., Rowe, D., Banerjee, S., & Datta, J. (2022). Evidence Estimation in Gaussian Graphical Models Using a Telescoping Block Decomposition of the Precision Matrix. arXiv preprint [arXiv:2205.01016](https://arxiv.org/abs/2205.01016).

[2] Chib, S. (1995). Marginal likelihood from the Gibbs output. Journal of the American Statistical
Association 90, 1313–1321.
