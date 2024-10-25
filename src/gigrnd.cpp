#include "graphical_evidence.h"


/*
 * Random sampling of the generalized inverse Gaussian distribution,
 * implementation from Devroye (2014) algorithm.
 */

double gigrnd(
  double lambda, 
  double a, 
  double b
) {

  /* Setup: sample from 2 parameter version of GIG(alpha, omega)  */
  double omega = sqrt(a * b);

  bool swap = false;
  if (lambda < 0) {
    lambda = -lambda;
    swap = true;
  }
  double alpha = sqrt(pow(omega, 2) + pow(lambda, 2)) - lambda;

  /* Find s and t */
  double x = -psi(1, alpha, lambda);
  double t;
  double s;
  if ((x >= 0.5) && (x <= 2)) {
    t = 1;
    s = 1;
  }
  else if (x > 2) {
    t = sqrt(2 / (alpha + lambda));
    s = sqrt(4 / (alpha * cosh(1) + lambda));
  }
  else {
    s = log(4 / (alpha + 2 * lambda));
    t = log(4 / (alpha + 2 * lambda));
  }

  /* Generation */
  double eta = -psi(t, alpha, lambda);
  double zeta = -dpsi(t, alpha, lambda);
  double theta = -psi(-s, alpha, lambda);
  double xi = dpsi(-s, alpha, lambda);
  double p = 1 / xi;
  double r = 1 / zeta;
  double td = t - (r * eta);
  double sd = s - (p * theta);
  double q = td + sd;

  /* Loop until sample is valid */
  double sample;
  bool done = false;
  while (!done) {

    double U = arma::randu();
    double V = arma::randu();
    double W = arma::randu();

    if (U < (q / (p + q + r))) {
      sample = -sd + (q * V);
    }
    else if (U < ((q + r) / (p + q + r))) {
      sample = td - (r * log(V));
    }
    else {
      sample = -sd + (p * log(V));
    }

    double f1 = exp(-eta - zeta * (sample - t));
    double f2 = exp(-theta + xi * (sample + s));
    if ((W * fun_g(sample, sd, td, f1, f2)) <=
        (exp(psi(sample, alpha, lambda)))) {
      done = true;
    }
  }
  sample = exp(sample) * (lambda / omega + sqrt(1 + pow(lambda / omega, 2)));
  if (swap) {
    sample = 1 / sample;
  }
  return(sample / sqrt(a / b));
}


/* Local functions for sampler  */
double psi(double x, double alpha, double lambda) {
  return(-alpha * (cosh(x) - 1) - lambda * (exp(x) - x - 1));
}

double dpsi(double x, double alpha, double lambda) {
  return(-alpha * sinh(x) - lambda * (exp(x) - 1));
}

double fun_g(double x, double sd, double td, double f1, double f2) {
  double rval = 0;
  if ((x >= -sd) && (x <= td)) {
    rval = 1;
  }
  else if (x > td) {
    rval = f1;
  }
  else if (x < -sd) {
    rval = f2;
  }
  return(rval);
}
