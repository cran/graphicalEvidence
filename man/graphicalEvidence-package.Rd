\name{graphicalEvidence-package}
\alias{graphicalEvidence-package}
\alias{graphicalEvidence}
\docType{package}
\title{
  Compute Marginal Likelihood in Gaussian Graphical Models
}
\description{
  This package allows estimation of marginal likelihood in Gaussian graphical 
  models through a novel telescoping block decomposition of the precision 
  matrix which allows estimation of model evidence via an application of 
  Chib's method. The currently implemented priors are: Bayesian graphical lasso 
  (BGL), Graphical horseshoe (GHS), Wishart, and G-Wishart. The top level 
  function used to estimate marginal likelihood is \code{\link{evidence}} which 
  expects the prior name, data, and relevant prior specific parameters. This 
  package also provides an MCMC prior sampler for the priors of BGL, GHS, and 
  G-Wishart, implemented in \code{\link{prior_sampling}}, which expects a prior 
  name and prior specific parameters. Both functions also expect the number of 
  burn-in iterations and the number of sampling iterations for the underlying 
  MCMC sampler.
  
  Bhadra, A., Sagar, K., Rowe, D., Banerjee, S., & Datta, J. (2022) 
  "Evidence Estimation in Gaussian Graphical Models Using a Telescoping Block 
  Decomposition of the Precision Matrix." <https://arxiv.org/abs/2205.01016>

  Chib, S. "Marginal likelihood from the Gibbs output." (1995) 
  <https://www.jstor.org/stable/2291521>
}
\details{
  This package implements marginal estimation for four priors, "Wishart"", 
  Bayesian Graphical Lasso ("BGL"), graphical horseshoe ("GHS"), and 
  "G-Wishart". An MCMC prior sampler is also provided for "BGL", "GHS", and
  "G-Wishart".
  
  For more information and a faster, less portable implementation, visit the 
  package repository on GitHub: 
  \url{https://github.com/dp-rho/graphicalEvidence}
}
\author{
Maintainer: David Rowe <david@rowe-stats.com>
}
\references{
  Bhadra, A., Sagar, K., Rowe, D., Banerjee, S., & Datta, J. (2022) 
  "Evidence Estimation in Gaussian Graphical Models Using a Telescoping Block 
  Decomposition of the Precision Matrix." <https://arxiv.org/abs/2205.01016>
  
  Chib, S. "Marginal likelihood from the Gibbs output." (1995) 
  <https://www.jstor.org/stable/2291521>
}
\keyword{ package }
\seealso{
  \code{\link{test_evidence}}: For basic example of functionality
  
  \code{\link{evidence}}: For top level estimation function
  
  \code{\link{prior_sampling}}: For the prior sampler function
}
\examples{
  test_results <- test_evidence(num_runs=3, prior_name='G_Wishart') 
}
