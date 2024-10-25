/* prototypes.h */

/* Functions in GammaSampler.cpp  */
double gamma_density(const double, const double, const double);

/* Functions in double looping_dmvnrm_arma.cpp */
double looping_dmvnrm_arma(
    arma::rowvec const&, arma::mat const&, arma::cube const&
);

void looping_mvpdf_process_iteration(
  arma::vec&, arma::rowvec const&, arma::cube const&, arma::mat const&, 
  arma::mat&, arma::rowvec&, arma::mat&, double const, arma::uword
);

void inplace_tri_mat_mult(
    arma::rowvec&, arma::mat const&
);

/* Functions in calc_eq_9.cpp */
double calc_eq_9(
  arma::uvec const&, arma::cube const&, arma::mat const&,
  arma::mat const&, arma::uword const, const unsigned int, 
  const unsigned int
);

/* Functions in solve_mu_reduced_hw.cpp */
arma::vec solve_mu_reduced_hw(
  const unsigned int, arma::uvec const&, arma::uvec const&,
  arma::mat const&
);

/* Functions in prior_sample_omega.cpp  */
void prior_sample_omega(
  const int, const int, const int, arma::vec&, arma::mat&,
  arma::mat&, arma::mat&, arma::cube&, arma::mat const&,
  arma::umat const&, std::vector<arma::uvec> const&,
  std::vector<arma::uvec> const&, arma::mat const&, arma::mat&
);

/* Functions in prior_sample_omega_rmatrix.cpp  */
void prior_sample_omega_rmatirx(
  const int, const int, const double, const int,
  const double, arma::vec&, arma::mat&, arma::mat&,
  arma::mat&, arma::cube&, arma::umat const&, arma::mat&,
  arma::mat&, arma::mat&
);

/* Functions in sample_omega_hw.cpp */
void sample_omega_hw(
  const int, const int, const int, arma::vec&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::mat&,
  arma::cube&, arma::mat const&, arma::mat const&, arma::umat const&,
  std::vector<arma::uvec> const&, std::vector<arma::uvec> const&,
  arma::mat const&, arma::mat const&, arma::mat&, const double,
  const double*
);

/* Functions in sample_omega_hw_rmatrix */
void sample_omega_hw_rmatrix(
  const int, const int, const int, const int, const double,
  arma::vec&, arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::mat&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::cube&, arma::mat const&,
  arma::umat const&, arma::mat const&, const double, const double*
);

/* Functions in initialize_indices.cpp  */
void initialize_indices(
  arma::umat&
);

void initialize_indices(
  arma::mat const&, arma::umat&, std::vector<arma::uvec>&,
  std::vector<arma::uvec>&
);

/* Functions in sample_omega_last_col.cpp */
void sample_omega_last_col(
  const unsigned int, const double, const double*, const double, arma::vec&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, std::vector<arma::uvec> const&,
  std::vector<arma::uvec> const&, arma::umat const&, arma::mat const&,
  arma::mat const&, arma::mat const&, arma::mat const&
);

/* Functions in sample_omega_last_col_rmatrix.cpp */
void sample_omega_last_col_rmatrix(
  const unsigned int, const double, const double*, const double, const double,
  const double, const int, arma::vec&, arma::mat&, arma::mat&, arma::mat&,
  arma::mat&, arma::mat&, arma::mat&, arma::umat const&, arma::mat const&,
  arma::mat const&, arma::mat const&
);

/* Functions in mcmc_last_col_rmatrix */
/* mcmc_last_col_rmatrix() not called internally  */
void last_col_calc_inv_omega_11_full(
  arma::mat&, arma::mat const&
);

double calc_gamma_subtractor(
  arma::vec const&, arma::mat const&
);

void last_col_prepare_sigma_reduced(
  arma::mat&, arma::mat const&, arma::vec const&, const double
);

void update_sigma_last_col(
  arma::mat&, arma::vec const&, const double
);

/* Functions in efficient_inv_omega_11_calc.cpp */
void efficient_inv_omega_11_calc(
  arma::mat&, arma::uvec const&, arma::mat const&,
  const unsigned int, const unsigned int
);

void efficient_inv_omega_11_calc_no_simd(
  arma::mat&, arma::uvec const&, arma::mat const&,
  const unsigned int, const unsigned int
);

/* Functions in update_sigma_inplace.cpp  */
void update_sigma_inplace(
  arma::mat&, arma::mat const&, double*, arma::uvec const&,
  const double, const unsigned int, const unsigned int
);

void update_sigma_inplace_no_simd(
  arma::mat&, arma::mat const&, double*, arma::uvec const&,
  const double, const unsigned int, const unsigned int
);

/* Functions in update_omega_inplace.cpp  */
void update_omega_inplace(
  arma::mat&, arma::mat const&, arma::vec const&, arma::uvec const&,
  const double, const unsigned int, const unsigned int
);

void update_omega_inplace_no_simd(
  arma::mat&, arma::mat const&, arma::vec const&, arma::uvec const&,
  const double, const unsigned int, const unsigned int
);

/* Functions in calc_eq_11.cpp  */
double calc_eq_11(
  const double, const double, const double, const unsigned int, 
  arma::vec const&
);

/* Functions in gigrnd.cpp  */
double gigrnd(
  double, double, double
);

double psi(
  double, double, double
);

double dpsi(
  double, double, double
);

double fun_g(
  double, double, double, double, double
);

/* Functions in mcmc_hw_rmatrix */
/* mcmc_hw_rmatrix() not called internally  */
void get_gamma_params_hw_rmatrix(
  double*, double*, const int, const int, const int,
  const double, arma::mat const&
);
