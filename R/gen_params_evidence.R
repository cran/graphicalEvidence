#' @title Generate Test Parameters
#' 
#' @description
#' Generates predetermined parameters for testing the functionality of
#' the graphical evidence method
#' 
#' @param prior_name The name of the prior for being tested with preexisting
#' test parameters, this is one of 'Wishart', 'BGL', 'GHS', 'G_Wishart'
#' 
#' @returns A list of matrices representing test parameters dependent on the 
#' prior specified in prior_name
#' 
#' @examples
#' # Generate test parameter matrices for G-Wishart prior 
#' gen_params_evidence('G_Wishart')
#' @export
gen_params_evidence <- function(
    prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart')    
) {

  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Get package directory
  package_dir <- system.file(package = "graphicalEvidence")
  
  # Read in prior specific test parameters
  params <- switch(
    prior_name,
    
    'Wishart' = list(
      x_mat = read.csv(file.path(
        package_dir, 'test_params', 'X_mat_Wishart_q_5_n_10_alpha_7.csv'
      ), header=FALSE),
      scale_mat = read.csv(file.path(
        package_dir, 'test_params', 'Scale_mat_Wishart_q_5_n_10_alpha_7.csv'
      ), header=FALSE)
    ),
    
    'BGL' = list(
      x_mat = read.csv(file.path(
        package_dir, 'test_params', 'X_mat_BGL_q_5_n_10_lambda_1.csv'
      ), header=FALSE)
    ),
    
    'GHS' = list(
      x_mat = read.csv(file.path(
        package_dir, 'test_params', 'X_mat_GHS_q_5_n_10_lambda_1.csv'
      ), header=FALSE)
    ),
    
    'G_Wishart' = list(
      g_mat = read.csv(file.path(
        package_dir, 'test_params', 'G_mat_q_5_n_10_delta_2.csv'
      ), header=FALSE),
      scale_mat = read.csv(file.path(
        package_dir, 'test_params', 'Scale_mat_q_5_n_10_delta_2.csv'
      ), header=FALSE),
      x_mat = read.csv(file.path(
        package_dir, 'test_params', 'X_mat_q_5_n_10_delta_2.csv'
      ), header=FALSE)
    )
  )
  
  return(params)
}