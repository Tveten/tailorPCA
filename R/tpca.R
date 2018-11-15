tpca <- function(cov_mat, 
                 cutoff = 0.99, 
                 n_sim  = 10^3,
                 print = TRUE, ...) {
  
  # TODO: Write structure.
  # TODO: How to handle too large covariance matrices and correlation changes
  # (due to slowness of nearPD())
  
  ## Testing inputs.
  assertthat::assert_that(is.numeric(cov_mat),
                          msg = 'cov_mat is not numeric.')
  assertthat::assert_that(class(cov_mat) == 'matrix',
                          msg = 'cov_mat is not of class "matrix". Construct your matrix by matrix() or convert your object by using as.matrix()."')
  assertthat::assert_that(isSymmetric(cov_mat), 
                          msg = 'cov_mat is not a symmetric matrix.')
  assertthat::assert_that(is_positive_definite(cov_mat),
                          msg = 'cov_mat is not positive definite (some eigenvalues are < 1e-8).')
  assertthat::assert_that(is_in_interval(cutoff, interval = c(0, 1)),
                          msg = 'cutoff must be a numeric between 0 and 1.')
  assertthat::assert_that(is_whole_number(n_sim), n_sim > 0,
                          msg = 'n_sim must be an integer larger than 0.')
  
  ## Make sure we are working on a correlation matrix from now on.
  if(!is_cor_mat(cov_mat)) {
    sd_vec <- sqrt(diag(cov_mat))
    sd_inv_mat <- diag(1/sd_vec)
    cor_mat <- sd_inv_mat %*% cov_mat %*% sd_inv_mat
  } else {
    cor_mat <- cov_mat
  }
  
  p <- ncol(cor_mat)
  
  # return_list <- list(axes, which_axes, summaries, hellinger_sims)
  # if (print) print(which_axes)
  # structure(return_list, class = 'tpca')
  cor_mat
}
