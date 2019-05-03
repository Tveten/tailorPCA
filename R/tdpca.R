#' Tailors the choice of principal components for change detection with lagged variables
#'
#' \code{tdpca} tailors the choice of principal components to keep when detection
#' of changepoints in the mean vector or covariance matrix is the aim.
#' It extends \code{\link{tpca}} by allowing the input to be a covariance matrix 
#' of a Hankel matrix (a data matrix with lagged variables stacked on top of 
#' eachother), and thus incorporate time dynamics.
#' Note that the dimension for the change distribution is the dimension of
#' the data without lagged variables.
#' See the documentation for \code{\link{tpca}} for more information.
#'
#' @param cov_mat A covariance matrix of lagged variables. Must be positive
#'   definite.
#' @param lag The number of lags used.
#' @param change_distr A string or a change distribution object. A string can be
#'   used to choose among a set of already implemented distributions:
#'   'full_uniform', 'mean_only', 'sd_only', 'cor_only'. Custom change
#'   distributions can be specified by using the \code{\link{set_uniform_cd}}
#'   function.
#' @param divergence A string specifying which divergence metric to use. 
#'   Available options: 'normal_hellinger', 'normal_KL' and 'normal_bhat'.
#' @param cutoff A numeric between 0 and 1 governing how many principal axes to
#'   retain.
#' @param max_axes An integer indicating the maximum number of axes that should be
#'   returned regardless of what the cutoff is.
#' @param n_sim An integer specifying the number of simulation runs.
#'
#' @return \code{tpca} returns an S3 object of class "tpca". This is a list with 
#' the following components:
#' \describe{
#'   \item{\code{axes}}{A matrix with the chosen principal axes as rows, 
#'   ordered in decreasing order of sensitivity.}
#'   \item{\code{which_axes}}{A vector indicating which principal axes that were
#'   chosen in decreasing order of sensitivity.}
#'   \item{\code{prop_axes_max}}{A vector with the proportion of simulations each axis
#'   was the most sensitive one.}
#'   \item{\code{divergence_sim}}{A matrix containing all the simulated draws 
#'   from the divergence metric along each principal axis. It is of dimension 
#'   data_dim x n_sim.}
#'   \item{\code{change_type}}{A character vector indicating the type of change 
#'   for each iteration of the simulation.}
#'   \item{\code{change_sparsity}}{A numeric vector indicating the sparsity of 
#'   the change for each iteration of the simulation.}
#' }
#'
#' @references 
#'
#' @examples 
#' 
#' @aliases tailorDPCA tailoredDPCA
#' 
#' @export

tdpca <- function(cov_mat, lag,
                  change_distr = 'full_uniform',
                  divergence   = 'normal_hellinger',
                  cutoff       = 0.99, 
                  max_axes     = ncol(cov_mat),
                  n_sim        = 10^3) {
  
  ## Make sure inputs are as expected. -----------------------------------------
  assert_cov_mat(cov_mat)
  
  # TODO: Error handling for lag.
  
  assert_class_length_noNA(cutoff, is.numeric, 1)
  assert_in_interval(cutoff, c(0, 1))
  
  assert_class_length_noNA(max_axes, is.numeric, 1)
  assert_in_interval(max_axes, c(0, Inf))
  
  assert_class_length_noNA(n_sim, is.numeric, 1)
  assert_natural_number(n_sim)
  
  
  ## MAIN ----------------------------------------------------------------------
  divergence_func <- get_divergence(divergence)
  
  data_dim <- ncol(cov_mat)
  unlagged_dim <- data_dim / (lag + 1)
  assert_natural_number(unlagged_dim)
  change_funcs <- get_change_distr(change_distr, unlagged_dim)
  
  # All changes are assumed to happen to standardized data.
  cor_mat_orig <- standardize_cov_mat(cov_mat)
  
  # This attribute is needed when changing correlations in case some dims ar ind.
  attr(cor_mat_orig, 'which_dims_cor') <- which_dims_cor(cor_mat_orig)
  
  pca_obj <- pca(cor_mat_orig)
  V <- pca_obj$vectors
  pre_mean_proj <- rep(0, data_dim)
  pre_sd_proj <- sqrt(pca_obj$values)
  
  change_type <- change_funcs$draw_types(n_sim)
  change_sparsity <- change_funcs$draw_sparsities(n_sim)
  divergence_sim <- vapply(1:n_sim, function(b) {
    post_param_orig <- draw_change_tdpca(cor_mat_orig, lag, change_funcs, 
                                         change_type[b], change_sparsity[b])
    post_mean_proj <- V %*% post_param_orig$mean
    post_sd_proj <- sqrt(diag(V %*% post_param_orig$cov_mat %*% t(V)))
    divergence_func(pre_mean_proj, pre_sd_proj, post_mean_proj, post_sd_proj)
  },
  numeric(data_dim))
  
  prop_axes_max <- prop_axes_max(divergence_sim)
  most_sensitive_axes <- which_axes(prop_axes_max, cutoff, max_axes)
  
  return_list <- list('axes'            = V[most_sensitive_axes, , drop = FALSE], 
                      'which_axes'      = most_sensitive_axes, 
                      'prop_axes_max'   = prop_axes_max,
                      'divergence_sim'  = divergence_sim,
                      'change_type'     = change_type,
                      'change_sparsity' = change_sparsity,
                      'divergence'      = divergence)
  invisible(structure(return_list, class = 'tpca'))
}
