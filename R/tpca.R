#' Tailors the choice of principal components for change detection
#'
#' \code{tpca} tailors the choice of principal components to keep when detection
#' of changepoints in the mean vector or covariance matrix is the aim. The
#' choice of principal axes to project data onto is based on a normal state
#' covariance matrix and a distribution over relevant change scenarios.
#'
#' This method is based on simulating changes to a distribution, followed by 
#' measuring the principal axes' sensitivity to each change by a statistical 
#' divergence. The most sensitive axis is recorded in each simulated change to 
#' estimate the probability of an axis being the most sensitive one over the
#' range of changes specified by a change distribution. 
#' 
#' Custom change distributions can be built by using the function 
#' \code{\link{set_uniform_cd}}. All components of the distribution are uniform,
#' but the probability/importance of each type of change can be specified, along
#' with the sparsity of the change, and all the sizes and directions of the 
#' changes. In each simulation run, after a change sparsity has been drawn,
#' which dimensions that are affected by a change is always randomized.
#' The more information about which changes that are of interest, the better and
#' less general the choice of axes will be.
#' 
#' See the references for more information.
#'
#' @param cov_mat A covariance matrix, i.e., a numeric matrix that is positive
#'   definite.
#' @param change_distr A string or a change distribution object. A string can be
#'   used to choose among a set of already implemented distributions:
#'   'full_uniform', 'mean_only', 'sd_only', 'cor_only'. Custom change
#'   distributions can be specified by using the \code{\link{set_uniform_cd}}
#'   function.
#' @param divergence A string specifying which divergence metric to use. 
#'   Available options: 'normal_hellinger'.
#' @param cutoff A numeric between 0 and 1 governing how many principal axes to
#'   retain.
#' @param max_axes An integer indicating the maximum number of axes that should be
#'   returned regardless of what the cutoff is.
#' @param n_sim An integer specifying the number of simulation runs.
#' @param print_which
#'
#' @return \code{tpca} returns an S3 object of class "tpca". This is a list with 
#' the following components:
#' \describe{
#'   \item{\code{axes}}{A matrix with the chosen principal axes as rows, 
#'   ordered in decreasing order of sensitivity.}
#'   \item{\code{which_axes}}{A vector indicating which principal axes that were
#'   chosen.}
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
#' @aliases tailorPCA tailoredPCA
#' 
#' @export

tpca <- function(cov_mat, 
                 change_distr = 'full_uniform',
                 divergence   = 'normal_hellinger',
                 cutoff       = 0.99, 
                 max_axes     = ncol(cov_mat),
                 n_sim        = 10^3) {
  
  ## Make sure inputs are as expected. -----------------------------------------
  assert_cov_mat(cov_mat)
  
  assert_class_length_noNA(cutoff, is.numeric, 1)
  assert_in_interval(cutoff, c(0, 1))
  
  assert_class_length_noNA(max_axes, is.numeric, 1)
  assert_in_interval(max_axes, c(0, Inf))
  
  assert_class_length_noNA(n_sim, is.numeric, 1)
  assert_natural_number(n_sim)
  
  assert_class_length_noNA(print_which, is.logical, 1)
  
  divergence_func <- get_divergence(divergence)
  
  data_dim <- ncol(cov_mat)
  change_funcs <- get_change_distr(change_distr, data_dim)
  
  ## MAIN ----------------------------------------------------------------------
  cor_mat_orig <- standardize_cov_mat(cov_mat)
  pca_obj <- pca(cor_mat_orig, eigen_values = TRUE)
  V <- pca_obj$vectors
  pre_mean_proj <- rep(0, data_dim)
  pre_sd_proj <- sqrt(pca_obj$values)
  
  change_type <- change_funcs$draw_types(n_sim)
  change_sparsity <- change_funcs$draw_sparsities(n_sim)
  divergence_sim <- vapply(1:n_sim, function(b) {
    post_param_orig <- draw_change(cor_mat_orig, change_funcs, 
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
                      'divergence_sim'  = divergence_sim,
                      'change_type'     = change_type,
                      'change_sparsity' = change_sparsity)
  invisible(structure(return_list, class = 'tpca'))
}

