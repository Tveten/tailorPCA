#' Tailors the choice of principal axes to a set of distributional changes
#'
#' Here goes the description.
#'
#' Detailed description.
#'
#' @param cov_mat A covariance matrix, i.e., a numeric matrix that is positive
#'   definite.
#' @param change_distr A string or a change distribution object. A string can be
#'   used to choose among a set of already implemented distributions:
#'   'full_uniform', 'mean_only', 'sd_only', 'cor_only'. Custom change
#'   distributions can be specified by using the \code{\link{set_uniform_cd()}}
#'   function.
#' @param divergence A string specifying which divergence metric to use.
#' @param cutoff A numeric between 0 and 1 governing how many principal axes to
#'   retain.
#' @param max_axes
#' @param n_sim
#' @param print_which
#'
#' @return
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
                 n_sim        = 10^3,
                 print_which  = TRUE) {
  
  ## Make sure inputs are as expected. -----------------------------------------
  #  cov_mat:
  assertthat::assert_that(is.numeric(cov_mat), !any(is.na(cov_mat)))
  assertthat::assert_that(class(cov_mat) == 'matrix',
                          msg = 'cov_mat is not of class "matrix".')
  assertthat::assert_that(isSymmetric(cov_mat), 
                          msg = 'cov_mat is not a symmetric matrix.')
  assertthat::assert_that(is_positive_definite(cov_mat),
                          msg = 'cov_mat is not positive definite (some eigenvalues are < 1e-8).')
  
  #  cutoff:
  assertthat::assert_that(is_in_interval(cutoff, interval = c(0, 1)),
                          msg = 'cutoff must be a numeric between 0 and 1.')
  
  #  n_sim:
  assertthat::assert_that(is_whole_number(n_sim), n_sim > 0,
                          msg = 'n_sim must be an integer larger than 0.')
  
  #  divergence:
  divergence_func <- get_divergence(divergence)
  
  #  change_distr:
  data_dim <- ncol(cov_mat)
  if (is.character(change_distr)) {
    change_funcs <- get_change_distr(change_distr, data_dim)
  } else if (class(change_distr) == 'change_distr') {
    change_funcs <- change_distr
  } else stop(paste0('change_distr must either be a character string or belong to class change_distr. See set_uniform_cd.'))
    
  change_type <- change_funcs$draw_types(n_sim)
  change_sparsity <- change_funcs$draw_sparsities(n_sim)
  
  ## Make sure we are working on a correlation matrix from now on. -------------
  if(!is_cor_mat(cov_mat)) {
    sd_vec <- sqrt(diag(cov_mat))
    sd_inv_mat <- diag(1/sd_vec)
    cor_mat_orig <- sd_inv_mat %*% cov_mat %*% sd_inv_mat
  } else {
    cor_mat_orig <- cov_mat
  }
  
  pca_obj <- pca(cor_mat_orig, eigen_values = TRUE)
  V <- pca_obj$vectors
  pre_mean_proj <- rep(0, data_dim)
  pre_sd_proj <- sqrt(pca_obj$values)
  
  divergence_sim <- vapply(1:n_sim, function(b) {
    post_param_orig <- draw_change(cor_mat_orig, change_funcs, 
                                   change_type[b], change_sparsity[b])
    post_mean_proj <- V %*% post_param_orig$mean
    post_sd_proj <- sqrt(apply(V, 1, function(v) v %*% post_param_orig$cov_mat %*% v))
    divergence_func(pre_mean_proj, pre_sd_proj, post_mean_proj, post_sd_proj)
  },
  numeric(data_dim))
  
  prop_max <- prop_axes_argmax(divergence_sim)
  most_sensitive_axes <- which_axes(prop_max, cutoff, max_axes)
  
  return_list <- list('axes'            = V[most_sensitive_axes, ], 
                      'which_axes'      = most_sensitive_axes, 
                      'divergence_sim'  = divergence_sim,
                      'change_type'     = change_type,
                      'change_sparsity' = change_sparsity)
  if (print_which) print(most_sensitive_axes)
  invisible(structure(return_list, class = 'tpca'))
}

