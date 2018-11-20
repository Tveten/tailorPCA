tpca <- function(cov_mat, 
                 change_distr = 'full_uniform',
                 divergence = 'hellinger',
                 cutoff = 0.99, 
                 n_sim  = 10^3,
                 print = TRUE) {
  # TODO: How to input change distribution?
  # TODO: How to handle too large covariance matrices and correlation changes
  # (due to slowness of nearPD())
  # TODO: Add handling of different change distributions.
  # TODO: Add handling of different divergence metrics.
  
  ## Make sure inputs are as expected. ----
  assertthat::assert_that(is.numeric(cov_mat),
                          msg = 'cov_mat is not numeric.')
  assertthat::assert_that(class(cov_mat) == 'matrix',
                          msg = 'cov_mat is not of class "matrix".')
  assertthat::assert_that(isSymmetric(cov_mat), 
                          msg = 'cov_mat is not a symmetric matrix.')
  assertthat::assert_that(is_positive_definite(cov_mat),
                          msg = 'cov_mat is not positive definite (some eigenvalues are < 1e-8).')
  assertthat::assert_that(is_in_interval(cutoff, interval = c(0, 1)),
                          msg = 'cutoff must be a numeric between 0 and 1.')
  assertthat::assert_that(is_whole_number(n_sim), n_sim > 0,
                          msg = 'n_sim must be an integer larger than 0.')
  
  ## Make sure we are working on a correlation matrix from now on. ----
  if(!is_cor_mat(cov_mat)) {
    sd_vec <- sqrt(diag(cov_mat))
    sd_inv_mat <- diag(1/sd_vec)
    cor_mat_orig <- sd_inv_mat %*% cov_mat %*% sd_inv_mat
  } else {
    cor_mat_orig <- cov_mat
  }
  
  data_dim <- ncol(cor_mat_orig)
  pca_obj <- pca(cor_mat_orig, eigen_values = TRUE)
  V <- pca_obj$vectors
  
  pre_mean_proj <- rep(0, data_dim)
  pre_sd_proj <- pca_obj$values
  
  change_funcs <- get_change_distr(change_distr, data_dim)
  change_type <- change_funcs$draw_types(n_sim)
  change_sparsity <- change_funcs$draw_sparsities(n_sim)
  
  # TODO: Make a FO that records results from draw_change in a list here,
  #       to avoid having to call anything from change_funcs here.
  # QUESTION: Add if-clauses for different change types to increase speed?
  hellinger_sim <- vapply(1:n_sim, function(b) {
    post_param_orig <- draw_change(cor_mat_orig, change_funcs, 
                                   change_type[b], change_sparsity[b])
    post_mean_proj <- V %*% post_param_orig$mean
    post_sd_proj <- sqrt(apply(V, 1, function(v) v %*% post_param_orig$cov_mat %*% v))
    calc_hellinger(pre_mean_proj, pre_sd_proj, post_mean_proj, post_sd_proj)
  },
  numeric(data_dim))
  # return_list <- list(axes, which_axes, summaries, hellinger_sim)
  # if (print) print(which_axes)
  # structure(return_list, class = 'tpca')
  hellinger_sim
}

