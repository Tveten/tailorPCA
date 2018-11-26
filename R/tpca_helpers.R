pca <- function(cov_mat, axes = 1:p, eigen_values = FALSE) {
  p <- nrow(cov_mat)
  if (eigen_values) {
    eigen_obj <- svd(cov_mat, nv = 0)
    V <- t(eigen_obj$u[, axes])
    lambda <- eigen_obj$d[axes]
    return(list('vectors' = V, 'values' = lambda))
  } else {
    eigen_obj <- svd(cov_mat, nv = 0)
    V <- t(eigen_obj$u[, axes])
    return(V)
  }
}

standardize_cov_mat <- function(cov_mat) {
  # Standardizes a covariance matrix to become a correlation matrix.
  if(!is_cor_mat(cov_mat)) {
    sd_vec <- sqrt(diag(cov_mat))
    sd_inv_mat <- diag(1/sd_vec)
    cor_mat_orig <- sd_inv_mat %*% cov_mat %*% sd_inv_mat
  } else {
    cor_mat_orig <- cov_mat
  }
}

prop_axes_max <- function(divergence_sim) {
  n_sim <- ncol(divergence_sim)
  data_dim <- nrow(divergence_sim)
  which_axes_max <- apply(divergence_sim, 2, which.max)
  prop_argmax <- table(which_axes_max) / n_sim
  prop_argmax_full <- rep(0, data_dim)
  ind <- as.integer(names(prop_argmax))
  prop_argmax_full[ind] <- prop_argmax
  prop_argmax_full
}

which_axes <- function(prop_max, keep_prop, max_axes) {
  order_axes <- order(prop_max, decreasing = TRUE)
  cum_prop <- cumsum(prop_max[order_axes])
  n_keep <- min(sum(cum_prop <= keep_prop) + 1, max_axes)
  order_axes[1:n_keep]
}