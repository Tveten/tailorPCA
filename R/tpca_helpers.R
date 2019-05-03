#' Principal component analysis
#' 
#' Performs efficient eigendecomposition of an input covariance matrix based on
#' which principal axes that are wanted. If all axes are wanted, \code{\link{svd}}
#' is used. \code{\link{RSpectra::eigs_sym}} is used if only the 
#' highest or lowest eigenvalues with corresponding eigenvectors are requested.
#' 
#' @param cov_mat A covariance matrix.
#' @param axes A vector indicating which principal axes are wanted.
#' 
#' @return \code{pca} returns an S3 object of class "pca". This is a list with 
#' the following components:
#' \describe{
#'   \item{\code{vectors}}{A matrix with the chosen principal axes/eigenvectors as rows.}
#'   \item{\code{values}}{A vector of the corresponding eigenvalues}
#' }
#'
#' @export
pca <- function(cov_mat, axes = 1:data_dim) {
  data_dim <- ncol(cov_mat)
  range_axes <- c(min(axes), max(axes))
  if (range_axes[1] >= data_dim/2) {
    pca_obj <- pca_small(cov_mat, axes)
  } else if (range_axes[2] <= data_dim/2) {
    pca_obj <- pca_large(cov_mat, axes)
  } else {
    pca_obj <- pca_all(cov_mat, axes)
  }
  structure(pca_obj, class = 'pca')
}

pca_all <- function(cov_mat, axes = 1:ncol(cov_mat)) {
  eigen_obj <- svd(cov_mat, nv = 0)
  eigen_vectors <- t(eigen_obj$u)[axes, , drop = FALSE]
  eigen_values <- eigen_obj$d[axes]
  list('vectors' = eigen_vectors, 'values' = eigen_values)
}

pca_small <- function(cov_mat, axes) {
  data_dim <- ncol(cov_mat)
  k <- data_dim - min(axes) + 1
  shifted_axes <- axes - min(axes) + 1
  eigen_obj <- RSpectra::eigs_sym(cov_mat, k = k, which = 'LM', sigma = 0)
  eigen_vectors <- t(eigen_obj$vectors)[shifted_axes, , drop = FALSE]
  eigen_values <- eigen_obj$values[shifted_axes]
  list('vectors' = eigen_vectors, 'values' = eigen_values)
}

pca_large <- function(cov_mat, axes) {
  k <- max(axes)
  eigen_obj <- RSpectra::eigs_sym(cov_mat, k = k, which = 'LM')
  eigen_vectors <- t(eigen_obj$vectors)[axes, , drop = FALSE]
  eigen_values <- eigen_obj$values[axes]
  list('vectors' = eigen_vectors, 'values' = eigen_values)
}

which_dims_cor <- function(cov_mat) {
  dim_is_ind <- function(x) {
    sum(x == 0) == length(x) - 1
  }

  data_dim <- ncol(cov_mat)
  cov_mat[abs(cov_mat) < sqrt(.Machine$double.eps)]
  ind_dims <- apply(cov_mat, 1, dim_is_ind)
  if (sum(ind_dims) == data_dim) {
    warning('cov_mat is a diagonal matrix, so trying to change correlations will result in error.')
    return(0)
  }
  (1:data_dim)[!ind_dims]
}

#' @export
standardize_cov_mat <- function(cov_mat) {
  # Standardizes a covariance matrix to become a correlation matrix.
  if(!is_cor_mat(cov_mat)) {
    sd_vec <- sqrt(diag(cov_mat))
    sd_inv_mat <- diag(1/sd_vec)
    cor_mat_orig <- sd_inv_mat %*% cov_mat %*% sd_inv_mat
  } else {
    cor_mat_orig <- cov_mat
  }
  cor_mat_orig
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
  n_keep <- min(sum(cum_prop < keep_prop) + 1, max_axes)
  order_axes[1:n_keep]
}

obs_needed <- function(cor_mat, alpha) {
  llr <- function(l1, l2) {
    log(l1) + log(l2) - 2 * log((l1 + l2) / 2) 
  }
  d <- ncol(cor_mat)
  l <- svd(cor_mat, nv = 0, nu = 0)$d
  l <- l[l < 1]
  n_tests <- length(l) - 1
  ind <- 1:n_tests
  log_V <- vapply(ind, function(i) llr(l[i], l[i + 1]), numeric(1))
  c <- stats::qchisq(1 - alpha / n_tests, 2)
  U <- - c * log_V + d
  max(U)
}
