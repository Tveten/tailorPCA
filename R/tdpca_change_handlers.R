draw_change_tdpca <- function(cor_mat, lag, change_funcs, change_type, change_sparsity) {
  data_dim <- ncol(cor_mat)
  unlagged_dim <- data_dim / (lag + 1)
  affected_dims <- change_funcs$draw_dims(change_sparsity)
  if (change_type == 'mean') {
    post_mean_orig_unlag <- rep(0, unlagged_dim)
    post_mean_orig_unlag[affected_dims] <- change_funcs$draw_mean(change_sparsity)
    post_mean_orig <- rep(post_mean_orig_unlag, lag + 1)
    post_cov_mat_orig <- cor_mat
  } else if (change_type == 'sd') {
    post_mean_orig <- rep(0, data_dim)
    post_cov_mat_orig <- change_cor_mat_tdpca(cor_mat, lag,
                                              affected_dims, 
                                              draw_sd = change_funcs$draw_sd)
  } else if (change_type == 'cor') {
    post_mean_orig <- rep(0, data_dim)
    post_cov_mat_orig <- change_cor_mat_tdpca(cor_mat, lag,
                                              affected_dims, 
                                              draw_cor = change_funcs$draw_cor)
  }
  list('mean' = post_mean_orig,
       'cov_mat' = post_cov_mat_orig)
}

#' Draw a changed correlation matrix for use in dynamic TPCA.
#' 
#' Changes an input correlation matrix that is based on lagged samples. 
#' This function assumes that the correlation matrix is based on
#' observations \eqn{X_t = vec(x_{t - l}, x_{t - l + 1}, \ldots, x_t)}.
#' 
#' The lagged observations \eqn{X_t} imposes a block structure on its
#' correlation matrix, where the l + 1 diagonal blocks are all estimates of the
#' correlation matrix of \eqn{x_t}. \code{change_cor_mat_tdpca} takes a set of 
#' indices indicating affected dimensions of \eqn{x_t}, and translates this into
#' corresponding changes in the input correlation matrix of \eqn{X_t}.
#' A change in one dimension of \eqn{x_t} (the unlagged, original data)
#' corresponds to changes in (lag + 1) dimensions of \eqn{X_t}.
#' See \code{\link{change_cor_mat}} for additional information.
#' 
#' @param cor_mat A correlation matrix of lagged observations to be changed.
#' @param lag The lag used to obtain the lagged observations.
#' @param affected_dims A vector specifying which unlagged dimensions 
#' (dimensions of \eqn{x_t}) that should be changed.
#' @param do_nearPD A logical indicating whether the Matrix::nearPD function should
#' be run on the changed correlation matrix to find the closest positve
#' definite matrix to it. Highly recommended, as the changes in
#' correlation are not guaranteed to result in a valid correlation matrix.
#' @param draw_cor A function to draw n (any natural number) changes in correlation from.
#' @param draw_sd A function to draw n (any natural number) changes in standard deviation from.
#' 
#' @return A changed correlation matrix, guaranteed to be positive definite if
#' do_nearPD = TRUE.
#' 
#' @export
change_cor_mat_tdpca <- function(cor_mat, lag, affected_dims, do_nearPD = TRUE,
                                 draw_cor = NULL, draw_sd = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_cor or draw_sigma
  #
  # Returns:
  #   Sigma2: The change covariance matrix.
  msg <- 'ERROR: Either a variance or a correlation change distribution must be specified'
  assertthat::assert_that(!is.null(draw_cor) || !is.null(draw_sd) , msg = msg)
  
  data_dim <- ncol(cor_mat)
  unlagged_dim <- data_dim / (lag + 1)
  change_sparsity <- length(affected_dims)
  
  n_blocks <- data_dim / unlagged_dim
  assert_natural_number(n_blocks)
  post_cov_mat <- cor_mat
  for (i in 1:n_blocks) {
    ind <- ((i - 1) * block_dim + 1):((i - 1) * block_dim + block_dim)
    post_cov_mat[ind, ind] <- change_cor_mat(cor_mat[ind, ind], affected_dims,
                                             draw_cor = draw_cor, draw_sd = draw_sd)
  }
  post_cov_mat
}

change_cor_mat_tdpca <- function(cor_mat, lag, affected_dims, do_nearPD = TRUE,
                                 draw_cor = NULL, draw_sd = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_cor or draw_sigma
  #
  # Returns:
  #   Sigma2: The changed covariance matrix.
  change_cor <- function(cor_mat, lag, draw_cor, sparsity) {
    if (length(affected_dims) < 2)
      stop('For changes in correlation, the number of affected dimensions must be >= 2')
    
    ind <- t(utils::combn(affected_dims, 2))
    change_factor <- draw_cor(nrow(ind))
    
    post_cor_mat <- cor_mat
    for (i in 1:(lag + 1)) {
      moved_ind <- (i - 1) * unlagged_dim + ind
      post_cor <- change_factor * cor_mat[moved_ind]
      post_cor_mat[moved_ind] <- post_cor
      post_cor_mat[moved_ind[, c(2, 1), drop = FALSE]] <- post_cor
    }
    
    if (do_nearPD) {
      maxit <- max(200 - data_dim, 0)
      post_cor_mat <- as.matrix(Matrix::nearPD(post_cor_mat,
                                               corr     = TRUE,
                                               maxit    = maxit,
                                               do2eigen = TRUE,
                                               posd.tol = 1e-8)$mat)
    }
    post_cor_mat
  }
  
  msg <- 'ERROR: Either a variance or a correlation change distribution must be specified'
  assertthat::assert_that(!is.null(draw_cor) || !is.null(draw_sd) , msg = msg)
  
  data_dim <- ncol(cor_mat)
  unlagged_dim <- data_dim / (lag + 1)
  change_sparsity <- length(affected_dims)
  
  post_cov_mat <- cor_mat
  
  # Correlation change handling
  if (!is.null(draw_cor)) {
    post_cov_mat <- change_cor(cor_mat, lag, draw_cor, change_sparsity)
  }
  
  # Variance change handling
  if (!is.null(draw_sd)) {
    post_sd_unlag <- rep(1, unlagged_dim)
    post_sd_unlag[affected_dims] <- draw_sd(change_sparsity)
    post_sd <- rep(post_sd_unlag, lag + 1)
    diag_post_sd <- diag(post_sd)
    post_cov_mat <- diag_post_sd %*% post_cov_mat %*% diag_post_sd
  }
  
  post_cov_mat
}
