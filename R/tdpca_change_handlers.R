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
    if (ncol(cor_mat) > 200)
      stop('The current implementation for changes in correlation is too slow for dimensions > 200')
    
    # cor_dims <- attr(cor_mat, 'which_dims_cor')
    # cor_mat_sparsity <- length(cor_dims)
    # if (cor_mat_sparsity < data_dim) {
    #   affected_dims <- sample(cor_dims, min(sparsity, cor_mat_sparsity))
    # }
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

duplicate_diag_block <- function(cov_mat, block_dim) {
  data_dim <- ncol(cov_mat)
  n_blocks <- data_dim / block_dim
  assert_natural_number(n_blocks)
  block_to_duplicate <- cov_mat[1:block_dim, 1:block_dim]
  adjusted_cov_mat <- cov_mat
  for (i in 2:n_blocks) {
    ind <- ((i - 1) * block_dim + 1):((i - 1) * block_dim + block_dim)
    adjusted_cov_mat[ind, ind] <- block_to_duplicate
  }
  adjusted_cov_mat
}
