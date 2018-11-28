draw_change <- function(cor_mat, change_funcs, change_type, change_sparsity) {
  data_dim <- ncol(cor_mat)
  affected_dims <- change_funcs$draw_dims(change_sparsity)
  if (change_type == 'mean') {
    post_mean_orig <- rep(0, data_dim)
    post_mean_orig[affected_dims] <- change_funcs$draw_mean(change_sparsity)
    post_cov_mat_orig <- cor_mat
  } else if (change_type == 'sd') {
    post_mean_orig <- rep(0, data_dim)
    post_cov_mat_orig <- change_cor_mat(cor_mat, 
                                        affected_dims, 
                                        draw_sd = change_funcs$draw_sd)
  } else if (change_type == 'cor') {
    post_mean_orig <- rep(0, data_dim)
    post_cov_mat_orig <- change_cor_mat(cor_mat, 
                                        affected_dims, 
                                        draw_cor = change_funcs$draw_cor)
  }
  list('mean' = post_mean_orig,
       'cov_mat' = post_cov_mat_orig)
}

change_cor_mat <- function(cor_mat, affected_dims, 
                           draw_cor = NULL, draw_sd = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_cor or draw_sigma
  #
  # Returns:
  #   Sigma2: The change covariance matrix.
  change_cor <- function(cor_mat, draw_cor, sparsity) {
    cor_dims <- attr(cor_mat, 'which_dims_cor')
    cor_mat_sparsity <- length(cor_dims)
    if (cor_mat_sparsity < data_dim) {
      affected_dims <- sample(cor_dims, min(sparsity, cor_mat_sparsity))
    }
    ind <- t(utils::combn(affected_dims, 2))
    change_factor <- draw_cor(nrow(ind))
    
    post_cor_mat <- cor_mat
    post_cor <- change_factor * cor_mat[ind]
    post_cor_mat[ind] <- post_cor
    post_cor_mat[ind[, c(2, 1), drop = FALSE]] <- post_cor
    
    maxit <- max(200 - data_dim, 0)
    as.matrix(Matrix::nearPD(post_cor_mat,
                             corr     = TRUE,
                             maxit    = maxit,
                             do2eigen = TRUE,
                             posd.tol = 1e-8)$mat)
  }
  
  msg <- 'ERROR: Either a variance or a correlation change distribution must be specified'
  assertthat::assert_that(!is.null(draw_cor) || !is.null(draw_sd) , msg = msg)
  
  data_dim <- ncol(cor_mat)
  change_sparsity <- length(affected_dims)
  post_cov_mat <- cor_mat
  
  # Correlation change handling
  if (!is.null(draw_cor)) {
    post_cov_mat <- change_cor(cor_mat, draw_cor, change_sparsity)
  }
  
  # Variance change handling
  if (!is.null(draw_sd)) {
    post_sd <- rep(1, data_dim)
    post_sd[affected_dims] <- draw_sd(change_sparsity)
    post_cov_mat <- diag(post_sd) %*% post_cov_mat %*% diag(post_sd)
  }
  
  post_cov_mat
}