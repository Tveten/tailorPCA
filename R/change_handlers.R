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
                           draw_rho = NULL, draw_sigma = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_rho or draw_sigma
  #
  # Returns:
  #   Sigma2: The change covariance matrix.
  change_cor <- function(cor_mat, draw_rho, n.affected) {
    affected_cor_dims <- sample(1:cor_mat_sparsity, min(n.affected, cor_mat_sparsity))
    ind <- combn(affected_cor_dims, 2)
    change_factor <- draw_rho(sum(1:length(ind)))
    post_cor_mat <- cor_mat
    for (i in 1:ncol(ind)) {
      post_cor_mat[ind[1, i], ind[2, i]] <- change_factor[i] * post_cor_mat[ind[1, i], ind[2, i]]
      post_cor_mat[ind[2, i], ind[1, i]] <- change.factor[i] * post_cor_mat[ind[2, i], ind[1, i]]
    }
    post_cor_mat <- as.matrix(nearPD(post_cor_mat, 
                                     corr     = TRUE, 
                                     maxit    = 20,
                                     do2eigen = TRUE,
                                     posd.tol = 1e-8)$mat)
  }
  
  msg <- 'ERROR: Either a variance or a correlation change distribution must be specified'
  assertthat::assert_that(all(is.null(c(draw_rho, draw_sigma))), msg = msg)
  
  data_dim <- ncol(cor_mat)
  cor_mat_sparsity <- sum(cor_mat[, 1] != 0)
  change_sparsity <- length(affected.dims)
  post_cov_mat <- cor_mat
  
  # Correlation change handling
  if (!is.null(draw_rho)) {
    ind <- 1:cor_mat_sparsity
    sub_cor_mat <- cor_mat[ind, ind, drop = FALSE]
    post_cov_mat[ind, ind] <- change_cor(sub_cor_mat, draw_rho, change_sparsity)
  }
  
  # Variance change handling
  if (!is.null(draw_sigma)) {
    post_sd <- rep(1, data_dim)
    post_sd[affected_dims] <- draw_sigma(change_sparsity)
    post_cov_mat <- diag(post_sd) %*% post_cov_mat %*% diag(post_sd)
  }
  
  post_cov_mat
}


