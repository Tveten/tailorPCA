library(microbenchmark)

test_is_cor_mat <- function() {
  if (!is_cor_mat(rcov_mat(10))) print('Yey')
  else print('Ney')
  
  if (is_cor_mat(rcor_mat(10))) print('Yey')
  else print('Ney')
}

test_tpca <- function(data_dim = 10, n_sim = 10^3, change_distr = 'full_uniform') {
  cov_mat1 <- rcov_mat(data_dim)
  tpca(cov_mat1, n_sim = n_sim, change_distr = change_distr)
}

test_tpca_custom <- function(data_dim = 10, n_sim = 10, ...) {
  cov_mat1 <- rcov_mat(data_dim)
  tpca(cov_mat1, n_sim = n_sim, change_distr = set_uniform_cd(data_dim, ...))
}

test_tpca_est <- function(data_dim = 10, n_sim = 10, change_distr = 'full_uniform') {
  cov_mat1 <- rcor_mat(data_dim, K0 = data_dim)
  hellinger_sims <- tpca(cov_mat1, n_sim = n_sim, change_distr = change_distr)
  hellinger_est <- rowMeans(hellinger_sims)
  plot(hellinger_est, type = 'l')
}

test_change_cor <- function(data_dim = 10, sparsity = data_dim/2) {
  set.seed(10)
  cor_mat1 <- rcor_mat(data_dim)
  affected_dims <- sample(1:data_dim, sparsity)
  draw_cor <- get_change_distr('cor_only', data_dim)$draw_cor
  post_cor_mat <- change_cor_mat(cor_mat1, affected_dims, draw_cor = draw_cor)
  post_cor_mat - cor_mat1
}

test_change_sd <- function(data_dim = 10, sparsity = data_dim/2) {
  set.seed(10)
  cor_mat1 <- rcor_mat(data_dim)
  affected_dims <- sample(1:data_dim, sparsity)
  draw_sd <- get_change_distr('sd_only', data_dim)$draw_sd
  change_cor_mat(cor_mat1, affected_dims, draw_sd = draw_sd)
}

test_duplicate_diag_block <- function() {
  cov_mat <- matrix(1:16, nrow = 4, ncol = 4)
  block_dim <- 2
  duplicate_diag_block(cov_mat, block_dim)
}

test_change_cormat_tdpca <- function(unlagged_dim = 2, change_type = 'sd') {
  # set.seed(10)
  lag <- 1
  data_dim <- (lag + 1) * unlagged_dim
  cor_mat <- rcor_mat(unlagged_dim)
  x <- gen_norm_data(3 * data_dim, rep(0, unlagged_dim), cor_mat)
  x_lagged <- lag_extend(x, lag)
  cor_mat_lagged <- 1 / (ncol(x_lagged) - 1) * x_lagged %*% t(x_lagged)
  cor_mat_lagged <- standardize_cov_mat(cor_mat_lagged)
  change_sparsity <- 2
  affected_dims <- sample(1:unlagged_dim, change_sparsity)
  print(affected_dims)
  if (change_type == 'sd') {
    draw_sd <- get_change_distr('sd_only', unlagged_dim)$draw_sd
    post_cor_mat <- change_cor_mat_tdpca(cor_mat_lagged, lag, affected_dims, draw_sd = draw_sd)
  } else if (change_type == 'cor') {
    draw_cor <- get_change_distr('cor_only', unlagged_dim)$draw_cor
    post_cor_mat <- change_cor_mat_tdpca(cor_mat_lagged, lag, affected_dims, draw_cor = draw_cor)
  }
  # print(cor_mat_lagged)
  # print(post_cor_mat)
  post_cor_mat / cor_mat_lagged
}

test_tdpca <- function(unlagged_dim = 3, n_sim = 10^2, change_distr = 'full_uniform') {
  lag <- 3
  data_dim <- (lag + 1) * unlagged_dim
  cor_mat <- rcor_mat(unlagged_dim)
  x <- gen_norm_data(3 * data_dim, rep(0, unlagged_dim), cor_mat)
  x_lagged <- lag_extend(x, lag)
  cor_mat_lagged <- 1 / (ncol(x_lagged) - 1) * x_lagged %*% t(x_lagged)
  eigen(cor_mat_lagged)
  tdpca_obj <- tdpca(cor_mat_lagged, lag = lag, n_sim = n_sim, change_distr = change_distr)
  tdpca_obj$prop_axes_max
}

gen_norm_data <- function(n, mu, Sigma) {
  t(MASS::mvrnorm(n, mu = mu, Sigma = Sigma))
}

lag_extend <- function(x, lag) {
  # Matrix: (x_{t - lag}, ...,  x_{t})^T, I.e., time t is in the lower rows.
  n <- ncol(x)
  p <- nrow(x)
  x_extended <- matrix(0, nrow = p * (lag + 1), ncol = n - lag)
  vapply((lag + 1):n, function(i) as.vector(x[, (i - lag):i]),
         numeric(p * (lag + 1)))
}

benchmark_change_func <- function() {
  microbenchmark(
    test_change_sd(100),
    test_change_cor(100)
  )
}

benchmark_tpca <- function() {
  microbenchmark(test_tpca(100, 1000))
}

compare_tpca <- function(data_dim = 10, n_sim = 10^3) {
  set.seed(20)
  cov_mat1 <- rcor_mat(data_dim, K0 = data_dim)
  hellinger_tpca <- tpca(cov_mat1, n_sim = n_sim, change_distr = 'full_uniform')
  hellinger_est_tpca <- rowMeans(hellinger_tpca)
  
  hellinger_est_old <- est_hellinger(cov_mat1, pca(cov_mat1, eigen_values = TRUE)) 
  print(hellinger_est_tpca - hellinger_est_old)
  plot(hellinger_est_tpca, type = 'l', col = 'red', ylim = c(0, 1))
  lines(hellinger_est_old, col = 'blue')
}

compare_pca <- function(data_dim) {
  # Conclusion: Use svd when all eigenvectors/eigenvalues are needed.
  
  cor_mat <- rcor_mat(data_dim)
  pca <- pca(cor_mat, axes = 1:(data_dim - 1))
  pca_small <- pca_small(cor_mat, data_dim - 1)
  pca_large <- pca_large(cor_mat, data_dim - 1)
  pca_irlba <- pca_large2(cor_mat, data_dim - 1)
  print(pca$values - pca_small$values)
  print(pca$values - pca_large$values)
  # Unstable for small singular values.
  print(pca$values - pca_irlba$values)
}

benchmark_RSpectra <- function(data_dim, k) {
  # Conclusion: RSpectra faster than irlba (2x for data_dim = 100 and k = 10).
  cor_mat <- rcor_mat(data_dim)
  microbenchmark::microbenchmark(
    pca_small(cor_mat, k),
    pca_large(cor_mat, k)
  )
}

test_obs_needed <- function(data_dim, n) {
  max(sapply(1:1000, function(i) {
    cor_mat_est <- rcor_mat_est(data_dim, n = n)
    obs_needed(cor_mat_est, 0.01)
  }))
}
