library(microbenchmark)

test_is_cor_mat <- function() {
  if (!is_cor_mat(generate_cov_mat(10))) print('Yey')
  else print('Ney')
  
  if (is_cor_mat(generate_cor_mat(10))) print('Yey')
  else print('Ney')
}

test_tpca <- function(data_dim = 10, n_sim = 10^3, change_distr = 'full_uniform') {
  cov_mat1 <- generate_cov_mat(data_dim)
  tpca(cov_mat1, n_sim = n_sim, change_distr = change_distr)
}

test_tpca_custom <- function(data_dim = 10, n_sim = 10, ...) {
  cov_mat1 <- generate_cov_mat(data_dim)
  tpca(cov_mat1, n_sim = n_sim, change_distr = set_uniform_cd(data_dim, ...))
}

test_tpca_est <- function(data_dim = 10, n_sim = 10, change_distr = 'full_uniform') {
  cov_mat1 <- generate_cor_mat(data_dim, K0 = data_dim)
  hellinger_sims <- tpca(cov_mat1, n_sim = n_sim, change_distr = change_distr)
  hellinger_est <- rowMeans(hellinger_sims)
  plot(hellinger_est, type = 'l')
}

test_change_cor <- function(data_dim = 10, sparsity = data_dim/2) {
  set.seed(10)
  cor_mat1 <- generate_cor_mat(data_dim)
  affected_dims <- sample(1:data_dim, sparsity)
  draw_cor <- get_change_distr('cor_only', data_dim)$draw_cor
  post_cor_mat <- change_cor_mat(cor_mat1, affected_dims, draw_cor = draw_cor)
  post_cor_mat - cor_mat1
}

test_change_sd <- function(data_dim = 10, sparsity = data_dim/2) {
  set.seed(10)
  cor_mat1 <- generate_cor_mat(data_dim)
  affected_dims <- sample(1:data_dim, sparsity)
  draw_sd <- get_change_distr('sd_only', data_dim)$draw_sd
  change_cor_mat(cor_mat1, affected_dims, draw_sd = draw_sd)
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
  cov_mat1 <- generate_cor_mat(data_dim, K0 = data_dim)
  hellinger_tpca <- tpca(cov_mat1, n_sim = n_sim, change_distr = 'full_uniform')
  hellinger_est_tpca <- rowMeans(hellinger_tpca)
  
  hellinger_est_old <- est_hellinger(cov_mat1, pca(cov_mat1, eigen_values = TRUE)) 
  print(hellinger_est_tpca - hellinger_est_old)
  plot(hellinger_est_tpca, type = 'l', col = 'red', ylim = c(0, 1))
  lines(hellinger_est_old, col = 'blue')
}

compare_pca <- function(data_dim) {
  # Conclusion: Use svd when all eigenvectors/eigenvalues are needed.
  
  cor_mat <- generate_cor_mat(data_dim)
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
  cor_mat <- generate_cor_mat(data_dim)
  microbenchmark::microbenchmark(
    pca_small(cor_mat, k),
    pca_large(cor_mat, k)
  )
}
