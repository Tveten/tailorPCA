library(clusterGeneration)

generate_cor_mat <- function(N, K0 = N) {
  # k : Sparsity level, number of correlated dimensions.
  
  # set.seed(652)
  Sigma <- genPositiveDefMat(K0, covMethod = 'onion', rangeVar = c(1, 1))$Sigma
  if (K0 != N) {
    identity.mat <- diag(rep(1, N - K0))
    zero.mat <- matrix(0, ncol = N - K0, nrow = K0)
    Sigma <- cbind(Sigma, zero.mat)
    Sigma <- rbind(Sigma, cbind(t(zero.mat), identity.mat))
  }
  return(Sigma)
}

generate_cov_mat <- function(N, K0 = N, range_sd = c(0.2, 5)) {
  R <- generate_cor_mat(N, K0)
  sigma <- diag(rep(runif(N, range_sd[1], range_sd[2])))
  sigma %*% R %*% sigma
}

test_is_cor_mat <- function() {
  if (!is_cor_mat(generate_cov_mat(10))) print('Yey')
  else print('Ney')
  
  if (is_cor_mat(generate_cor_mat(10))) print('Yey')
  else print('Ney')
}
