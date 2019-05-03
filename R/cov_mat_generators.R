#' Generate a random correlation matrix
#' 
#' A wrapper for \code{\link{clusterGeneration::rcorrmatrix}}, but with the possibility of
#' setting some dimensions to being independent of the rest. 
#' 
#' @param d An integer specifying the dimension of the correlation matrix.
#' @param k0 An integer. d - k0 are the number of independent dimensions.
#' @param alphad A positive real number. See \code{\link{clusterGeneration::rcorrmatrix}}.
#' 
#' @return A correlation matrix with an attribute 'which_dims_cor', indicating 
#' which dimensions that are the correlated ones.
#' 
#' @export
rcor_mat <- function(d, k0 = d, alphad = 1) {
  # K0: Sparsity level, number of correlated dimensions.
  
  if (k0 == 0) return(diag(rep(1, d)))
  
  Sigma <- clusterGeneration::rcorrmatrix(k0, alphad = alphad)
  if (k0 != d) {
    identity.mat <- diag(rep(1, d - k0))
    zero.mat <- matrix(0, ncol = d - k0, nrow = k0)
    Sigma <- cbind(Sigma, zero.mat)
    Sigma <- rbind(Sigma, cbind(t(zero.mat), identity.mat))
  }
  structure(Sigma, 'which_dims_cor' = 1:k0)
}

#' Generate a random covariance matrix
#' 
#' Generates a random correlation matrix by \eqn{R} \code{\link{rcor_mat}}, then
#' scales it to a covariance matrix by \eqn{CRC}, where
#' C is a diagonal matrix of uniformly generated standard deviations in
#' a specified range.
#' 
#' @param d An integer specifying the dimension of the correlation matrix.
#' @param k0 An integer. d - k0 are the number of independent dimensions.
#' @param range_sd A vector of length 2 specifying the lower and upper bound
#' of the uniform distribution for generating standard deviations.
#' @param ... Other arguments to \code{\link{rcor_mat}}.
#' 
#' @return A covariance matrix with an attribute 'which_dims_cor', indicating 
#' which dimensions that are the correlated ones.
#' 
#' @export
rcov_mat <- function(d, k0 = d, range_sd = c(0.2, 5), ...) {
  R <- rcor_mat(d, k0, ...)
  sigma <- diag(rep(runif(d, range_sd[1], range_sd[2])), nrow = d)
  structure(sigma %*% R %*% sigma, 'which_dims_cor' = 1:k0)
}

#' Generate a random correlation matrix estimate
#' 
#' Generates a random correlation matrix by \eqn{R} \code{\link{rcor_mat}} 
#' before drawing a Wishart matrix with parameters \eqn{R} and n.
#' Lastly, the Wishart matrix is standardized to a correlation matrix.
#' 
#' @param d An integer specifying the dimension of the correlation matrix.
#' @param k0 An integer. d - k0 are the number of independent dimensions.
#' @param n An integer specifying the number of samples the estimate is based on.
#' @param ... Other arguments to \code{\link{rcor_mat}}.
#' 
#' @return A correlation matrix estimate with an attribute 'n_obs', indicating 
#' how many samples the estimate is based on.
#' 
#' @export
rcor_mat_est <- function(d, k0 = d, n = 2 * d, ...) {
  # d:  data dimension
  # k0: d - k0 is the number of completely independent dimensions.
  # n:  Number of observations that the estimated cor_mat is based on.
  Sigma <- rcor_mat(d, k0, ...)
  Sigma_est <- 1 / (n - 1) * rWishart(1, n, Sigma)[, , 1]
  structure(standardize_cov_mat(Sigma_est), 'n_obs' = n)
}
