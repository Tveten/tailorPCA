is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is_in_interval <- function(x, interval = c(0, 1)) {
  if(!is.numeric(x)) return(FALSE)
  if (x < interval[1] || x > interval[2]) return(FALSE)
  else return(TRUE)
}

is_cor_mat <- function(cov_mat) {
  isTRUE(all.equal(diag(cov_mat), rep(1, ncol(cov_mat))))
}

is_positive_definite <- function(cov_mat, tol = 1e-8) {
  eigen_values <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigen_values < tol)) return(FALSE)
  else return(TRUE)
}