is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is_in_interval <- function(x, interval = c(0, 1)) {
  # Input:
  #   x: A single numeric.
  if(!is.numeric(x)) return(FALSE)
  if (x < interval[1] || x > interval[2]) return(FALSE)
  else return(TRUE)
}

is_interval <- function(x) {
  condition1 <- length(x == 2)
  if (length(x == 2)) {
    return(x[1] <= x[2])
  } else return(FALSE)
}

is_cor_mat <- function(cov_mat) {
  isTRUE(all.equal(diag(cov_mat), rep(1, ncol(cov_mat))))
}

is_positive_definite <- function(cov_mat, tol = 1e-8) {
  eigen_values <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigen_values < tol)) return(FALSE)
  else return(TRUE)
}

is_prob <- function(p) {
  # p: A vector  with probability weights.
  condition1 <- all(vapply(p, is_in_interval, logical(1), c(0, 1)))
  condition2 <- isTRUE(all.equal(sum(p), 1))
  if (condition1 && condition2) return(TRUE)
  else return(FALSE)
}

interval_msg <- function(obj) {
  obj_name <- deparse(substitute(obj))
  paste0(obj_name, ' must be a vector of length 2 where the first element is smaller than or equal two the second.')
}