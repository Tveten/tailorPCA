assert_cov_mat <- function(cov_mat) {
  cov_mat_str <- deparse(substitute(cov_mat))
  na_msg <- paste0(cov_mat_str, ' cannot contain NAs.')
  assertthat::assert_that(!any(is.na(cov_mat)), msg = na_msg)
  numeric_msg <- paste0(cov_mat_str, ' must be numeric.')
  assertthat::assert_that(is.numeric(cov_mat), msg = numeric_msg)
  matrix_msg <- paste0(cov_mat_str, ' must be of class "matrix".')
  assertthat::assert_that(class(cov_mat) == 'matrix', msg = matrix_msg)
  symmetric_msg <- paste0(cov_mat_str, ' is not a symmetric matrix.')
  assertthat::assert_that(isSymmetric(cov_mat), msg = symmetric_msg)
  posdef_msg <- paste0(cov_mat_str, ' is not a positive definite matrix (some eigenvalues are < 1e-8).')
  assertthat::assert_that(is_positive_definite(cov_mat), msg = posdef_msg)
}

assert_natural_number <- function(n) {
  msg <- paste0(n, 'must be an integer larger than 0.')
  assertthat::assert_that(is_whole_number(n), n > 0, msg = msg)
}

assert_in_interval <- function(x, interval) {
  interval_msg <- paste0(x, ' must be a numeric between ', interval[1], 
                         ' and ', interval[2], '.')
  assertthat::assert_that(is_in_interval(x, interval), msg = interval_msg)
}

assert_class_length_noNA <- function(x, is_class, l) {
  length_msg <- paste0(x, 'must have length ', l, '.')
  assertthat::assert_that(length(x) == l, msg = length_msg)
  na_msg <- paste0(x, 'cannot be NA.')
  assertthat::assert_that(!is.na(x), msg = na_msg)
  is_class_str <- deparse(substitute(is_class))
  class_msg <- paste0(x, ' must be ', is_class_str)
  assertthat::assert_that(is_class(x), msg = class_msg)
}

is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is_in_interval <- function(x, interval) {
  # Input:
  #   x: A single numeric.
  x >= interval[1] && x <= interval[2]
}

is_interval <- function(x) {
  if (length(x == 2)) {
    return(x[1] <= x[2])
  } else return(FALSE)
}

is_cor_mat <- function(cov_mat) {
  isTRUE(all.equal(diag(cov_mat), rep(1, ncol(cov_mat))))
}

is_positive_definite <- function(cov_mat, tol = 1e-8) {
  eigen_values <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  all(eigen_values >= tol)
}

is_prob <- function(p) {
  # p: A vector  with probability weights.
  condition1 <- all(vapply(p, is_in_interval, logical(1), c(0, 1)))
  condition2 <- isTRUE(all.equal(sum(p), 1))
  condition1 && condition2
}

interval_msg <- function(obj) {
  obj_name <- deparse(substitute(obj))
  paste0(obj_name, ' must be a vector of length 2 where the first element is smaller than or equal two the second.')
}