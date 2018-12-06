context('tpca and tpca_helpers')

test_that('tpca returns sensible output', {
  N <- 10
  n_sim <- 10^2
  cor_mat <- rcov_mat(N, N/2)
  cutoffs <- c(0, 0.5, 0.8, 0.9, 0.99, 1)
  for (j in seq_along(cutoffs)) {
    tpca_obj <- tpca(cor_mat, cutoff = cutoffs[j], n_sim = n_sim)
    if (j == 1) expect_equal(length(tpca_obj$which_axes), 1)
    if (j == length(cutoffs)) expect_equal(length(tpca_obj$which_axes), N)
    expect_true(all(dim(tpca_obj$divergence_sim) == c(N, n_sim)))
    for (i in seq_along(tpca_obj)) {
      expect_true(!any(is.na(tpca_obj[[i]])))
      expect_true(!is.null(tpca_obj[[i]]))
      expect_true(!any(is.nan(tpca_obj[[i]])))
    }
  }
})

expect_identical_vectors <- function(x, y) {
  expect_true(all(x == y))
}

test_that('which_dims_cor return correct dimensions', {
  N <- 10
  K0 <- c(0, 2, 5, 10)
  cor_mats <- lapply(K0, rcor_mat, N = N)
  which_dims_list <- lapply(cor_mats, which_dims_cor)
  expected_output <- list(0, 1:2, 1:5, 1:10)
  Map(expect_identical_vectors, which_dims_list, expected_output)
})

test_that('which_axes returns correctly', {
  prop_max <- c(0.01, 0.04, 0.05, 0.1, 0.3, 0.5)
  keep_prop <- c(0, 0.5, 0.9, 0.93, 1)
  max_axes <- c(1, 2, 5)
  expect_list <- list(6, 6, 6, 
                      6, 6, 6, 
                      6, 6:5, 6:4,
                      6, 6:5, 6:3,
                      6, 6:5, 6:2)
  for (i in seq_along(keep_prop)) {
    for (j in seq_along(max_axes)) {
      expect_identical_vectors(which_axes(prop_max, keep_prop[i], max_axes[j]),
                               expect_list[[(i - 1) * length(max_axes) + j]])
    }
  }
})