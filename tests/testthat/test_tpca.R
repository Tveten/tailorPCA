context('tpca')

test_that('tpca returns sensible output', {
  N <- 10
  n_sim <- 10^2
  cor_mat <- generate_cov_mat(N, N/2)
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