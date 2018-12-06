context('Handling changes to the correlation matrix')

test_that('Correct number of correlations are changed', {
  for (i in 1:100) {
    N <- 20
    K0 <- sample(2:N, 1)
    cor_mat <- rcor_mat(N, K0)
    attr(cor_mat, 'which_dims_cor') <- which_dims_cor(cor_mat)
    change_funcs <- set_uniform_cd(N)
    change_sparsity <- change_funcs$draw_sparsities(1)
    affected_dims <- change_funcs$draw_dims(change_sparsity)
    draw_cor <- change_funcs$draw_cor
    cor_mat_changed <- change_cor_mat(cor_mat, affected_dims,
                                      do_nearPD = FALSE, draw_cor = draw_cor)
    expect_equal(sum(cor_mat - cor_mat_changed != 0), 
                 2 * sum(1:(min(change_sparsity, K0) - 1)))
  }
})

test_that('Attributes of changed correlation matrix is correct', {
  N <- 20
  cor_mat <- rcor_mat(N, N - 1)
  attr(cor_mat, 'which_dims_cor') <- which_dims_cor(cor_mat)
  change_funcs <- set_uniform_cd(N)
  change_sparsity <- change_funcs$draw_sparsities(1)
  affected_dims <- change_funcs$draw_dims(change_sparsity)
  draw_cor <- change_funcs$draw_cor
  draw_sd <- change_funcs$draw_sd
  cor_mat_changed <- change_cor_mat(cor_mat, affected_dims,
                                    draw_cor = draw_cor)
  expect_identical(is.numeric(cor_mat), is.numeric(cor_mat_changed))
  expect_identical(class(cor_mat), class(cor_mat_changed))
  expect_identical(dim(cor_mat), dim(cor_mat_changed))
})

test_that('Correct standard deviations changed', {
  for (i in 1:100) {
    N <- 10
    K0 <- sample(2:N, 1)
    cor_mat <- rcor_mat(N, K0)
    attr(cor_mat, 'which_dims_cor') <- which_dims_cor(cor_mat)
    change_funcs <- set_uniform_cd(N)
    change_sparsity <- change_funcs$draw_sparsities(1)
    affected_dims <- change_funcs$draw_dims(change_sparsity)
    draw_sd <- change_funcs$draw_sd
    cor_mat_changed <- change_cor_mat(cor_mat, affected_dims, draw_sd = draw_sd)
    dims_changed <- (1:N)[which(diag(cor_mat) - diag(cor_mat_changed) != 0)]
    expect_true(all(dims_changed == sort(affected_dims)))
  }
})
