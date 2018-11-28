library(testthat)
context('Handling changes to the correlation matrix')

for (i in 1:100) {
  N <- 20
  K0 <- sample(2:N, 1)
  cor_mat <- generate_cor_mat(N, K0)
  attr(cor_mat, 'which_dims_cor') <- which_dims_cor(cor_mat)
  change_funcs <- set_uniform_cd(N)
  change_sparsity <- change_funcs$draw_sparsities(1)
  affected_dims <- change_funcs$draw_dims(change_sparsity)
  draw_cor <- change_funcs$draw_cor
  cor_mat_changed <- change_cor_mat(cor_mat, affected_dims,
                                    do_nearPD = FALSE, draw_cor = draw_cor)
  test_that('Correct number of correlations are changed', {
    expect_equal(sum(cor_mat - cor_mat_changed != 0), 
                 2 * sum(1:(min(change_sparsity, K0) - 1)))
  })
}
