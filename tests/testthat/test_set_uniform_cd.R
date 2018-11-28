context('Basic use of a change distribution.')

# data_dim    <- 50
# prob        <-  c(0.2, 0.5, 0.3)
# sparsities  <-  2:(data_dim / 2)
# mean_int    <-  c(-1, 1)
# sd_int      <-  c(3^(-1), 2)
# sd_inc_prob <-  0.8
# cor_int     <-  c(0.25, 0.75)
# change_distr <- set_uniform_cd(data_dim, prob, sparsities, mean_int, sd_int,
#                                sd_inc_prob, cor_int)

test_that('draw_types works as intended', {
  data_dim <- 50
  change_distr <- set_uniform_cd(data_dim)
  expect_equal(length(change_distr$draw_types(3)), 3)
  expect_equal(length(change_distr$draw_types(100)), 100)
  expect_true(!any(set_uniform_cd(data_dim, prob = c(0, 0.5, 0.5))$draw_types(100) == 'mean'))
  expect_true(!any(set_uniform_cd(data_dim, prob = c(0.5, 0, 0.5))$draw_types(100) == 'sd'))
  expect_true(!any(set_uniform_cd(data_dim, prob = c(0.5, 0.5, 0))$draw_types(100) == 'cor'))
})

test_that('Illegal inputs for data_dim generates errors', {
  expect_error(set_uniform_cd(-2))
  expect_error(set_uniform_cd(NA))
  expect_error(set_uniform_cd(NULL))
  expect_error(set_uniform_cd(NaN))
  expect_error(set_uniform_cd(Inf))
  expect_error(set_uniform_cd(1:5))
  expect_error(set_uniform_cd('a'))
  expect_error(set_uniform_cd(TRUE))
  expect_error(set_uniform_cd(1.5))
})

test_that('Illegal inputs for prob generates errors', {
  expect_error(set_uniform_cd(10, prob = 0.2))
  expect_error(set_uniform_cd(10, prob = c(0.5, 0.5)))
  expect_error(set_uniform_cd(10, prob = c(0.5, 0.2, 0.5)))
  expect_error(set_uniform_cd(10, prob = c(-0.1, 0.1, 1)))
  expect_error(set_uniform_cd(10, prob = c(1.1, -0.1, 0)))
  expect_error(set_uniform_cd(10, prob = c(0.5, 0.5, 0, 0)))
  expect_error(set_uniform_cd(10, prob = as.character(c(0.5, 0.5, 0))))
})

test_that('Illegal inputs for sparsity generates errors', {
  expect_error(set_uniform_cd(10, sparsities = 2:11))
  expect_error(set_uniform_cd(10, sparsities = 2:11))
  expect_error(set_uniform_cd(10, sparsities = 1:10))
  expect_error(set_uniform_cd(10, sparsities = rep(1/10, 10)))
  expect_error(set_uniform_cd(10, sparsities = c(NA, 2:9)))
  expect_error(set_uniform_cd(10, sparsities = c(Inf, 2:9)))
})

test_that('Illegal inputs for mean_int generates errors', {
  expect_error(set_uniform_cd(10, mean_int = 2))
  expect_error(set_uniform_cd(10, mean_int = c(NA, 3)))
  expect_error(set_uniform_cd(10, mean_int = c(NaN, 3)))
  expect_error(set_uniform_cd(10, mean_int = 1:3))
  expect_error(set_uniform_cd(10, mean_int = as.character(1:2)))
})

test_that('Illegal inputs for sd_int generates errors', {
  expect_error(set_uniform_cd(10, sd_int = 2))
  expect_error(set_uniform_cd(10, sd_int = c(NA, 3)))
  expect_error(set_uniform_cd(10, sd_int = c(NaN, 3)))
  expect_error(set_uniform_cd(10, sd_int = 1:3))
  expect_error(set_uniform_cd(10, sd_int = as.character(1:2)))
  expect_error(set_uniform_cd(10, sd_int = c(-1, 3)))
})

test_that('Illegal inputs for sd_inc_prob generates errors', {
  expect_error(set_uniform_cd(10, sd_inc_prob = 2))
  expect_error(set_uniform_cd(10, sd_inc_prob = NA))
  expect_error(set_uniform_cd(10, sd_inc_prob = NULL))
  expect_error(set_uniform_cd(10, sd_inc_prob = Inf))
  expect_error(set_uniform_cd(10, sd_inc_prob = NaN))
  expect_error(set_uniform_cd(10, sd_inc_prob = 'a'))
  expect_error(set_uniform_cd(10, sd_inc_prob = c(0.5, 0.5)))
})

test_that('Illegal inputs for cor_int generates errors', {
  expect_error(set_uniform_cd(10, cor_int = c(0, 1.5)))
  expect_error(set_uniform_cd(10, cor_int = c(-0.5, 0.5)))
  expect_error(set_uniform_cd(10, cor_int = c(NaN, 1)))
  expect_error(set_uniform_cd(10, cor_int = c(NULL, 1)))
  expect_error(set_uniform_cd(10, cor_int = c(NA, 1)))
  expect_error(set_uniform_cd(10, cor_int = 1:3))
})
