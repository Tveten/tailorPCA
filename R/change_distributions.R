set_uniform_cd <- function(data_dim,
                           prob        = rep(1/3, 3), 
                           sparsities  = 2:data_dim,
                           mean_int    = c(-1.5, 1.5), 
                           sd_int      = c(2.5^(-1), 2.5), 
                           sd_inc_prob = 0.5,
                           cor_int     = c(0, 1)) {
  # Input:
  #   data_dim:    Dimensionality of the data.
  #   prob:        Vector with probabilities of each change type (mean, sd, cor)
  #   sparsities:  Vector with the relevant change sparsities, e.g. 2:data_dim
  #   mean_int:    Lower and upper bound for a change in mean.
  #   sd_int:      Lower and upper bound for a change in standard deviation.
  #   sd_inc_prob: The probability of 
  #   cor_int:     Lower and upper bound for a change in correlation. Must be 
  #                between 0 and 1.
  
  # TODO: A lot of error handling.
  ## ERROR HANDLING ------------------------------------------------------------
  #  prob:
  assertthat::assert_that(is.numeric(prob), !any(is.na(prob)))
  prob_length_msg = paste0('prob =', prob, 'has to be a vector of length 3.')
  assertthat::assert_that(length(prob) == 3, msg = prob_msg)
  prob_msg = paste0('prob =', prob, 
                    'is not a probability (summing to one and elements between 0 and 1).')
  assertthat::assert_that(is_prob(prob), msg = prob_msg)
  
  #  sparsities:
  sparsity.msg <- 'sparsities must be an integer vector between 2 (minimum sparsity for changes in the correlation) and data_dim.'
  assertthat::assert_that(is.numeric(sparsities), !any(is.na(sparsities)))
  assertthat::assert_that(all(is_whole_number(sparsities)), 
                          all(sparsities >= 2),
                          all(sparsities <= data_dim),
                          msg = sparsity.msg)
  
  #  mean_int:
  assertthat::assert_that(is.numeric(mean_int), !any(is.na(mean_int)))
  assertthat::assert_that(is_interval(mean_int), msg = interval_msg(mean_int))
  
  #  sd_int:
  assertthat::assert_that(is.numeric(sd_int), !any(is.na(sd_int)))
  assertthat::assert_that(is_interval(sd_int), msg = interval_msg(sd_int))
  
  #  sd_inc_prob:
  assertthat::assert_that(is.numeric(sd_inc_prob), !is.na(sd_inc_prob))
  assertthat::assert_that(sd_inc_prob >= 0, sd_inc_prob <= 1,
                          msg = 'sd_inc_prob must be a value between 0 and 1.')
  
  #  cor_int:
  assertthat::assert_that(is.numeric(cor_int), !any(is.na(cor_int)))
  assertthat::assert_that(is_interval(cor_int), msg = interval_msg(cor_int))
  
  ## MAIN ----------------------------------------------------------------------
  return_list <- list(
  'draw_types'      = function(n_sim) {
                        change_types <- c('mean', 'sd', 'cor')
                        sample(change_types, n_sim, prob = prob, replace = TRUE)
                      },
  'draw_sparsities' = function(n_sim) sample(sparsities, n_sim, replace = TRUE),
  'draw_dims'       = function(n) sample(1:data_dim, n),
  'draw_mean'       = function(n) runif(n, mean_int[1], mean_int[2]),
  'draw_sd'         = function(n) {
                        sd_increases <- runif(round(sd_inc_prob * n), 1, sd_int[2])
                        sd_decreases <- runif(n - length(sd_increases), sd_int[1], 1)
                        sample(c(sd_increases, sd_decreases), n)
                      },
  'draw_cor'        = function(n) runif(n, cor_int[1], cor_int[2])
  )
  structure(return_list, class = 'change_distr')
}

change_distr_env <- new.env(hash = FALSE)

change_distr_env$full_uniform <- function(data_dim) {
  set_uniform_cd(data_dim,
                 prob        = rep(1/3, 3), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))
}

change_distr_env$mean_only <- function(data_dim) {
  set_uniform_cd(data_dim,
                 prob        = c(1, 0, 0), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))
}

change_distr_env$sd_only <- function(data_dim) {
  set_uniform_cd(data_dim,
                 prob        = c(0, 1, 0), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))
}

change_distr_env$cor_only <- function(data_dim) {
  set_uniform_cd(data_dim,
                 prob        = c(0, 0, 1), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))
}

get_change_distr <- function(change_distr_str, data_dim) {
  assertthat::assert_that(is.character(change_distr_str))
  
  change_distr <- change_distr_env[[change_distr_str]]
  # msg = paste0('The supplied change distribution ("', change_distr_str,'") is not implemented. Use ', names(change_distr_env), ', or make your own by using set_uniform_cd()')
  msg = paste0('The supplied change distribution is not implemented. Use ', 
               paste0(names(change_distr_env), collapse = ', '), 
               ', or make your own by using set_uniform_cd().')
  assertthat::assert_that(!is.null(change_distr), msg = msg)
  
  change_distr(data_dim)
}