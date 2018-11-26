#' Title
#'
#' Description
#'
#' Detailed description.
#'
#' @param data_dim
#' @param prob
#' @param sparsities An integer vector between 2 (minimum sparsity for changes in the correlation) and data_dim.
#' @param mean_int
#' @param sd_int
#' @param sd_inc_prob
#' @param cor_int 
#'
#' @return
#'
#' @examples 
#' 
#' @export

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
  
  ## ERROR HANDLING ------------------------------------------------------------
  assert_class_length_noNA(prob, is.numeric, 3)
  assert_prob(prob)
  
  assert_class_length_noNA(sparsities, is.numeric)
  assert_integer_in_interval(sparsities, c(2, data_dim))
  
  assert_class_length_noNA(mean_int, is.numeric, 2)
  assert_interval(mean_int)
  
  assert_class_length_noNA(sd_int, is.numeric, 2)
  assert_interval(sd_int)
  
  assert_class_length_noNA(sd_inc_prob, is.numeric, 1)
  assert_in_interval(sd_inc_prob, c(0, 1))
  
  assert_class_length_noNA(cor_int, is.numeric, 2)
  assert_interval(cor_int)
  
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

get_change_distr <- function(change_distr, data_dim) {
  if (class(change_distr) == 'change_distr') 
    return(change_distr)
  else if (is.character(change_distr)) {
    assert_class_length_noNA(change_distr, is.character, 1)
    change_distr_func <- change_distr_env[[change_distr_str]]
    msg = paste0('The supplied change distribution is not implemented. Use ', 
                 paste0(names(change_distr_env), collapse = ', '), 
                 ', or make your own by using the function set_uniform_cd.')
    assertthat::assert_that(!is.null(change_funcs), msg = msg)
    return(change_distr_func(data_dim))
  } else 
    stop(paste0('change_distr must either be a character string or belong to class change_distr. See set_uniform_cd.'))
    
}