#' Set a uniform change distribution.
#'
#' This function is used to specify a change distribution to be used in
#' conjuction with the function \code{tpca}. All components of the distribution
#' (marginal and conditional distributions) are uniform, but the
#' probability/importance of each type of change can be specified, along with
#' the range of sparsity of the change, as well as ranges for the sizes and
#' directions of the different change types. In each simulation run, after a
#' change sparsity has been drawn, which dimensions that are affected by a
#' change is always randomized. Choices for the distribution should reflect
#' prior knowledge about which changes that are of interest to detect.
#'
#' See references.
#' 
#' @param data_dim An integer specifying the dimension of the data.
#' @param prob A numeric vector of length 3 specifying the probability of a
#' change in the mean, variance or correlation, respectively.
#' @param sparsities An integer vector containing values between 2 (minimum 
#' sparsity for changes in the correlation) and data_dim.
#' @param mean_int A vector of length 2 specifying the lower and upper bound
#' of the interval that changes in the mean components are drawn uniformly from.
#' @param sd_int A vector of length 2 specifying the lower and upper bound
#' of the interval that multiplicative changes in the standard deviations
#' are drawn from. Must be between \eqn{\epsilon} and \eqn{\epsilon^{-1}}, where 
#' \eqn{\epsilon = } \code{sqrt(.Machine$double.eps)}.
#' @param sd_inc_prob A numeric between 0 and 1 indicating the probability of
#' a change in the standard deviation being an increase.
#' @param cor_int A vector of length 2 specifying the lower and upper bound
#' of the interval that multiplicative changes in the correlations are drawn 
#' from. Only values between 0 and 1 are allowed, i.e. a jump towards
#' 0 of the correlation coefficients.
#
#' @return \code{set_uniform_cd} returns an S3 object of class "change_distr". 
#' This is a list containing the following functions:
#' \describe{
#'   \item{\code{draw_types(n_sim)}}{A function that draws \code{n_sim} change 
#'   types among ('mean', 'sd', 'cor') with probabilities \code{prob}.}
#'   \item{\code{draw_sparsities(n_sim)}}{A function that draws \code{n_sim} 
#'   change sparsities uniformly from \code{sparsities}.}
#'   \item{\code{draw_dims(k)}}{A function that uniformly draws a size k subset 
#'   of affected dimensions.}
#'   \item{\code{draw_mean(k)}}{A function that draws k values uniformly from
#'   \code{mean_int}.}
#'   \item{\code{draw_sd(k)}}{A function that draws an increase in the 
#'   standard deviation uniformly from [0, \code{sd_int[2]}] with probability 
#'   \code{sd_inc_prob}, and a decrease in standard deviation uniformly
#'   from [\code{sd_int[1]}, 1] with probability 1 - \code{sd_inc_prob}.}
#'   \item{\code{draw_cor(k)}}{A function that draws k values uniformly from
#'   \code{cor_int}}
#' }
#' The functions with argument \code{n_sim} are called with the number of 
#' simulations as arguments, while the ones with argument \code{k} are called
#' in each simulation run with a number depending on the change sparsity.
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