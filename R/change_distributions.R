set_uniform_cd <- function(prob        = rep(1/3, 3), 
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
  # msg = 
  # assertthat::assert_that(, msg = msg)
  
  sparsities_expr <- substitute(sparsities)
  
  function(data_dim) {
    # TODO: Error handling of sparsities.
    sparsities <- eval(sparsities_expr)
    list(
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
  }
}

change_distr_env <- new.env(hash = FALSE)

change_distr_env$full_uniform <- set_uniform_cd(prob        = rep(1/3, 3), 
                                                sparsities  = 2:data_dim,
                                                mean_int    = c(-1.5, 1.5), 
                                                sd_int      = c(2.5^(-1), 2.5),
                                                sd_inc_prob = 0.5,
                                                cor_int     = c(0, 1))

change_distr_env$mean_only <- 
  set_uniform_cd(prob        = c(1, 0, 0), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))

change_distr_env$sd_only <- 
  set_uniform_cd(prob        = c(0, 1, 0), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))

change_distr_env$cor_only <- 
  set_uniform_cd(prob        = c(0, 0, 1), 
                 sparsities  = 2:data_dim,
                 mean_int    = c(-1.5, 1.5), 
                 sd_int      = c(2.5^(-1), 2.5),
                 sd_inc_prob = 0.5,
                 cor_int     = c(0, 1))

get_change_distr <- function(change_distr_str, data_dim) {
  assertthat::assert_that(is.character(change_distr_str))
  change_distr <- change_distr_env[[change_distr_str]](data_dim)
  
  msg = paste0("The supplied change distribution ('", change_distr_str,"') is not implemented.")
  assertthat::assert_that(!is.null(change_distr), msg = msg)
  
  change_distr
}

add_change_distr <- function() {
  
}


# print(environment())
# foo2 <- function(a = b) {
#   a_expr <- substitute(a)
#   function(b) {
#     print(eval(a_expr))
#   }
# }
# 
# foo <- function(a = b) {
#   if (!exists('b'))
#   default_env <- new.env(hash = FALSE)
#   delayedAssign('a', b, eval.env = default_env)
#   function(b) {
#     default_env$b <- b
#     print(a)
#   }
# }
