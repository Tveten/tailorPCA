change_distr_env <- new.env(hash = FALSE)

change_distr_env$full_uniform <- function() {
  
}
  
)

change_distr_env$uniform <- function()
  p <- rep(1/3, 3)  # Probability of each type of change.
  draw_k <- function(n) {
    sample(2:(N/2), n, replace = TRUE)
  }
  draw_mu <- function(n) {
    runif(n, -1.5, 1.5)
  }
  draw_sigma <- function(n) {
    change.factors <- runif(n, 1, 2.5)
    ind.decrease <- sample(c(TRUE, FALSE), n, replace = TRUE)
    change.factors[ind.decrease] <- sapply(change.factors[ind.decrease],
                                           function(x) 1/x)
    change.factors
  }
  draw_rho <- function(n) {
    runif(n, 0, 1)
  }
  affected_dims <- function(k) {
    sample(1:N, k)
  }
  
)


get_change_distr <- function(change_distr_str) {
  assertthat::assert_that(is.character(change_distr_str))
  change_distr <- change_distr_env[[change_distr_str]]
  
  msg = paste0("The supplied change distribution ('", change_distr_str,"') is not implemented.")
  assertthat::assert_that(!is.null(change_distr), msg = msg)
  
  change_distr
}

add_change_distr <- function() {
  
}