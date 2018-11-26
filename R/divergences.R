divergence_env <- new.env(hash = FALSE)

divergence_env$normal_hellinger <- 
  function(mu1, sigma1, mu2, sigma2) {
    sqrt(1 - sqrt(2 * sigma1 * sigma2 / (sigma1^2 + sigma2^2)) *
           exp(-1/4 * (mu1 - mu2)^2 / (sigma1^2 + sigma2^2)))
  }

get_divergence <- function(divergence_name) {
  assert_class_length_noNA(divergence_name, is.character, 1)
  divergence <- divergence_env[[divergence_name]]
  msg = paste0('The supplied divergence is not implemented. Use ', 
               paste0(names(divergence_env), collapse = ', '), 
               ', or add your own by using add_divergence().')
  assertthat::assert_that(!is.null(divergence), msg = msg)
  divergence
}

add_divergence <- function(name, divergence) {
  divergence_env[[name]] <- divergence
}