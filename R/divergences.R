divergence_env <- new.env(hash = FALSE)

divergence_env$normal_hellinger <- function(mu1, sigma1, mu2, sigma2) {
  sqrt(1 - sqrt(2 * sigma1 * sigma2 / (sigma1^2 + sigma2^2)) *
         exp(-1/4 * (mu1 - mu2)^2 / (sigma1^2 + sigma2^2)))
}

divergence_env$normal_KL <- function(mu1, sigma1, mu2, sigma2) {
  1 / 2 * (sigma2^2 / sigma1^2 + (mu1 - mu2)^2 / sigma1 -
             1 + log(sigma1^2 / sigma2^2))
}

divergence_env$normal_bhat <- function(mu1, sigma1, mu2, sigma2) {
  - 1 / 2 * log(2 * sigma1 * sigma2 / (sigma1^2 + sigma2^2)) + 
    (mu1 - mu2)^2 / (4 * (sigma1^2 + sigma2^2))
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

#' Add a statistical divergence function.
#'
#' This function is used to specify a divergence function beyond the three
#' already implemented ones.
#' 
#' After calling this function, your implemented divergence function can used in
#' \code{\link{tpca}} by the argument \code{divergence = name}.
#' 
#' @param name A character string with the name of the divergence.
#' @param divergence A function with arguments (pre_change_mean, pre_change_sd,
#' post_change_mean, post_change_sd), that returns the measured divergence.
#' 
#' @examples
#' simple_divergence <- function(mu1, sigma1, mu2, sigma2) {
#'   (mu2 - mu1)^2 + sigma2^2 / sigma1^2
#' }
#' add_divergence('simple', simple_divergence)
#' tpca(generate_cov_mat(10), divergence = 'simple')
#
#' @export
add_divergence <- function(name, divergence) {
  divergence_env[[name]] <- divergence
}