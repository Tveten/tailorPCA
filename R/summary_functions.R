#' @export
summary.tpca <- function(tpca_obj) {
  cat('Probability of projection j being the most sensitive (the projections are ordered decreasingly w.r.t variance): \n')
  prop <- tpca_obj$prop_axes_max
  names(prop) <- paste0("j = ", 1:length(prop))
  print(prop)
  cat('\n')
  cat('Chosen axes/projections in decreasing order of sensitivity: \n')
  cat(tpca_obj$which_axes)
}