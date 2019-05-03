#' @export
get_example_cor_mat <- function() {
  d <- 20
  set.seed(7)
  rcor_mat(d)
}

p_cor_greater_than <- function(a, cor_mat) {
  d <- ncol(cor_mat)
  correlations <- cor_mat - diag(rep(1, d))
  n_correlations <- sum(correlations != 0)
  sum(abs(correlations) >= a) / n_correlations
}

example_tpca_figure <- function(show = FALSE) {
  cor_mat <- get_example_cor_mat()
  tpca_obj <- tpca(cor_mat, 'semisparse_uniform', n_sim = 10^4)
  type_plot <- ggplot_types(tpca_obj, title = '(a)')
  sparsity_plot <- ggplot_sparsities(tpca_obj, title = '(b)')
  singles_plot <- ggplot_singles(tpca_obj, title = '(c)')
  
  if (show) gridExtra::grid.arrange(type_plot, sparsity_plot, singles_plot, nrow = 1)
  gridExtra::arrangeGrob(type_plot, sparsity_plot, singles_plot, nrow = 1)
}

save_example_figure <- function() {
  save_figure(example_tpca_figure(), 'ex_hellinger-d20')
}

example_prop_tpca_figure <- function(show = FALSE) {
  cor_mat <- get_example_cor_mat()
  tpca_obj <- tpca(cor_mat, 'semisparse_uniform', n_sim = 10^4)
  overall_plot <- ggplot_prop(tpca_obj, title = 'Overall')
  mean_plot <- ggplot_prop(tpca_obj, 'mean', title = 'Mean')
  sd_plot <- ggplot_prop(tpca_obj, 'sd', title = 'Variance')
  cor_plot <- ggplot_prop(tpca_obj, 'cor', title = 'Correlation')
  
  if (show) gridExtra::grid.arrange(overall_plot, mean_plot, sd_plot, cor_plot, nrow = 1)
  invisible(gridExtra::arrangeGrob(overall_plot, mean_plot, sd_plot, cor_plot, nrow = 1))
}

save_example_prop_figure <- function() {
  figure <- example_prop_tpca_figure()
  name <- 'ex_hellinger_prop-d20'
  save_figure(figure, name, base_width = 2.16, base_height = 2.6)
}

#' Reproduces the figures in Section 2
#' 
#' This function does the following:
#' 1) Runs the tpca-algorithm on the randomly drawn correlation 
#' matrix with 10^4 change simulations. 
#' (Call \code{get_example_cor_mat()} to get the correlation matrix.)
#' 2) Creates ggplot objects for the sensitivity results and arranges them.
#' 3) Saves two files to your current directory: ex_hellinger-d20.png and
#' ex_hellinger_prop-d20.png.
#' 
#' @export
reproduce_tpca_example <- function() {
  save_example_figure()
  save_example_prop_figure()
}