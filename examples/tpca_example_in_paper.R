library(tpca)

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
  tpca_obj <- tpca(cor_mat, 'halfsparse_uniform', n_sim = 10^4)
  singles_plot <- ggplot_singles(tpca_obj, title = '(a)')
  type_plot <- ggplot_types_mean(tpca_obj, title = '(b)')
  sparsity_plot <- ggplot_sparsity_mean(tpca_obj, title = '(c)')
  
  if (show) gridExtra::grid.arrange(singles_plot, type_plot, sparsity_plot, nrow = 1)
  gridExtra::arrangeGrob(singles_plot, type_plot, sparsity_plot, nrow = 1)
}

save_example_figure <- function() {
  save_figure(example_tpca_figure(), 'ex_hellinger-d20-3')
}

example_prop_tpca_figure <- function(show = FALSE) {
  cor_mat <- get_example_cor_mat()
  tpca_obj <- tpca(cor_mat, 'halfsparse_uniform', n_sim = 10^4)
  overall_plot <- ggplot_prop_max(tpca_obj, title = 'Overall', print_p = TRUE)
  mean_plot <- ggplot_prop_max(tpca_obj, 'mean', title = 'Mean', print_p = TRUE)
  sd_plot <- ggplot_prop_max(tpca_obj, 'sd', title = 'Variance', print_p = TRUE)
  cor_plot <- ggplot_prop_max(tpca_obj, 'cor', title = 'Correlation', print_p = TRUE)
  
  if (show) gridExtra::grid.arrange(overall_plot, mean_plot, sd_plot, cor_plot, nrow = 1)
  invisible(gridExtra::arrangeGrob(overall_plot, mean_plot, sd_plot, cor_plot, nrow = 1))
}

save_example_prop_figure <- function() {
  figure <- example_prop_tpca_figure()
  name <- 'ex_hellinger_prop-d20'
  save_figure(figure, name, base_width = 2.16, base_height = 2.6)
}