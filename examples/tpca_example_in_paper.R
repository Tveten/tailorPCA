library(tpca)

get_example_cor_mat <- function() {
  d <- 20
  k0 <- d / 2
  set.seed(7)
  rcor_mat(d, k0)
}

example_tpca_figure <- function(show = FALSE) {
  cor_mat <- get_example_cor_mat()
  tpca_obj <- tpca(cor_mat, 'halfsparse_uniform', n_sim = 10^4)
  singles_plot <- ggplot_singles(tpca_obj)
  type_plot <- ggplot_types_mean(tpca_obj)
  sparsity_plot <- ggplot_sparsity_mean(tpca_obj)
  
  if (show) gridExtra::grid.arrange(singles_plot, type_plot, sparsity_plot, nrow = 1)
  gridExtra::arrangeGrob(singles_plot, type_plot, sparsity_plot, nrow = 1)
}

save_example_figure <- function() {
  save_figure(example_tpca_figure(), 'ex_hellinger-d20-k010-3')
}
