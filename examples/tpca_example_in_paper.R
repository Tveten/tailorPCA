library(tpca)

get_example_cor_mat <- function() {
  d <- 20
  k0 <- d / 2
  set.seed(7)
  rcor_mat(d, k0)
}

example_tpca_figure <- function(save = FALSE) {
  cor_mat <- get_example_cor_mat()
  tpca_obj <- tpca(cor_mat, 'halfsparse_uniform', n_sim = 10^4)
  singles_plot <- ggplot_singles(tpca_obj)
  type_plot <- ggplot_types_mean(tpca_obj)
  sparsity_plot <- ggplot_sparsity_mean(tpca_obj)
  gridExtra::grid.arrange(singles_plot, type_plot, sparsity_plot, nrow = 1)
  all_plots <- gridExtra::arrangeGrob(singles_plot, type_plot, 
                                      sparsity_plot, nrow = 1)
}

save_figure <- function(ggplot_obj, name, extension = 'png') {
  plot_grid <- dim(ggplot_obj)
  A4_width <- 8.27  # inches
  A4_height <- 11.69  # inches
  base_width <- 2.88
  base_height <- 2.3
  if (is.null(plot_grid)) {
    width = base_width
    height = base_height
  } else {
    width = plot_grid[2] * base_width
    height = plot_grid[1] * base_height
  }
  dir <- '/mn/sarpanitu/ansatte-u6/martintv/Documents/tpca paper/Text/Fig'
  ggplot2::ggsave(paste0(dir, '/', name, '.', extension), ggplot_obj,
                  width = width, height = height, units = 'in')
}

save_example_figure <- function() {
  save_figure(example_tpca_figure(), 
              'ex_hellinger-d20-k010-3')
}
