subset_sims <- function(tpca_obj, type = unique(sim_type), 
                        sparsity = unique(sim_sparsity)) {
  sim <- tpca_obj$divergence_sim
  sim_type <- tpca_obj$change_type
  sim_sparsity <- tpca_obj$change_sparsity
  sim[, (sim_type %in% type) & (sim_sparsity %in% sparsity)]
}

save_figure <- function(ggplot_obj, name, 
                        base_width  = NULL, 
                        base_height = NULL, 
                        extension = 'png') {
  plot_grid <- dim(ggplot_obj)
  A4_width <- 8.27  # inches
  A4_height <- 11.69  # inches
  if (is.null(base_width)) base_width <- 2.88
  if (is.null(base_height)) base_height <- 2.3
  if (is.null(plot_grid)) {
    width = base_width
    height = base_height
  } else {
    width = plot_grid[2] * base_width
    height = plot_grid[1] * base_height
  }
  ggplot2::ggsave(paste0(name, '.', extension), ggplot_obj,
                  width = width, height = height, units = 'in')
}

merge_tpca_list <- function(tpca_list) {
  div_sim_list <- list()
  type_list <- list()
  sparsity_list <- list()
  for (b in seq_along(tpca_list)) {
    div_sim_list[[b]] <- tpca_list[[b]]$divergence_sim
    type_list[[b]] <- tpca_list[[b]]$change_type
    sparsity_list[[b]] <- tpca_list[[b]]$change_sparsity
  }
  tpca_obj <- list()
  tpca_obj$divergence_sim <- do.call('cbind', div_sim_list)
  tpca_obj$change_type <- do.call('c', type_list)
  tpca_obj$change_sparsity <- do.call('c', sparsity_list)
  tpca_obj
}

first_up <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get_ylab <- function(divergence) {
  if (divergence == 'normal_hellinger') ylab <- 'Hellinger dist.'
  else if (divergence == 'normal_KL') ylab <- 'KL divergence'
  else if (divergence == 'normal_bhat') ylab <- 'Bhattacharyya dist.'
  else ylab <- 'Divergence'
  ylab
}