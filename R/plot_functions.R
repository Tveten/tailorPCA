plot.tpca <- function(tpca_obj, ...) {
  par(ask = TRUE)
  plot_funcs <- list(ggplot_prop_max, 
                     ggplot_ci, 
                     ggplot_types_mean, 
                     ggplot_sparsity_mean)
  capture.output(lapply(plot_funcs, function(f) f(tpca_obj)))
  invisible(NULL)
}

ggplot_prop_max <- function(tpca_obj, 
                            type     = unique(tpca_obj$change_type),
                            sparsity = unique(tpca_obj$change_sparsity),
                            title    = NULL,
                            xlab     = 'Principal axis nr. n',
                            ylab     = 'Proportion n is most sensitive') {
  sims <- subset_sims(tpca_obj, type = type, sparsity = sparsity)
  prop_axes_max <- prop_axes_max(sims)
  data_dim <- length(prop_axes_max)
  prop_max_df <- data.frame('axes' = 1:data_dim,
                            'prop_max' = prop_axes_max)
  ggplot(prop_max_df, aes(x = axes, y = prop_max)) +
    geom_bar(stat = 'identity') +
    theme_light() +
    ggtitle(title) +
    labs(x = xlab, y = ylab) +
    scale_y_continuous(limits = c(0, 1))
}

ggplot_types_mean <- function(tpca_obj, 
                              types = unique(tpca_obj$change_type),
                              title = NULL,
                              xlab  = 'Principal axis nr.',
                              ylab  = 'Hellinger distance') {
  data_dim <- nrow(tpca_obj$divergence_sim)
  type_subsets <- lapply(types, function(type) subset_sims(tpca_obj, type = type))
  type_means <- lapply(type_subsets, rowMeans)
  
  shown_names <- list('mean' = 'Mean', 'sd' = 'Variance', 'cor' = 'Correlation')
  names(type_means) <- unlist(shown_names[types])
  type_means$Overall <- rowMeans(tpca_obj$divergence_sim)
  
  means_df <- melt(type_means, value.name = 'divergence')
  names(means_df)[2] <- 'type'
  ordered_levels <- c('Mean', 'Variance', 'Correlation', 'Overall')
  means_df$type <- factor(means_df$type, levels = ordered_levels)
  means_df$axis <- rep(1:data_dim, length(types) + 1)
  
  line_sizes <- c(0.3, 0.3, 0.3, 0.6)
  line_cols <- c('#FF0000FF', '#0000FFFF', '#00FF00FF', 'black')
  names(line_cols) <- ordered_levels
  names(line_cols) <- ordered_levels
  
  ggplot(means_df, aes(x = axis, y = divergence, color = type)) +
    geom_line(aes(size = type)) +
    theme_light() +
    labs(x = xlab, y = ylab) +
    ggtitle(title) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual('Change \n type', values = line_cols) +
    scale_size_manual(values = line_sizes, guide = FALSE)
}

