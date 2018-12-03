plot.tpca <- function(tpca_obj, ...) {
  par(ask = TRUE)
  titles <- list('Estimated probability of axis n being the most sensitive',
                 'Pointwise 0.25 and 0.975 quantiles of the divergences',
                 'Estimated divergence for different change types',
                 'Estimated divergence for different change sparsities')
  plot_funcs <- list(ggplot_prop_max, 
                     ggplot_quantiles, 
                     ggplot_types_mean, 
                     ggplot_sparsity_mean)
  capture.output(Map(function(f, title) f(tpca_obj, title = title),
                     plot_funcs, titles))
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

ggplot_sparsity_mean <- function(tpca_obj, 
                                 sparsities = unique(tpca_obj$change_sparsity),
                                 title      = NULL,
                                 xlab       = 'Principal axis nr.',
                                 ylab       = 'Hellinger distance') {
  max.sparsities <- 100
  sparsities = sort(unique(tpca_obj$change_sparsity))
  if (length(sparsities) > max.sparsities) {
    ind <- round(seq(1, length(sparsities), length.out = max.sparsities))
    sparsities <- sparsities[ind]
  }
  
  data_dim <- nrow(tpca_obj$divergence_sim)
  sparsity_subsets <- lapply(sparsities, function(sparsity) {
    subset_sims(tpca_obj, sparsity = sparsity)
  })
  sparsity_means <- lapply(sparsity_subsets, rowMeans)
  sparsity_means[[length(sparsities) + 1]] <-  rowMeans(tpca_obj$divergence_sim)
  sparsity_means <- do.call('cbind', sparsity_means)
  
  colnames(sparsity_means) <- c(sparsities, max(sparsities) + 1)
  means_df <- melt(sparsity_means, 
                   varnames   = c('axis', 'K'), 
                   value.name = 'divergence')
  means_df$K <- as.factor(means_df$K)
  
  label_ind <- seq.int(1, length(sparsities), length.out = min(length(sparsities), 5))
  legend_labels <- vapply(sparsities[label_ind], 
                          function(x) paste('K =', x), 
                          character(1))
  legend_labels <- c(legend_labels, 'Overall')
  col <- c(colorRampPalette(c('orange', 'blue'))(length(sparsities)), 'black')
  line_sizes <- c(rep(0.3, length(col) - 1), 0.6)
  
  ggplot(means_df, aes(x = axis, y = divergence, color = K)) +
    geom_line(aes(size = K)) +
    theme_light() +
    labs(x = xlab, y = ylab) +
    ggtitle(title) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual('Change \n sparsity', values = col,
                       breaks = c(sparsities[label_ind], max(sparsities) + 1),
                       labels = legend_labels) +
    scale_size_manual(values = line_sizes, guide = FALSE)
}

ggplot_quantiles <- function(tpca_obj, 
                             quantiles = c(0.25, 0.975),
                             title     = NULL,
                             xlab      = 'Principal axis nr.',
                             ylab      = 'Hellinger distance') {
  divergence_sims <- tpca_obj$divergence_sim
  data_dim <- nrow(divergence_sims)
  mean_divergence <- rowMeans(divergence_sims)
  divergence_quantiles <- apply(divergence_sims, 1, quantile, probs = quantiles)
  plot_df <- data.frame('axis' = 1:data_dim,
                        'mean' = mean_divergence)
  plot_df <- cbind(plot_df, t(divergence_quantiles))
  plot_df <- melt(plot_df, id.vars = 'axis')
  ggplot(plot_df, aes(x = axis, y = value, linetype = variable)) +
    geom_line() +
    theme_light() +
    labs(x = xlab, y = ylab) +
    ggtitle(title) +
    scale_linetype_manual(guide = FALSE, 
                          values = c(1, rep(2, nrow(divergence_quantiles))))
}
