#' @export
plot.tpca <- function(tpca_obj, ...) {
  par(ask = TRUE)
  on.exit(par(ask = FALSE))
  titles <- list('Estimated probability of axis n being the most sensitive',
                 'Pointwise 0.25 and 0.975 quantiles of the divergences',
                 'Estimated divergence for different change types',
                 'Estimated divergence for different change sparsities')
  plot_funcs <- list(ggplot_prop, 
                     ggplot_quantiles, 
                     ggplot_types, 
                     ggplot_sparsities)
  capture.output(Map(function(f, title) f(tpca_obj, title = title),
                     plot_funcs, titles))
  invisible(NULL)
}

#' @export
ggplot_prop <- function(tpca_obj, 
                        type     = unique(tpca_obj$change_type),
                        sparsity = unique(tpca_obj$change_sparsity),
                        title    = NULL,
                        print_p  = FALSE) {
  sims <- subset_sims(tpca_obj, type = type, sparsity = sparsity)
  prop_axes_max <- prop_axes_max(sims)
  data_dim <- length(prop_axes_max)
  prop_max_df <- data.frame('axes' = 1:data_dim,
                            'prop_max' = prop_axes_max)
  if (print_p) {
    print(title)
    print(prop_max_df)
  }
  
  xlab <- 'Projection j'
  ylab <- 'Prob. j is most sensitive'
  
  ggplot2::ggplot(prop_max_df, ggplot2::aes(x = axes, y = prop_max)) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::theme_light() +
    ggplot2::ggtitle(title) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::scale_y_continuous(limits = c(0, 1))
}

#' @export
ggplot_types <- function(tpca_obj, 
                         types = unique(tpca_obj$change_type),
                         title = NULL) {
  data_dim <- nrow(tpca_obj$divergence_sim)
  type_subsets <- lapply(types, function(type) subset_sims(tpca_obj, type = type))
  type_means <- lapply(type_subsets, rowMeans)
  
  shown_names <- list('mean' = 'Mean', 'sd' = 'Var', 'cor' = 'Cor')
  names(type_means) <- unlist(shown_names[types])
  type_means$All <- rowMeans(tpca_obj$divergence_sim)
  
  means_df <- reshape2::melt(type_means, value.name = 'divergence')
  names(means_df)[2] <- 'type'
  ordered_levels <- c('Mean', 'Var', 'Cor', 'All')
  means_df$type <- factor(means_df$type, levels = ordered_levels)
  means_df$axis <- rep(1:data_dim, length(types) + 1)
  
  line_sizes <- c(0.3, 0.3, 0.3, 0.6)
  line_cols <- c('#FF0000FF', '#0000FFFF', '#00FF00FF', 'black')
  names(line_sizes) <- ordered_levels
  names(line_cols) <- ordered_levels
  xlab <- 'Projection j'
  ylab <- paste0('E[ ', get_ylab(tpca_obj$divergence), ' ]')
  
  ggplot2::ggplot(means_df, ggplot2::aes(x = axis, y = divergence, color = type)) +
    ggplot2::geom_line(ggplot2::aes(size = type)) +
    ggplot2::theme_light() +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual('Change \n type', values = line_cols) +
    ggplot2::scale_size_manual(values = line_sizes, guide = FALSE)
}

#' @export
ggplot_sparsities <- function(tpca_obj, 
                              sparsities = unique(tpca_obj$change_sparsity),
                              title      = NULL) {
  set_color <- function(sparsities) {
    n_legends <- min(length(sparsities), 5)
    label_ind <- floor(seq(1, length(sparsities), length.out = n_legends))
    legend_labels <- vapply(sparsities[label_ind], 
                            function(x) paste('K =', x), 
                            character(1))
    legend_labels <- c(legend_labels, 'All')
    col <- c(colorRampPalette(c('orange', 'blue'))(length(sparsities)), 'black')
    legend_breaks <- c(sparsities[label_ind], max(sparsities) + 1)
    line_sizes <- c(rep(0.3, length(col) - 1), 0.6)
    list('col'           = col, 
         'legend_breaks' = legend_breaks, 
         'legend_labels' = legend_labels)
  }
  
  set_shown_sparsities <- function() {
    max_sparsities <- 100
    sparsities <- sort(unique(tpca_obj$change_sparsity))
    if (length(sparsities) > max_sparsities) {
      ind <- round(seq(1, length(sparsities), length.out = max_sparsities))
      sparsities <- sparsities[ind]
    }
    sparsities
  }
  
  sparsities <- set_shown_sparsities()
  
  data_dim <- nrow(tpca_obj$divergence_sim)
  sparsity_subsets <- lapply(sparsities, function(sparsity) {
    subset_sims(tpca_obj, sparsity = sparsity)
  })
  sparsity_means <- lapply(sparsity_subsets, rowMeans)
  sparsity_means[[length(sparsities) + 1]] <-  rowMeans(tpca_obj$divergence_sim)
  sparsity_means <- do.call('cbind', sparsity_means)
  colnames(sparsity_means) <- c(sparsities, max(sparsities) + 1)
  
  means_df <- reshape2::melt(sparsity_means, 
                             varnames   = c('axis', 'K'), 
                             value.name = 'divergence')
  means_df$K <- as.factor(means_df$K)
  
  col_obj <- set_color(sparsities)
  line_sizes <- c(rep(0.3, length(col_obj$col) - 1), 0.6)
  xlab <- 'Projection j'
  ylab <- paste0('E[ ', get_ylab(tpca_obj$divergence), ' ]')
  
  ggplot2::ggplot(means_df, ggplot2::aes(x = axis, y = divergence, color = K)) +
    ggplot2::geom_line(ggplot2::aes(size = K)) +
    ggplot2::theme_light() +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_color_manual('Change \n sparsity', 
                                values = col_obj$col,
                                breaks = col_obj$legend_breaks, 
                                labels = col_obj$legend_labels) +
    ggplot2::scale_size_manual(values = line_sizes, guide = FALSE)
}

#' @export
ggplot_quantiles <- function(tpca_obj, 
                             quantiles = c(0.025, 0.975),
                             title     = NULL) {
  divergence_sims <- tpca_obj$divergence_sim
  data_dim <- nrow(divergence_sims)
  mean_divergence <- rowMeans(divergence_sims)
  divergence_quantiles <- apply(divergence_sims, 1, quantile, probs = quantiles)
  plot_df <- data.frame('axis' = 1:data_dim,
                        'mean' = mean_divergence)
  plot_df <- cbind(plot_df, t(divergence_quantiles))
  plot_df <- reshape2::melt(plot_df, id.vars = 'axis')
  
  xlab <- 'Projection j'
  ylab <- get_ylab(tpca_obj$divergence)
  
  ggplot2::ggplot(plot_df, ggplot2::aes(x = axis, y = value, linetype = variable)) +
    ggplot2::geom_line() +
    ggplot2::theme_light() +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_linetype_manual(guide = FALSE, 
                          values = c(1, rep(2, nrow(divergence_quantiles))))
}

#' @export
ggplot_singles <- function(tpca_obj, 
                           n = NULL,
                           title = NULL) {
  data_dim <- nrow(tpca_obj$divergence_sim)
  n_sim <- ncol(tpca_obj$divergence_sim)
  types = unique(tpca_obj$change_type)
  if (is.null(n)) {
    all_type_inds <- lapply(types, function(type) which(tpca_obj$change_type == type))
    random_type_ind <- vapply(all_type_inds, function(ind) sample(ind, 1), numeric(1))
    change_sparsities <- tpca_obj$change_sparsity[random_type_ind]
    # names(type_sample) <- types
    divergence_subset <- tpca_obj$divergence_sim[, random_type_ind]
    shown_names <- list('mean' = paste0('Mean, \n K = ', change_sparsities[types == 'mean']),
                        'sd' = paste0('Var, \n K = ', change_sparsities[types == 'sd']), 
                        'cor' = paste0('Cor, \n K = ', change_sparsities[types == 'cor']))
    colnames(divergence_subset) <- unlist(shown_names[types])
    divergence_df <- reshape2::melt(divergence_subset)
    names(divergence_df) <- c('axis', 'type', 'divergence')
    ordered_levels <- shown_names[c('mean', 'sd', 'cor')]
    divergence_df$type <- factor(divergence_df$type, levels = ordered_levels)
    line_sizes <- c(0.3, 0.3, 0.3)
    line_cols <- c('#FF0000FF', '#0000FFFF', '#00FF00FF')
    names(line_sizes) <- ordered_levels
    names(line_cols) <- ordered_levels
    xlab <- 'Projection j'
    ylab <- get_ylab(tpca_obj$divergence)
    
    ggplot2::ggplot(divergence_df, ggplot2::aes(x = axis, y = divergence, color = type)) +
      ggplot2::geom_line(ggplot2::aes(size = type)) +
      ggplot2::theme_light() +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::ggtitle(title) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_color_manual('Change \n type', values = line_cols) +
      ggplot2::scale_size_manual(values = line_sizes, guide = FALSE)
  } else {
    random_ind <- sample(1:n_sim, n)
    change_sparsities <- tpca_obj$change_sparsity[random_ind]
    divergence_subset <- tpca_obj$divergence_sim[, random_ind]
    shown_names <- lapply(1:n, function(i) {
      ind <- random_ind[i]
      paste0(i, ':  ', tpca_obj$change_type[ind], ', K = ', tpca_obj$change_sparsity[ind])
    })
    colnames(divergence_subset) <- unlist(shown_names)
    divergence_df <- reshape2::melt(divergence_subset)
    names(divergence_df) <- c('axis', 'type', 'divergence')
    
    ggplot2::ggplot(divergence_df, ggplot2::aes(x = axis, y = divergence, color = type)) +
      ggplot2::geom_line() +
      ggplot2::theme_light() +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::ggtitle(title) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::guides(color = ggplot2::guide_legend(title = 'Change \n type, sparsity'))
  }
}