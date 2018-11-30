plot.tpca <- function(tpca_obj, ...) {
  par(ask = TRUE)
  for (i in 1:n_plots) {
    
  }
}

subset_sims <- function(tpca_obj, type = unique(sim_type), 
                        sparsity = unique(sim_sparsity)) {
  sim <- tpca_obj$divergence_sim
  sim_type <- tpca_obj$change_type
  sim_sparsity <- tpca_obj$change_sparsity
  sim[, (sim_type %in% type) & (sim_sparsity %in% sparsity)]
}

ggplot_prop_max <- function(tpca_obj, 
                            type = unique(tpca_obj$change_type),
                            sparsity = unique(tpca_obj$change_sparsity),
                            title = NULL,
                            xlab = 'Principal axis nr. n',
                            ylab = 'Proportion n is most sensitive') {
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

ggplot_types_mean <- function(tpca_obj) {
  full_mean <- rowMeans(tpca_obj$divergence_sim)
  
  types <- sort(unique(tpca_obj$change_type))
  type_subsets <- lapply(types, function(type) subset_sims(tpca_obj, type = type))
  type_means <- lapply(sim_subsets, rowMeans)
  if (any(types == 'mean')) {
    ind <- which(types == 'mean')
    name()
  }
}

ggplot_types_avg <- function(tpca_obj) {
  n.cor.mat <- length(h.obj.list)
  types <- sort(unique(h.obj.list[[1]]$type))
  N <- length(h.obj.list[[1]]$est)
  h.est.array <- array(0, dim = c(N, n.cor.mat, length(types)))
  for (i in 1:n.cor.mat) {
    for (j in 1:length(types)) {
      h.est.array[, i, j] <- get_type_est(h.obj.list[[i]], types[j])
    }
  }
  h.est.type <- apply(h.est.array, c(1, 3), mean)
  h.est.avg <- rowMeans(h.est.type)
  
  h.est.df <- data.frame(PA.nr = integer(), 
                         h.est = double(), 
                         type = character())
  cols <- c(rainbow(3), 'black')
  col <- rep(NA, length(types) + 1)
  if (any(types == 'mu')) {
    ind <- which(types == 'mu')
    col[ind] <- cols[1]
    mean.part <- data.frame(PA.nr = 1:nrow(h.est.type),
                            h.est = h.est.type[, ind], 
                            type  = rep('Mean', nrow(h.est.type)))
    h.est.df <- rbind(h.est.df, mean.part)
  }
  if (any(types == 'sigma')) {
    ind <- which(types == 'sigma')
    col[ind] <- cols[2]
    var.part <- data.frame(PA.nr = 1:nrow(h.est.type),
                           h.est = h.est.type[, ind], 
                           type  = rep('Variance', nrow(h.est.type)))
    h.est.df <- rbind(h.est.df, var.part)
  }
  if (any(types == 'rho')) {
    ind <- which(types == 'rho')
    col[ind] <- cols[3]
    cor.part <- data.frame(PA.nr = 1:nrow(h.est.type),
                           h.est = h.est.type[, ind], 
                           type  = rep('Correlation', nrow(h.est.type)))
    h.est.df <- rbind(h.est.df, cor.part)
  }
  
  avg.part <- data.frame(PA.nr = 1:nrow(h.est.type),
                         h.est = h.est.avg,
                         type = rep('Average', nrow(h.est.type)))
  h.est.df <- rbind(h.est.df, avg.part)
  col[length(col)] <- cols[length(cols)]
  
  ggplot(h.est.df, aes(x = PA.nr, y = h.est, color = type)) +
    geom_line(aes(size = type)) +
    theme_light() +
    labs(x = 'Principal axis nr.', y = 'Hellinger distance') +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual('Change \n type', values = col) +
    scale_size_manual(values = c(rep(0.3, length(col) - 1), 0.6),
                      guide = FALSE)
}

