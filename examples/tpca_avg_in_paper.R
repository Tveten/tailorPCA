library(doParallel)

avg_hellinger_sim <- function(d, k0, n_cov_mat, n_sim, cov_mat_type,
                              change_distr = 'full_uniform') {
  init_log_file <- function() {
    write('+++++++++++++++++++++++++++++++++++++++++++++', 
          file = log_file, append = TRUE)
    write(paste0('d = ', d, ', k0 = ', k0, 
                 ', n_cov_mats = ', n_cov_mat, ', n_sim = ', n_sim, '.'),
          file = log_file, append = TRUE)
    write('+++++++++++++++++++++++++++++++++++++++++++++', 
          file = log_file, append = TRUE)
  }
  log_current <- function() {
    time_used <- proc.time()[3]/60
    write(paste0('Sim nr. ', b, '. Time used: ', time_used, ' min.'),
          file = log_file, append = TRUE)
  }
  
  dir <- './examples/'
  log_file <- paste0(dir, 'avg_hellinger_log_file.txt')
  init_log_file()
  
  c1 <- makeCluster(4, outfile = '', type = 'PSOCK')
  registerDoParallel(c1)
  tpca_list <- foreach(b = 1:n_cov_mat) %dopar% {
    log_current()
    cor_mat <- tpca::rcor_mat(d, k0)
    tpca::tpca(cor_mat, change_distr = change_distr, n_sim = n_sim)
  }
  stopCluster(c1)
  R_file <- paste0(dir, 'tpca_list_', cov_mat_type, '.RData')
  save(tpca_list, file = R_file)
  invisible(tpca_list)
}

test_hellinger_sim <- function(d, k0) {
  n_cov_mat <- 10
  n_sim <- 100
  avg_hellinger_sim(d, k0, n_cov_mat, n_sim, 'test')
}

run_avg_hellinger_sim <- function() {
  run_sim <- function(k0, cov_mat_type) {
    d <- 100
    n_cov_mat <- 10^3
    n_sim <- 10^3
    avg_hellinger_sim(d, k0, n_cov_mat, n_sim, cov_mat_type)
  }
  
  k0s <- c(10, 50, 100)
  cov_mat_types <- c('sparse', 'halfsparse', 'dense')
  Map(run_sim, k0s, cov_mat_types)
}

get_avg_figure <- function(cov_mat_type, title = FALSE, show = FALSE) {
  get_title <- function(cov_mat_type) {
    if (cov_mat_type == 'halfsparse') return('Half-sparse')
    else return(first_up(cov_mat_type))
  }
  
  dir <- './examples/'
  # Always returns object named tpca_list.
  load(paste0(dir, 'tpca_list_', cov_mat_type, '.RData'))
  tpca_obj <- merge_tpca_list(tpca_list)
  if (title) type_plot <- ggplot_types_mean(tpca_obj, title = title_str)
  else type_plot <- ggplot_types_mean(tpca_obj)
  sparsity_plot <- ggplot_sparsity_mean(tpca_obj)
  
  if (show) gridExtra::grid.arrange(type_plot, sparsity_plot, nrow = 2)
  list('type' = type_plot, 'sparsity' = sparsity_plot)
}

avg_tpca_figure <- function() {
  cov_mat_types <- c('dense', 'halfsparse', 'sparse')
  figures <- lapply(cov_mat_types, get_avg_figure)
  names(figures) <- cov_mat_types
  gridExtra::arrangeGrob(figures$dense$type, figures$halfsparse$type,
                         figures$sparse$type, figures$dense$sparsity, 
                         figures$halfsparse$sparsity, figures$sparse$sparsity,
                         nrow = 2)
}

avg_tpca_figure_dense <- function() {
  figures <- get_avg_figure('dense')
  gridExtra::arrangeGrob(figures$type, figures$sparsity, nrow = 1)
}

save_avg_figure_dense <- function() {
  figure <- avg_tpca_figure_dense()
  name <- 'avg_hellinger_dense-d100'
  save_figure(figure, name)
}

save_avg_figure <- function() {
  figure <- avg_tpca_figure()
  name <- 'avg_hellinger-d100'
  save_figure(figure, name)
}
