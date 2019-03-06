# Other runs:
#   * Larger dimension, d.
#   * Other change distributions (change sizes, sparsity (already done by sparsity plot)).
#      - Smaller sizes, larger sizes.
#   * Show variance, quantile plots.
#   * All affected streams changes by the same size.

library(doParallel)

avg_hellinger_sim2 <- function(d, n_cov_mat, n_sim, change_distr) {
  init_log_file <- function() {
    write('+++++++++++++++++++++++++++++++++++++++++++++', 
          file = log_file, append = TRUE)
    write(paste0('d = ', d, ', change_distr = ', change_distr, 
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
    cor_mat <- tpca::rcor_mat(d)
    tpca::tpca(cor_mat, change_distr = change_distr, n_sim = n_sim)
  }
  stopCluster(c1)
  R_file <- paste0(dir, 'tpca_list_', change_distr, '_D', d, '.RData')
  save(tpca_list, file = R_file)
  invisible(tpca_list)
}

run_Hsim_other_cd <- function() {
  run_sim <- function(change_distr) {
    d <- 100
    n_cov_mat <- 10^3
    n_sim <- 10^3
    avg_hellinger_sim2(d, n_cov_mat, n_sim, change_distr)
  }
  
  change_distr <- c('full_uniform_equal', 'full_uniform_large', 'full_uniform_small')
  lapply(change_distr, run_sim)
}

run_Hsim_uniform_D200 <- function() {
  d <- 200
  n_cov_mat <- 10^3
  n_sim <- 10^3
  change_distr <- 'full_uniform'
  avg_hellinger_sim2(d, n_cov_mat, n_sim, change_distr)
}

run_all_additional_sims <- function() {
  run_Hsim_other_cd()
  run_Hsim_uniform_D200()
}

get_avg_figure <- function(change_distr, title = FALSE, show = FALSE) {
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
