avg_hellinger_sim <- function(d, n_cov_mat, n_sim, change_distr) {
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
  
  dir.create('examples', showWarnings = FALSE)
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
    avg_hellinger_sim(d, n_cov_mat, n_sim, change_distr)
  }
  
  change_distr <- c('full_uniform', 'full_uniform_equal', 
                    'full_uniform_large', 'full_uniform_small')
  lapply(change_distr, run_sim)
}

run_Hsim_uniform_D200 <- function() {
  d <- 200
  n_cov_mat <- 10^3
  n_sim <- 10^3
  change_distr <- 'full_uniform'
  avg_hellinger_sim(d, n_cov_mat, n_sim, change_distr)
}

run_all_additional_sims <- function() {
  run_Hsim_other_cd()
  run_Hsim_uniform_D200()
}

get_avg_figure <- function(tpca_obj, title = FALSE, show = FALSE) {
  get_title <- function(file_id) {
    file_id
  }
  
  if (title) type_plot <- ggplot_types(tpca_obj, title = get_title())
  else type_plot <- ggplot_types(tpca_obj)
  sparsity_plot <- ggplot_sparsities(tpca_obj)
  
  if (show) show(gridExtra::grid.arrange(type_plot, sparsity_plot, nrow = 1))
  list('type' = type_plot, 'sparsity' = sparsity_plot)
}

get_prop_figure <- function(tpca_obj, show = FALSE) {
  prop_plot <- ggplot_prop(tpca_obj)
  
  if (show) show(prop_plot)
  prop_plot
}

tpca_summary_figure <- function(file_id, show = FALSE) {
  # This function assumes that avg_hellinger_sim have been run, such that
  # an .RData-object has been saved in the examples directory.
  
  dir <- './examples/'
  # Always loads object named tpca_list.
  load(paste0(dir, 'tpca_list_', file_id, '.RData'))
  tpca_obj <- merge_tpca_list(tpca_list)
  avg_figures <- get_avg_figure(tpca_obj)
  prop_figure <- get_prop_figure(tpca_obj)
  
  if (show) show(gridExtra::grid.arrange(avg_figures$type, avg_figures$sparsity, prop_figure, nrow = 1))
  invisible(gridExtra::arrangeGrob(avg_figures$type, avg_figures$sparsity, prop_figure, nrow = 1))
}

save_summary_figure <- function(file_id) {
  figure <- tpca_summary_figure(file_id)
  file_name <- paste0('hellinger_summary_', file_id)
  save_figure(figure, file_name)
}

save_all_additional_figures <- function() {
 file_ids <- c('full_uniform_equal_D100', 'full_uniform_large_D100', 
               'full_uniform_small_D100', 'full_uniform_D200') 
  invisible(lapply(file_ids, save_summary_figure))
}