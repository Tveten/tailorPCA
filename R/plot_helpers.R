subset_sims <- function(tpca_obj, type = unique(sim_type), 
                        sparsity = unique(sim_sparsity)) {
  sim <- tpca_obj$divergence_sim
  sim_type <- tpca_obj$change_type
  sim_sparsity <- tpca_obj$change_sparsity
  sim[, (sim_type %in% type) & (sim_sparsity %in% sparsity)]
}

