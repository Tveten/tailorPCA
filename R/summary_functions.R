summary.tpca <- function(tpca.obj) {
  # For future use.
  if (any(change_type == 'mean')) {
    h_sim_mean <- hellinger_sim[, change_type == 'mean']
    prop_max_mean <- get_prop(h_sim_mean)
  }
  if (any(change_type == 'sd')) {
    h_sim_sd <- hellinger_sim[, change_type == 'sd']
    prop_max_sd <- get_prop(h_sim_sd)
  }
  if (any(change_type == 'cor')) {
    h_sim_mean <- hellinger_sim[, change_type == 'mean']
    prop_max_mean <- get_prop(h_sim_meanmu)
  }
  prop_max_mean <- get_prop(h_sim_meanmu)
  prop.top.sigma <- get_prop(h.sim.sigma)
  prop.top.rho <- get_prop(h.sim.rho)
}
    