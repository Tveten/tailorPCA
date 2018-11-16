get_type_est <- function(tpca_obj, type) {
  if (type == 'avg') return(h.obj$est)
  else {
    h.sim.type <- h.obj$sims[, h.obj$type == type, drop = FALSE]
    rowMeans(h.sim.type, na.rm = TRUE)
  }
}

get_sparsity_est <- function(tpca_obj, k) {
  h.sim.sparsity <- h.obj$sims[, h.obj$sparsity == k, drop = FALSE]
  sparsity.est <- rowMeans(h.sim.sparsity, na.rm = TRUE)
  # if (any(is.nan(rowMeans(h.sim.sparsity, na.rm = TRUE)))) 
  #   print(h.sim.sparsity)
  if (any(is.nan(sparsity.est))) sparsity.est <- rep_len(NA, length(sparsity.est))
  sparsity.est
}
