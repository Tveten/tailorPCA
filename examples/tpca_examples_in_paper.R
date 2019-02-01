N <- 20
k <- N / 2
set.seed(7)
Sigma1.x <- generate_cor_mat(N, K0 = N/2)

tpca_obj <- tpca()

# Entirely fixed change (also affected dims)
fixed_k <- function(n) {
  rep(N - 50, n)
}
fixed_mu <- function(n) {
  rep(1, n)
}
fixed_sigma <- function(n) {
  rep(2, n)
}
fixed_rho <- function(n) {
  rep(0, n)
}
fixed_dims <- function(k) {
  set.seed(15)
  sample(1:N, k)
}

est.hellinger.obj <- est_hellinger(Sigma1.x, 
                                   PCA.obj, 
                                   n.sim = 10^2, 
                                   return.all = TRUE,
                                   k_distr = fixed_k,
                                   mu_distr = fixed_mu,
                                   sigma_distr = fixed_sigma,
                                   rho_distr = fixed_rho,
                                   affected_dims = fixed_dims)
type.plot.ex <- ggplot_est_type_avg(list(est.hellinger.obj))
type.plot.ex
ggsave('Example100HellingerTypeFixed8.pdf', type.plot.ex,
       width = 4, height = 2.6)