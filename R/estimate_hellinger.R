# setwd('C:/Users/matve/OneDrive/Dokumenter/Studier/PhD/DimRedPaper/R')
# source('GenerateData.R')
# source('DimRedFunctions.R')
# setwd('C:/Users/matve/OneDrive/Dokumenter/Studier/PhD/DimRedPaper/R/Output')
library(todor)

# TODO: 

pca <- function(Sigma, p_axes = 1:N, eigen.values = FALSE) {
  N <- nrow(Sigma)
  if (eigen.values) {
    eigen.obj <- svd(Sigma, nv = 0)
    A <- t(eigen.obj$u[, i.PC])
    lambda <- eigen.obj$d[i.PC]
    return(list('vectors' = A, 'values' = lambda))
  } else {
    eigen.obj <- svd(Sigma, nv = 0)
    A <- t(eigen.obj$u[, i.PC])
    return(A)
  }
}

calc_hellinger <- function(mu1, sigma1, mu2, sigma2) {
  # sigma_i: standard deviation
  sqrt(1 - sqrt(2 * sigma1 * sigma2 / (sigma1^2 + sigma2^2)) *
         exp(-1/4 * (mu1 - mu2)^2 / (sigma1^2 + sigma2^2)))
}


change_cor_mat <- function(R, affected.dims, 
                           draw_rho = NULL, draw_sigma = NULL) {
  # At least one of the NULL-arguments must be supplied:
  #   functions draw_rho or draw_sigma
  #
  # Returns:
  #   Sigma2: The change covariance matrix.
  
  change_cor <- function(R.sub, draw_rho, n.affected) {
    affected.cor.dims <- sample(1:K0, min(n.affected, K0))
    ind <- combn(affected.cor.dims, 2)
    change.factor <- draw_rho(sum(1:length(ind)))
    R.changed <- R.sub
    for (i in 1:ncol(ind)) {
      R.changed[ind[1, i], ind[2, i]] <- change.factor[i] * R.changed[ind[1, i], ind[2, i]]
      R.changed[ind[2, i], ind[1, i]] <- change.factor[i] * R.changed[ind[2, i], ind[1, i]]
    }
    R.changed <- as.matrix(nearPD(R.changed, 
                                  corr = TRUE, 
                                  maxit = 20,
                                  do2eigen = TRUE,
                                  posd.tol = 1e-8)$mat)
  }
  
  N <- ncol(R)
  K0 <- sum(R[, 1] != 0)
  K <- length(affected.dims)
  Sigma2 <- R
  
  if (all(is.null(c(draw_rho, draw_sigma))))
    stop('ERROR: Either a variance or a correlation change distribution must be specified')
  
  # Correlation change handling
  if (!is.null(draw_rho)) {
    R.sub <- R[1:K0, 1:K0, drop = FALSE]
    Sigma2[1:K0, 1:K0] <- change_cor(R.sub, draw_rho, K)
  }
  
  # Variance change handling
  if (!is.null(draw_sigma)) {
    sigma2.vec <- rep(1, N)
    sigma2.vec[affected.dims] <- draw_sigma(K)
    Sigma2 <- diag(sigma2.vec) %*% Sigma2 %*% diag(sigma2.vec)
  }
  
  return(Sigma2)
}


est_hellinger <- function(Sigma1.x, PCA.obj, n.sim = 10^3, 
                          return.all = FALSE,  ...) {
  # Optional arguments:
  # p: Probabilities of (mean change, variance change, correlation change)
  # k.distr: Change distribution of the sparsity.
  # mu.distr: Change distribution of the change in mean.
  # sigma.distr: Change distribution of the change in variance.
  # rho.distr: Change distribution of the change in correlation.
  
  # Setup
  N <- ncol(Sigma1.x)
  U <- PCA.obj$vectors
  
  # Default change distributions
  p <- rep(1/3, 3)  # Probability of each type of change.
  draw_k <- function(n) {
    sample(2:(N/2), n, replace = TRUE)
  }
  draw_mu <- function(n) {
    runif(n, -1.5, 1.5)
  }
  draw_sigma <- function(n) {
    change.factors <- runif(n, 1, 2.5)
    ind.decrease <- sample(c(TRUE, FALSE), n, replace = TRUE)
    change.factors[ind.decrease] <- sapply(change.factors[ind.decrease],
                                           function(x) 1/x)
    change.factors
  }
  draw_rho <- function(n) {
    runif(n, 0, 1)
  }
  affected_dims <- function(k) {
    sample(1:N, k)
  }
  
  # Optional arguments for overriding the default distributions.
  args <- list(...)
  if (length(args) > 0) {
    items <- names(args)
    if (any(items == 'p')) p <- args$p
    if (any(items == 'k_distr')) draw_k <- args$k_distr
    if (any(items == 'mu_distr')) draw_mu <- args$mu_distr
    if (any(items == 'sigma_distr')) draw_sigma <- args$sigma_distr
    if (any(items == 'rho_distr')) draw_rho <- args$rho_distr
    if (any(items == 'affected_dims')) affected_dims <- args$affected_dims
  }
  
  # Initialization
  change.types <- c('mu', 'sigma', 'rho')
  change.type <- sample(change.types, n.sim, prob = p, replace = TRUE)
  k <- draw_k(n.sim)
  
  mu1 <- rep(0, N)
  sigma1 <- sqrt(PCA.obj$values)
  hellinger.sim <- matrix(NA, nrow = N, ncol = n.sim)
  
  for (b in 1:n.sim) {
    affected.dims <- affected_dims(k[b])
    if (change.type[b] == 'mu') {
      mu2.x <- rep(0, N)
      mu2.x[affected.dims] <- draw_mu(k[b])
      mu2 <- U %*% mu2.x
      hellinger.sim[, b] <- round(calc_hellinger(mu1, sigma1, mu2, sigma1), 6)
    } 
    if (change.type[b] == 'sigma') {
      Sigma2.x <- change_cor_mat(Sigma1.x, affected.dims, draw_sigma = draw_sigma)
      sigma2 <- sqrt(apply(U, 1, function(u) u %*% Sigma2.x %*% u))
      hellinger.sim[, b] <- round(calc_hellinger(mu1, sigma1, mu1, sigma2), 6)
    } 
    if (change.type[b] == 'rho') {
      Sigma2.x <- change_cor_mat(Sigma1.x, affected.dims, draw_rho = draw_rho)
      sigma2 <- sqrt(apply(U, 1, function(u) u %*% Sigma2.x %*% u))
      hellinger.sim[, b] <- round(calc_hellinger(mu1, sigma1, mu1, sigma2), 6)
    }
  }
  
  mean.h <- rowMeans(hellinger.sim, na.rm = TRUE)
  
  if (return.all) {
    # quantile.h <- apply(hellinger.sim, 1, quantile, probs = 1:19/20, na.rm = TRUE)
    # sd.h <- apply(hellinger.sim, 1, sd, na.rm = TRUE)
    # mean.h.ci <- matrix(c(mean.h - 1.96 * sd.h/sqrt(n.sim), 
    #                       mean.h + 1.96 * sd.h/sqrt(n.sim)),
    #                       nrow = length(mean.h))
    return(list('est'   = mean.h,
                # 'quantile' = quantile.h,
                # 'est.ci' = mean.h.ci,
                # 'se' = sd.h / sqrt(n.sim),
                'sims' = hellinger.sim,
                'sparsity' = k,
                'type' = change.type))
  } else {
    return(mean.h)
  }
}

which_pa <- function(est.h.obj, 
                     keep.prop = 1, 
                     max.axes = length(est.h.obj$est),
                     return.all = TRUE) {
  get_prop <- function(h.sim) {
    top.PA <- apply(h.sim, 2, which.max)
    prop.top <- table(top.PA) / length(top.PA)
    prop.top.full <- rep(0, nrow(h.sim))
    ind <- as.integer(names(prop.top))
    prop.top.full[ind] <- prop.top
    prop.top.full
  }
  
  get_axes <- function(prop.top, keep.prop) {
    order.axes <- order(prop.top, decreasing = TRUE)
    cum.prop <- cumsum(prop.top[order.axes])
    n.keep <- min(sum(cum.prop <= keep.prop) + 1, max.axes)
    order.axes[1:n.keep]
  }
  
  h.sim <- est.h.obj$sims
  prop.top.avg <- get_prop(h.sim)
  top.axes.avg <- get_axes(prop.top.avg, keep.prop)
  avg.df <- data.frame('axes' = top.axes.avg,
                       'prop.max' = prop.top.avg[top.axes.avg])
  
  # To avoid unnecessary computing if not everything is to be returned.
  if (return.all) {
    type <- est.h.obj$type
    if (any(type == 'mu')) h.sim.mu <- h.sim[, type == 'mu']
    if (any(type == 'sigma')) h.sim.sigma <- h.sim[, type == 'sigma']
    if (any(type == 'rho')) h.sim.rho <- h.sim[, type == 'rho']
    
    prop.top.mu <- get_prop(h.sim.mu)
    prop.top.sigma <- get_prop(h.sim.sigma)
    prop.top.rho <- get_prop(h.sim.rho)
    if (keep.prop == 1) {
      return(data.frame('mu' = prop.top.mu,
                        'sigma' = prop.top.sigma,
                        'rho' = prop.top.rho,
                        'avg' = prop.top.avg))
    } else {
      top.axes.mu <- get_axes(prop.top.mu, keep.prop)
      top.axes.sigma <- get_axes(prop.top.sigma, keep.prop)
      top.axes.rho <- get_axes(prop.top.rho, keep.prop)
      mu.df <- data.frame('axes' = top.axes.mu,
                          'prop.max' = prop.top.mu[top.axes.mu])
      sigma.df <- data.frame('axes' = top.axes.sigma,
                             'prop.max' = prop.top.sigma[top.axes.sigma])
      rho.df <- data.frame('axes' = top.axes.rho,
                           'prop.max' = prop.top.rho[top.axes.rho])
      return(list('mu' = mu.df,
                  'sigma' = sigma.df,
                  'rho' = rho.df,
                  'avg' = avg.df))
    }
  } else {
    if (keep.prop == 1) {
      return(data.frame('axes' = 1:nrow(h.sim),
                        'prop.max' = prop.top.avg))
    } else return(avg.df)
  }
}
