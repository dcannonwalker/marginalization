# get data set
data(kidney, package = 'SimSeq')

# simulate a data set
sim_data_set <- function(params, G, N, p_diff_exp = 0.2, p_up = 0.5, basefc = 1.5, 
                         sample_effects = TRUE, sig_u = 0.3) {
  idx <- sample(1:length(params$total$mean), size = G)
  mu1 <- params$total$mean[idx]
  mu2 <- mu1
  n_diff_exp <- round(G * p_diff_exp)
  upidx <- 1:round(n_diff_exp * p_up)
  dnidx <- (max(upidx) + 1):n_diff_exp
  fc <- basefc + rexp(n = n_diff_exp, rate = 1)
  mu2[upidx] <- mu2[upidx]*fc[upidx]
  mu2[dnidx] <- mu2[dnidx]/fc[dnidx]
  disp <- params$total$disp
  phi <- disp[idx]
  
  counts <- matrix(nrow = G, ncol = N)
  if (sample_effects == TRUE) {
    S <- runif(N, min = 0.7, max = 1.3)
  } else {
    S <- rep(1, N)
  }
  
  pairs <- rep(1:(N / 2), 2)
  z_g <- model.matrix(~0 + as.factor(pairs))
  u_g <- mvtnorm::rmvnorm(G, mean = rep(0, N / 2), sigma = diag(rep(sig_u^2, N / 2)))
  ranef <- t(z_g %*% t(u_g))
  
  for (i in 1:G) {
    mutemp <- rep(c(mu1[i], mu2[i]), each = N / 2) * S * exp(ranef[i, ])
    sizetemp <- rep(1 / phi[i], N)
    counts[i, ] <- rnbinom(N, size = sizetemp, mu = mutemp)
  }
  list(counts = counts, mu1 = mu1, mu2 = mu2, phi = phi, S = S, ranef = ranef, sig_u = sig_u, z_g = z_g, fc = fc, u_g = u_g)
}

# kparams <- get_params(counts = kidney$counts, groups = kidney$treatment, mean_count_filter = 10)
kparams <- readRDS("cs_byn_ranef/data/kparams.rds")
N_g <- 10
G <- 1000
x_g <- rep(c(0, 1), each = N_g / 2)
pi0 <- 0.8
sig_S <- 0.2
sample_design_g <- factor(1:N_g)
ds <- sim_data_set(params = kparams, G = G, N = N_g, basefc = 1.5)

sim_list <- list(
  N_g = N_g,
  G = G,
  y_g = ds$counts,
  x_g = x_g,
  pi0 = pi0,
  K = N_g / 2,
  z_g = ds$z_g,
  sig_u = ds$sig_u,
  sig_S = sig_S,
  sample_design_g = sample_design_g,
  ds = ds,
  kparams = kparams
)
saveRDS(sim_list, "cs_byn_ranef/data/sim_list.rds")