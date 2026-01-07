library(rjags)
load.module("mix")
set.seed(3)
model_string <- "model {
# likelihood
for (i in 1:N) {
  y[i] ~ dpois(lambda[i])
  mu[i] <- b0[G_i[i]] + D[G_i[i]] * b1[G_i[i]] * x[i] + z[i, ] %*% u
  lambda[i] <- exp(mu[i])
}

# priors

tau_b1 <- 1 / sig_b1^2
for (i in 1:G) {
  D[i] ~ dbern(pi0)
  b0[i] ~ dnorm(mu_b0, 1 / sig_b0^2)
  b1[i] ~ dnormmix(c(mu_b11, mu_b12), c(tau_b1, tau_b1), c(0.5, 0.5))
}

for (i in 1:N_u) {
  u[i] ~ dnorm(0, 1 / sig_u^2)
}

mu_b0 ~ dnorm(0, 1 / 10)
mu_b11 ~ dnorm(0, 1 / 10)
mu_b12 ~ dnorm(0, 1 / 10)
sig_b0 ~ dunif(0, 1000)
sig_b1 ~ dunif(0, 1000)
# expects y, pi0, x, G, G_i, N, sig_u
}"

sim_list <- readRDS("cs9/data/sim_list.rds")
jags_names <- c("y", "x", "pi0", "N", "G", "G_i", "sig_u", "z", "N_u")
data_list <- sim_list[jags_names]
model <- jags.model(file = textConnection(model_string), data = data_list)
update(model, n.iter = 1e3)
vn <- c("b0", "b1", "D", "sig_b0", "sig_b1", "u", "mu_b0", "mu_b11", "mu_b12")
post <- coda.samples(model, variable.names = vn, n.iter = 1e3)
saveRDS(post, "cs10/data/jags_fit.rds")

