library(rjags)
set.seed(3)
model_string <- "model {
# likelihood
for (i in 1:N) {
  y[i] ~ dnorm(mu[i], 1 / sig^2)
  mu[i] <- b0[G_i[i]] + D[G_i[i]] * b1[G_i[i]] * x[i] + z[i, ] %*% u
}

# priors
for (i in 1:G) {
  D[i] ~ dbern(pi0)
  b0[i] ~ dnorm(0, 1 / sig_b0^2)
  b1[i] ~ dnorm(0, 1 / sig_b1^2)
}

for (i in 1:N_u) {
  u[i] ~ dnorm(0, 1 / sig_u^2)
}
sig ~ dunif(0, 1000)
sig_b0 ~ dunif(0, 1000)
sig_b1 ~ dunif(0, 1000)
# expects y, pi0, x, G, G_i, N, sig_u
}"

sim_list <- readRDS("cs8/data/sim_list.rds")
jags_names <- c("y", "x", "pi0", "N", "G", "G_i", "sig_u", "z", "N_u")
data_list <- sim_list[jags_names]
model <- jags.model(file = textConnection(model_string), data = data_list)
update(model, n.iter = 1e4)
vn <- c("b0", "b1", "D", "sig", "sig_b0", "sig_b1", "u")
post <- coda.samples(model, variable.names = vn, n.iter = 1e3)
saveRDS(post, "cs8/data/jags_fit.rds")

