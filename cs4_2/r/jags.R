library(rjags)
set.seed(3)
model_string <- "model {
# likelihood
for (i in 1:N) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- b0[G_i[i]] + b1[G_i[i]] * x[i]
}

# priors
for (i in 1:G) {
  b0[i] ~ dnorm(0, 1 / 10)
  b1[i] ~ dnorm(0, 1 / 10)
}
# expects y, tau, x, G, G_i, N
}"

sim_list <- readRDS("cs4_2/data/sim_list.rds")
jags_names <- c("y", "tau", "x", "N", "G", "G_i")
data_list <- sim_list[jags_names]
model <- jags.model(file = textConnection(model_string), data = data_list)
update(model, n.iter = 1e4)
vn <- c("b0", "b1")
post <- coda.samples(model, variable.names = vn, n.iter = 1e3)
saveRDS(post, "cs4_2/data/jags_fit.rds")

