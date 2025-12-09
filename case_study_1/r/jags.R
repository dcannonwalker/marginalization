library(rjags)
set.seed(3)
model_string <- "model {
# likelihood
for (i in 1:N) {
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- b0 + D[i] * b1 * x[i]
}

# priors
b0 ~ dnorm(0, 1)
b1 ~ dnorm(0, 1)
for (i in 1:N) {
  D[i] ~ dbern(pi0)
}
# expects y, tau, pi0, x
}"

sim_list <- readRDS("case_study_1/data/sim_list.rds")
jags_names <- c("y", "tau", "x", "pi0", "N")
data_list <- sim_list[jags_names]
model <- jags.model(file = textConnection(model_string), data = data_list)
update(model, n.iter = 1000)
vn <- c("b0", "b1", "D")
# vn <- "D"
post <- coda.samples(model, variable.names = vn, n.iter = 1000)
saveRDS(post, "case_study_1/data/jags_fit.rds")

