library(MASS)
G <- 20
N <- 10
K <- N / 2
b0 <- rnorm(G, mean = 4)
b1 <- rnorm(G, sd = 2)
u <- mvtnorm::rmvnorm(G, mean = rep(0, K))
x <- cbind(1, rep(c(0, 1), each = 5))
pairs <- factor(rep(1:K, 2))
z <- model.matrix(~0 + pairs)
log_S <- kronecker(t(rep(1, G)), rnorm(N))
S <- t(exp(log_S))
log_mu <- t(x %*% t(b)) + t(z %*% t(u)) + 
  t(model.matrix(~0 + factor(1:N)) %*% log_S)
mu <- exp(log_mu)
phi <- runif(G, 0.5, 2)
counts_nb <- t(sapply(1:G, function(g) {
  rnbinom(N, size = phi[g], mu = mu[g, ])
}))
counts_pois <- t(sapply(1:G, function(g) {
  rpois(N, lambda = mu[g, ])
}))

fit_nb_nb <- lapply(1:G, function(g) {
  ff <- glm.nb(counts_nb[g, ] ~ x[, 2] + pairs + offset(log_S[, 1]))
  list(fit = ff, aic = ff$aic)
})
names(fit_nb_nb) <- paste0("gene", 1:G)

r <- counts_nb[3, ]
offset(S[1, ])
glm.nb(r ~ x[, 2] + pairs + offset(log_S[, 1]))
