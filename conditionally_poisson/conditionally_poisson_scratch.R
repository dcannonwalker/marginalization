library(MASS)
library(dplyr)
set.seed(1)
G <- 100
N <- 10
K <- N / 2
b0 <- rnorm(G, mean = 6)
b1 <- rnorm(G, sd = 1)
b <- cbind(b0, b1)
u <- sqrt(0.5) * mvtnorm::rmvnorm(G, mean = rep(0, K))
x <- cbind(1, rep(c(0, 1), each = 5))
pairs <- factor(rep(1:K, 2))
z <- model.matrix(~0 + pairs)
log_S <- kronecker(t(rep(1, G)), rnorm(N, sd = 0.25))
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

# nb simulation
fit_nb_nb <- lapply(1:G, function(g) {
  ff <- glm.nb(counts_nb[g, ] ~ x[, 2] + pairs + offset(log_S[, 1]))
  list(fit = ff, aic = ff$aic)
})
names(fit_nb_nb) <- paste0("gene", 1:G)

fit_pois_nb <- lapply(1:G, function(g) {
  ff <- glm(counts_nb[g, ] ~ x[, 2] + pairs + offset(log_S[, 1]), family = poisson)
  list(fit = ff, aic = ff$aic)
})
names(fit_pois_nb) <- paste0("gene", 1:G)

aic_nb_nb <- sapply(fit_nb_nb, function(ff) ff$aic)
aic_pois_nb <- sapply(fit_pois_nb, function(ff) ff$aic)

augment_aic <- function(aic1, aic2, cnames) {
  probs <- t(apply(cbind(aic1, aic2), 1, function(aic) exp(-(aic - min(aic)) / 2)))
  colnames(probs) <- paste0("prob_", cnames)
  out <- cbind(aic1, aic2)
  colnames(out) <- cnames
  out <- as_tibble(cbind(out, probs)) 
}
aic_nb <- augment_aic(aic_nb_nb, aic_pois_nb, c("nb", "pois")) %>% 
  mutate(b0 = b0, b1 = b1, phi = phi,
         preferred_model = factor(if_else(prob_nb > prob_pois, "nb", "pois"))) %>%
  arrange(prob_pois)
summary(aic_nb$preferred_model)

# pois simulation
fit_nb_pois <- lapply(1:G, function(g) {
  ff <- glm.nb(counts_pois[g, ] ~ x[, 2] + pairs + offset(log_S[, 1]))
  list(fit = ff, aic = ff$aic)
})
names(fit_nb_pois) <- paste0("gene", 1:G)

fit_pois_pois <- lapply(1:G, function(g) {
  ff <- glm(counts_pois[g, ] ~ x[, 2] + pairs + offset(log_S[, 1]), family = poisson)
  list(fit = ff, aic = ff$aic)
})
names(fit_pois_pois) <- paste0("gene", 1:G)

aic_nb_pois <- sapply(fit_nb_pois, function(ff) ff$aic)
aic_pois_pois <- sapply(fit_pois_pois, function(ff) ff$aic)

augment_aic <- function(aic1, aic2, cnames) {
  probs <- t(apply(cbind(aic1, aic2), 1, function(aic) exp(-(aic - min(aic)) / 2)))
  colnames(probs) <- paste0("prob_", cnames)
  out <- cbind(aic1, aic2)
  colnames(out) <- cnames
  out <- as_tibble(cbind(out, probs)) 
}
aic_pois <- augment_aic(aic_nb_pois, aic_pois_pois, c("nb", "pois")) %>% 
  mutate(b0 = b0, b1 = b1, phi = phi,
         preferred_model = factor(if_else(prob_nb > prob_pois, "nb", "pois"))) %>%
  arrange(prob_pois)
summary(aic_pois$preferred_model)
saveRDS(list(
  fit_nb_nb = fit_nb_nb,
  fit_nb_pois = fit_nb_pois,
  fit_pois_nb = fit_pois_nb,
  fit_pois_pois = fit_pois_pois,
  aic_nb = aic_nb,
  aic_pois = aic_pois
), "aic_w_pairs.rds")

# ======== simplified model ==========
# nb simulation
fit0_nb_nb <- lapply(1:G, function(g) {
  ff <- glm.nb(counts_nb[g, ] ~ x[, 2] + offset(log_S[, 1]))
  list(fit = ff, aic = ff$aic)
})
names(fit0_nb_nb) <- paste0("gene", 1:G)

fit0_pois_nb <- lapply(1:G, function(g) {
  ff <- glm(counts_nb[g, ] ~ x[, 2] + offset(log_S[, 1]), family = poisson)
  list(fit = ff, aic = ff$aic)
})
names(fit0_pois_nb) <- paste0("gene", 1:G)

aic0_nb_nb <- sapply(fit0_nb_nb, function(ff) ff$aic)
aic0_pois_nb <- sapply(fit0_pois_nb, function(ff) ff$aic)

aic0_nb <- augment_aic(aic0_nb_nb, aic0_pois_nb, c("nb", "pois")) %>% 
  mutate(b0 = b0, b1 = b1, phi = phi,
         preferred_model = factor(if_else(prob_nb > prob_pois, "nb", "pois"))) %>%
  arrange(prob_pois)
summary(aic0_nb$preferred_model)

# pois simulation
fit0_nb_pois <- lapply(1:G, function(g) {
  ff <- glm.nb(counts_pois[g, ] ~ x[, 2] + offset(log_S[, 1]))
  list(fit = ff, aic = ff$aic)
})
names(fit0_nb_pois) <- paste0("gene", 1:G)

fit0_pois_pois <- lapply(1:G, function(g) {
  ff <- glm(counts_pois[g, ] ~ x[, 2] + offset(log_S[, 1]), family = poisson)
  list(fit = ff, aic = ff$aic)
})
names(fit0_pois_pois) <- paste0("gene", 1:G)

aic0_nb_pois <- sapply(fit0_nb_pois, function(ff) ff$aic)
aic0_pois_pois <- sapply(fit0_pois_pois, function(ff) ff$aic)

aic0_pois <- augment_aic(aic0_nb_pois, aic0_pois_pois, c("nb", "pois")) %>% 
  mutate(b0 = b0, b1 = b1, phi = phi,
         preferred_model = factor(if_else(prob_nb > prob_pois, "nb", "pois"))) %>%
  arrange(prob_pois)
summary(aic0_pois$preferred_model)

saveRDS(list(
  fit0_nb_nb = fit0_nb_nb,
  fit0_nb_pois = fit0_nb_pois,
  fit0_pois_nb = fit0_pois_nb,
  fit0_pois_pois = fit0_pois_pois,
  aic0_nb = aic0_nb,
  aic0_pois = aic0_pois
), "aic_wo_pairs.rds")
