set.seed(1)

# general
G <- 250 
N_g <- 10

# fixed effects and related
x_g <- rep(c(0, 1), each = N_g / 2)
pi0 <- 0.8 # remember, this is the probability that D_g is 0
# i.e., that there is not trt effect for g
# M <- 0.7
D_g <- sample(c(0, 1), G, replace = TRUE, prob = c(pi0, 1 - pi0))
b0_g <- rnorm(G, mean = 4)
b1_g <- sample(c(-1, 1), size = G, replace = TRUE) * log(1.5 + rexp(G))
x <- rep(x_g, G)
D <- rep(D_g, each = N_g)
b1 <- rep(b1_g, each = N_g)
b0 <- rep(b0_g, each = N_g)

# random effects and related
pairs <- rep(1:(N_g / 2), 2)
z_g <- model.matrix(~0 + as.factor(pairs))
sig_u <- 1 
u_g <- mvtnorm::rmvnorm(G, mean = rep(0, N_g / 2), sigma = diag(rep(sig_u^2, N_g / 2)))
z <- kronecker(diag(rep(1, G)), z_g)
u <- c(t(u_g))
# check that the ordering is correct: 
z %*% u # looks right

# normalization factors
sig_S <- 0.2
S <- rnorm(N_g, sd = sig_S)
sample <- factor(1:10)
sample_design_g <- model.matrix(~0+sample)
sample_design <- kronecker(rep(1, G), sample_design_g)
sample_effects_g <- sample_design_g %*% S
sample_effects <- sample_design %*% S

# mean
mu_g <- cbind(b0_g, b0_g + D_g * x_g * b1_g)
mu_no_sample <- b0 + x * D * b1 + z %*% u
mu <- mu_no_sample + sample_effects

head(mu_no_sample)
head(mu)

# bcv (sqrt dispersion, or size = 1 / bcv^2)
bcv_g <- runif(G, 0.01, 0.4)
bcv <- rep(bcv_g, each = N_g)


# response
y  <- rnbinom(N_g * G, size = 1 / bcv^2, mu = exp(mu))
y_g <- matrix(nrow = G, ncol = N_g)
G_i <- rep(seq(1, G), each = N_g)
S_i <- rep(1:10, G)
for (i in 1:G) {
  y_g[i, ] <- y[G_i == i]
}

sim_list <- list(y = y, y_g = y_g, x = x, x_g = x_g, pi0 = pi0, N_g = N_g, G = G, 
                 G_i = G_i, 
                 N = G * N_g, 
                 D = D, 
                 b0 = b0, b1 = b1,
                 D_g = D_g, b0_g = b0_g, b1_g = b1_g,
                 u = u, u_g = u_g, z_g = z_g, z = z, sig_u = sig_u, N_u = G * N_g / 2,
                 N_u_g = N_g / 2,
                 mu_no_sample = mu_no_sample,
                 sample_effects_g = sample_effects_g,
                 sample_effects = sample_effects,
                 sample_design_g = sample_design_g,
                 sample_design = sample_design,
                 S_i = S_i, 
                 sig_S = sig_S,
                 bcv_g = bcv_g,
                 bcv = bcv)
saveRDS(sim_list, "cs_eR5/data/sim_list.rds")
