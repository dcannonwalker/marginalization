set.seed(1)
N_g <- 10
G <- 100
x_g <- rep(c(0, 1), each = N_g / 2)

pi0 <- 0.8 # remember, this is the probability that D_g is 0
# i.e., that there is not trt effect for g
D_g <- sample(c(0, 1), G, replace = TRUE, prob = c(pi0, 1 - pi0))
b0_g <- rnorm(G, mean = 4)
b1_g <- sample(c(-1, 1), size = G, replace = TRUE) * log(1.5 + rexp(G))
x <- rep(x_g, G)
D <- rep(D_g, each = N_g)
b1 <- rep(b1_g, each = N_g)
b0 <- rep(b0_g, each = N_g)

# normalization factors
sig_S <- 0.2
S <- rnorm(N_g, sd = sig_S)
sample <- factor(1:N_g)
sample_design_g <- model.matrix(~0+sample)
sample_design <- kronecker(rep(1, G), sample_design_g)
sample_effects_g <- sample_design_g %*% S
sample_effects <- sample_design %*% S

# mean
mu_g <- cbind(b0_g, b0_g + D_g * x_g * b1_g)
mu_no_sample <- b0 + x * D * b1 
mu <- mu_no_sample + sample_effects

head(mu_no_sample)
head(mu)

# response
y  <- rpois(N_g * G, lambda = exp(mu))
y_g <- matrix(nrow = G, ncol = N_g)
G_i <- rep(seq(1, G), each = N_g)
S_i <- rep(1:N_g, G)
for (i in 1:G) {
  y_g[i, ] <- y[G_i == i]
}

sim_list <- list(y = y, y_g = y_g, x = x, x_g = x_g, pi0 = pi0, N_g = N_g, G = G, 
                 G_i = G_i, 
                 N = G * N_g, 
                 D = D, 
                 b0 = b0, b1 = b1,
                 D_g = D_g, b0_g = b0_g, b1_g = b1_g,
                 mu_no_sample = mu_no_sample,
                 sample_effects_g = sample_effects_g,
                 sample_effects = sample_effects,
                 sample_design_g = sample_design_g,
                 sample_design = sample_design,
                 S_i = S_i, 
                 sig_S = sig_S)
saveRDS(sim_list, "cs_s1/data/sim_list.rds")