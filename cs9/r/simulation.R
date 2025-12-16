set.seed(1)

# general
G <- 50 
N_g <- 10

# fixed effects and related
x_g <- rep(c(0, 1), each = N_g / 2)
pi0 <- 0.5
M <- 4
D_g <- sample(c(0, 1), G, replace = TRUE, prob = c(pi0, 1 - pi0))
b0_g <- rnorm(G, mean = 4)
b1_g <- sample(c(-1, 1), size = G, replace = TRUE) * (M + rnorm(G, sd = 1))
x <- rep(x_g, G)
D <- rep(D_g, each = N_g)
b1 <- rep(b1_g, each = N_g)
b0 <- rep(b0_g, each = N_g)

# random effects and related
pairs <- rep(1:(N_g / 2), 2)
z_g <- model.matrix(~0 + as.factor(pairs))
sig_u <- 1 
u_g <- mvtnorm::rmvnorm(G, mean = rep(0, N_g / 2), sigma = diag(rep(sig_u, N_g / 2)))
z <- kronecker(diag(rep(1, G)), z_g)
u <- c(t(u_g))
# check that the ordering is correct: 
z %*% u # looks right

# mean
mu_g <- cbind(b0_g, b0_g + D_g * x_g * b1_g)
mu <- b0 + x * D * b1 + z %*% u

# response
y  <- rpois(N_g * G, lambda = exp(mu))
y_g <- matrix(nrow = G, ncol = N_g)
G_i <- rep(seq(1, G), each = N_g)
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
                 N_u_g = N_g / 2)
saveRDS(sim_list, "cs9/data/sim_list.rds")
