set.seed(1)
G <- 50 
N_g <- 10
x_g <- rep(c(0, 1), each = N_g / 2)
pi0 <- 0.5
D_g <- sample(c(0, 1), G, replace = TRUE, prob = c(pi0, 1 - pi0))
b0_g <- rnorm(G)
b1_g <- sample(c(-1, 1), size = G, replace = TRUE) * (4 + rnorm(G, sd = 0.2))
tau <- 10 
mu_g <- cbind(b0_g, b0_g + D_g * x_g * b1_g)
x <- rep(x_g, G)
D <- rep(D_g, each = N_g)
b1 <- rep(b1_g, each = N_g)
b0 <- rep(b0_g, each = N_g)
mu <- b0 + x * D * b1
y  <- rnorm(N_g * G, mean = mu, sd = sqrt(1 / tau))
y_g <- matrix(nrow = G, ncol = N_g)
G_i <- rep(seq(1, G), each = N_g)
for (i in 1:G) {
  y_g[i, ] <- y[G_i == i]
}

sim_list <- list(y = y, y_g = y_g, tau = tau, x = x, x_g = x_g, pi0 = pi0, N_g = N_g, G = G, 
                 G_i = G_i, 
                 N = G * N_g, 
                 sig = sqrt(1 / tau), D = D, 
                 b0 = b0, b1 = b1,
                 D_g = D_g, b0_g = b0_g, b1_g = b1_g)
saveRDS(sim_list, "cs5/data/sim_list.rds")
