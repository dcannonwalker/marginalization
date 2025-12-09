set.seed(1)
N <- 1e3
x <- rep(c(0, 1), each = N / 2)
D <- sample(c(0, 1), N, replace = TRUE, prob = c(0.5, 0.5))
b0 <- 0
b1 <- 2
tau <- 1
mu <- b0 + D * x * b1
y  <- rnorm(N, mean = mu, sd = sqrt(1 / tau))
sim_list <- list(y = y, tau = tau, x = x, pi0 = 0.5, N = N, sig = sqrt(1 / tau), D = D, b0 = b0, b1 = b1)
saveRDS(sim_list, "case_study_1/data/sim_list.rds")