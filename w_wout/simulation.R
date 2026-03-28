library(mvtnorm)
seedn <- 101
set.seed(seedn)
N <- 10
G <- 1000
P <- N / 2
x <- rep(c(0, 1), each = N / 2)
z <- model.matrix(~0 + factor(rep(1:P, 2)))
ma <- 6
sa <- 1
mb <- 0
sb <- 2
sp <- 1
su <- 0.5
alpha <- rnorm(G, ma, sa)
# beta <- rnorm(G, mb, sb)
updown <- sample(c(-1, 1), G, replace = TRUE)
beta <- updown * log(1.5 + rexp(G))
hist(beta)
bcv <- abs(rnorm(G, 0, sp))
phi <- 1 / bcv^2
prob <- 0.8
u <- rmvnorm(G, mean = rep(0, P), sigma = su * diag(P))
y <- matrix(nrow = G, ncol = N)
eta <- matrix(nrow = G, ncol = N)
D <- sample(c(0, 1), size = G, replace = TRUE, prob = c(1 - prob, prob)) # is gene i null?

for (i in 1:G) {
  # i <- 1
  eta[i, ] <- alpha[i] + (1 - D[i]) * x * beta[i] + z %*% u[i, ]
  y[i, ] <- rnbinom(N, mu = exp(eta[i, ]), size = phi[i])
}

standata <- list(
  N = N,
  P = P,
  G = G,
  x = x,
  z = z,
  po = 0,
  y = y,
  mao = ma,
  prob = prob
)

saveRDS(list(
  simlist = list(
    y = y,
    N = N,
    G = G,
    P = P,
    x = x,
    z = z,
    ma = ma,
    seedn = seedn,
    ma = ma, sa = sa, 
    mb = mb, sb = sb,
    sp = sp,
    su = su,
    alpha = alpha,
    beta = beta,
    bcv = bcv,
    phi = phi,
    u = u,
    D = D,
    eta = eta,
    y = y
  ),
  standata = standata), "w_wout/data/sim1.rds")

