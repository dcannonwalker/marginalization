library(cmdstanr)

mwout <- cmdstan_model("priors/wout.stan")
data(kidney, package = "SimSeq")
str(kidney)
keep <- edgeR::filterByExpr(kidney$counts, min.count = 100)
counts <- kidney$counts[keep, ]
G <- 1000
N <- 10
P <- N / 2
x <- as.numeric(kidney$treatment[1:N]) - 1
z <- model.matrix(~ 0 + factor(kidney$replic[1:N]))
po <- 0
y <- counts[1:G, 1:N]
mao <- log(mean(y))

standata <- list(
  G = G,
  N = N,
  P = P,
  x = x,
  z = z,
  po = po,
  y = y,
  mao = mao
)

fit <- mwout$sample(data = standata, chains = 1)
