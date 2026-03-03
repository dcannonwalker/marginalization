library(dplyr)
library(posterior)
m <- cmdstanr::cmdstan_model("w_wout/wout.stan")
names(m$variables()$data)
d <- list(
  N = 10, 
  G = 100,
  P = 5,
  x = rep(c(0, 1), each = 5),
  z = model.matrix(~0 + factor(rep(1:5, 2))),
  po = 1,
  y = numeric(0),
  mao = 4
)
po <- m$sample(data = d, chains = 1)
dr <- po$draws()
s <- summary(dr)
n <- paste0("ysim[", rep(1:d$G, d$N), ",", rep(1:d$N, each = d$G), "]")
nphi <- paste0("bcv[", 1:d$G, "]")
hist(s %>%
       filter(variable %in% nphi) %>%
       pull(mean), breaks = 20)
bcv <- extract_variable_array(dr, "bcv")
ysim <- extract_variable_array(dr, "ysim")
ysim1 <- ysim[1, 1, , ]

dl <- d
dl$po <- 0
dl$y <- ysim1
l <- m$sample(data = dl, chains = 1)
