library(cmdstanr)

simlist <- readRDS("w_wout_test2/data/sim_list.rds")
mwout <- cmdstan_model("w_wout_test2/stan/wout_fixed_s.stan")
mw <- cmdstan_model("w_wout_test2/stan/w_fixed_s.stan")
names(mw$variables()$data)
temp <- edgeR::DGEList(counts = simlist$y_g)
temp <- edgeR::normLibSizes(temp)
s <- log(temp$samples$lib.size * temp$samples$norm.factors)
standata <- list(
  N = simlist$N_g,
  G = simlist$G,
  P = simlist$N_g / 2,
  x = simlist$x_g,
  z = simlist$z_g,
  po = 0,
  y = simlist$y_g,
  mao = 0,
  s = s,
  prob = 0.9,
  bss = 50
)

fit_w <- mw$sample(data = standata, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
fit_w$save_object("w_wout_test2/data/fitw.rds")
fit_wout <- mwout$sample(data = standata, chains = 1, iter_warmup = 1000, iter_sampling = 1000)
fit_wout$save_object("w_wout_test2/data/fitwout.rds")

