library(cmdstanr)
sim1data <- readRDS("w_wout/data/sim1.rds")
mw <- cmdstan_model("w_wout/w.stan")
mwo <- cmdstan_model("w_wout/wout.stan")

mw_fit <- mw$sample(data = sim1data$standata, chains = 2, parallel_chains = 2)
mwout_fit <- mwout$sample(data = sim1data$standata, chains = 2, parallel_chains = 2)

