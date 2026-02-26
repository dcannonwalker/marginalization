library(bayesplot)
library(dplyr)
fit <- readRDS("cs_BYN1/data/nb_fit_8.rds")
smry <- readRDS("cs_BYN1/data/nb_smry_8.rds")
b1_names <- smry %>%
  filter(grepl("^b1\\[", variable)) %>%
  pull(variable)
sim_list <- readRDS("cs_BYN1/data/sim_list.rds")
draws <- fit$draws()
for (i in 1:(length(b1_names) / 50)) {
  i <- 7
  mcmc_areas(draws, pars = b1_names[(50 * i - 49):(50 * i)])
}
mcmc_areas(draws, pars = b1_names[1:50])
