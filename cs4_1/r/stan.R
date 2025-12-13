library(cmdstanr)
sim_list <- readRDS("cs4_1/data/sim_list.rds")
model <- cmdstan_model("cs4_1/stan/model.stan")
stan_names <- names(model$variables()$data)
data_list <- sim_list[stan_names]

post <- model$sample(data = data_list, seed = 2, parallel_chains = 1, chains = 1)
post$save_object("cs4_1/data/stan_fit.rds")
