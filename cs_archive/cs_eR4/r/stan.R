library(cmdstanr)
sim_list <- readRDS("cs_eR4/data/sim_list.rds")
model <- cmdstan_model("cs_eR4/stan/model.stan")
stan_names <- names(model$variables()$data)
data_list <- sim_list[stan_names[stan_names %in% names(sim_list)]]
data_list$two_comp_mu_b1 <- 0

post <- model$sample(data = data_list, seed = 2, parallel_chains = 1, chains = 1, iter_warmup = 5e3)
post$save_object("cs_eR4/data/stan_fit.rds")
stan_smry <- post$summary()
saveRDS(stan_smry, "cs_eR4/data/stan_smry.rds")
