library(cmdstanr)
sim_list <- readRDS("cs_eR5/data/sim_list.rds")
model <- cmdstan_model("cs_eR5/stan/nb.stan")
stan_names <- names(model$variables()$data)
data_list <- sim_list[stan_names[stan_names %in% names(sim_list)]]
data_list$two_comp_mu_b1 <- 0
phi <- 1 / sim_list$bcv_g^2
phi_bounds <- c(0.5 * min(phi), 2 * max(phi))
data_list$phi_bounds <- phi_bounds

post <- model$sample(data = data_list, seed = 2, parallel_chains = 1, chains = 1, iter_warmup = 5e3)
post$save_object("cs_eR5/data/nb_fit.rds")
stan_smry <- post$summary()
saveRDS(stan_smry, "cs_eR5/data/nb_smry.rds")
