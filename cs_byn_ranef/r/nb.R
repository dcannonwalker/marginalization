library(cmdstanr)
model <- cmdstan_model("cs_byn_ranef/stan/nb.stan")
stan_names <- names(model$variables()$data)
sim_list <- readRDS("cs_byn_ranef/data/sim_list.rds")
data_list <- sim_list[stan_names[stan_names %in% names(sim_list)]]
data_list$two_comp_mu_b1 <- 0
data_list$sample_design_g <- model.matrix(~0+sim_list$sample_design_g)
phi <- 1 / sim_list$ds$phi
phi_bounds <- c(0.5 * min(phi), 2 * max(phi))
data_list$phi_bounds <- phi_bounds
# data_list$pi0 <- 0.8
post <- model$sample(data = data_list, seed = 2, parallel_chains = 1, chains = 1, iter_warmup = 1e3)
post$save_object("cs_byn_ranef/data/nb_fit_8.rds")
smry <- post$summary()
saveRDS(smry, "cs_byn_ranef/data/nb_smry_8.rds")
