library(cmdstanr)
m <- cmdstan_model("cbp_model.stan")
names(m$variables()$data)
