library(posterior)
library(dplyr)
library(ggplot2)
library(rjags)
library(cmdstanr)

sim_list <- readRDS("case_study_1/data/sim_list.rds")
stan_post <- readRDS("case_study_1/data/stan_fit.rds")
jags_post <- readRDS("case_study_1/data/jags_fit.rds")

stan_smry <- stan_post$summary()
jags_smry <- summary(jags_post)

stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(mean, sd)
jags_p <- jags_smry$statistics[1:1000, 1:2]
colnames(jags_p) <- c("mean", "sd")
p_df <- rbind(stan_p, jags_p) %>%
  mutate(model = rep(c("stan", "jags"), each = 1000),
         variable = rep(paste0("p[", 1:1000, "]"), 2),
         true = rep(sim_list$D, 2),
         y = rep(sim_list$y, 2))

stan_b <- stan_smry %>%
  filter(grepl("^b", variable)) %>%
  select(mean, sd)
jags_b <- jags_smry$statistics[1001:1002, 1:2]
colnames(jags_b) <- c("mean", "sd")
b_df <- rbind(stan_b, jags_b) %>%
  mutate(model = rep(c("stan", "jags"), each = 2),
         variable = rep(c("b0", "b1"), 2))

gg1 <- ggplot(p_df, aes(y, mean, color = model)) + geom_point(alpha = 0.2) + 
  theme_minimal() + 
  ggtitle("Comparison of class probabilities", subtitle = "JAGS model without marginalization, Stan model with marginalization") + 
  xlab("Observed value") +
  ylab("Estimated class probability")
gg1
saveRDS(gg1, "case_study_1/data/gg1.rds")

