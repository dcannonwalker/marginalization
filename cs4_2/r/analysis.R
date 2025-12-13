library(posterior)
library(dplyr)
library(ggplot2)
library(rjags)
library(cmdstanr)

sim_list <- readRDS("cs4_2/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs4_2/data/stan_fit.rds")
jags_post <- readRDS("cs4_2/data/jags_fit.rds")

stan_smry <- stan_post$summary()
jags_smry <- summary(jags_post)

stan_b <- stan_smry %>%
  filter(grepl("^b", variable)) %>%
  select(variable, mean, sd)
jags_b <- as.data.frame(jags_smry$statistics[, 1:2])
colnames(jags_b) <- c("mean", "sd")
jags_b$variable <- rownames(jags_b)
b_df <- rbind(stan_b, jags_b) %>%
  mutate(model = rep(c("stan", "jags"), each = 2 * G),
         variable = rep(rep(c("b0", "b1"), each = G), 2),
         true = rep(c(sim_list$b0_g, sim_list$b1_g), 2)) 

gg2 <- ggplot(b_df, aes(true, mean, color = model)) + 
  geom_point(alpha = 0.2) + 
  theme_minimal() + 
  facet_wrap(.~variable) + 
  geom_abline() + 
  ggtitle("Comparison of true and estimated regression parameters") +
  xlab("True parameter value") +
  ylab("Estimated parameter value")
gg2
saveRDS(gg2, "cs4_2/data/gg2.rds")