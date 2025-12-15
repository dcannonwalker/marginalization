library(posterior)
library(dplyr)
library(ggplot2)
library(rjags)
library(cmdstanr)

sim_list <- readRDS("cs5_2/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs5_2/data/stan_fit.rds")
jags_post <- readRDS("cs5_2/data/jags_fit.rds")

stan_smry <- stan_post$summary()
jags_smry <- summary(jags_post)

stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(mean, median, sd)
jags_p <- cbind(jags_smry$statistics[1:G, 1:2], jags_smry$quantiles[1:G, 3])
colnames(jags_p) <- c("mean", "sd", "median")
p_df <- rbind(stan_p, jags_p) %>%
  mutate(model = rep(c("stan", "jags"), each = G),
         variable = rep(paste0("p[", 1:G, "]"), 2),
         true = rep(sim_list$D_g, 2),
         b1 = rep(sim_list$b1_g * sim_list$D_g, 2))

stan_b <- stan_smry %>%
  filter(grepl("^b", variable)) %>%
  select(variable, mean, sd)
jags_b <- as.data.frame(jags_smry$statistics[(G + 1):(3 * G), 1:2])
colnames(jags_b) <- c("mean", "sd")
jags_b$variable <- rownames(jags_b)
b_df <- rbind(stan_b, jags_b) %>%
  mutate(model = rep(c("stan", "jags"), each = 2 * G),
         variable = rep(rep(c("b0", "b1"), each = G), 2),
         true_class = rep(sim_list$D_g, 4),
         est_class = c(rep(p_df$mean[1:G], 2), rep(p_df$mean[(G + 1):(2 * G)], 2)), 
         true = rep(c(sim_list$b0_g, sim_list$b1_g * sim_list$D_g), 2)) %>%
  mutate(mean = if_else(variable == "b1", mean * est_class, mean))

gg1 <- ggplot(p_df, aes(b1, mean, color = model)) + geom_point(alpha = 0.2) + 
  theme_minimal() + 
  ggtitle("Comparison of class probabilities", subtitle = "JAGS model without marginalization, Stan model with marginalization") + 
  xlab("True treatment effect") +
  ylab("Estimated class probability") + 
  ylim(c(0, 1)) + 
  facet_wrap(.~model)
gg1
gg2 <- ggplot(b_df, aes(true, mean, color = model)) + 
  geom_point(alpha = 0.2) + 
  theme_minimal() + 
  facet_wrap(.~variable) + 
  geom_abline() + 
  ggtitle("Comparison of true and estimated regression parameters") +
  xlab("True parameter value") +
  ylab("Estimated parameter value")
gg2
saveRDS(gg1, "cs5_2/data/gg1.rds")
saveRDS(gg2, "cs5_2/data/gg2.rds")
