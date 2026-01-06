library(posterior)
library(dplyr)
library(ggplot2)
library(rjags)
library(cmdstanr)

sim_list <- readRDS("cs12/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs12/data/stan_fit.rds")
jags_post <- readRDS("cs12/data/jags_fit.rds")

stan_smry <- stan_post$summary()
jags_smry <- summary(jags_post)

# ============= various sig_ ================
jags_sig_b0 <- jags_post[[1]][, which(dimnames(jags_post[[1]])[[2]] == "sig_b0")]
stan_sig_b0 <- extract_variable_matrix(stan_post, "sig_b0")
plot(density(jags_sig_b0))
plot(density(stan_sig_b0))

jags_sig_b1 <- jags_post[[1]][, which(dimnames(jags_post[[1]])[[2]] == "sig_b1")]
stan_sig_b1 <- extract_variable_matrix(stan_post, "sig_b1")
plot(density(jags_sig_b1))
plot(density(stan_sig_b1))

# ============= various mu ================
jags_mu_b0 <- jags_post[[1]][, which(dimnames(jags_post[[1]])[[2]] == "mu_b0")]
stan_mu_b0 <- extract_variable_matrix(stan_post, "mu_b0")
plot(density(jags_mu_b0))
plot(density(stan_mu_b0))

jags_mu_b11 <- jags_post[[1]][, which(dimnames(jags_post[[1]])[[2]] == "mu_b11")]
stan_mu_b11 <- extract_variable_matrix(stan_post, "mu_b11")
plot(density(jags_mu_b11))
plot(density(stan_mu_b11))

jags_mu_b12 <- jags_post[[1]][, which(dimnames(jags_post[[1]])[[2]] == "mu_b12")]
stan_mu_b12 <- extract_variable_matrix(stan_post, "mu_b12")
plot(density(jags_mu_b12))
plot(density(stan_mu_b12))

# ============= S ================
jags_S <- jags_post[[1]][, which(grepl("S", dimnames(jags_post[[1]])[[2]]))]
stan_S <- extract_variable_array(stan_post, variable = "S")[, 1, ]
plot(density(jags_S[, 5]))
plot(density(stan_S[, 8]))

# ================ p ===================

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

# =============== b ====================
stan_b <- stan_smry %>%
  filter(grepl("b1\\[|b0\\[", variable)) %>%
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

# ============== u ====================
stan_u <- stan_smry %>%
  filter(grepl("^u_g", variable)) %>%
  select(variable, mean, sd) %>%
  mutate(model = "stan")
jags_u <- as.data.frame(jags_smry$statistics[(3 * G + 6):(8 * G + 5), 1:2])
colnames(jags_u) <- c("mean", "sd")
u_name <- paste0("u_g[", rep(1:50, each = 5), ",", rep(1:5, 50), "]")
jags_u$variable <- u_name
jags_u <- tibble(jags_u) %>%
  mutate(model = "jags")
true_u <- tibble(true = sim_list$u, variable = u_name)
u_df <- rbind(stan_u, jags_u) %>%
  left_join(true_u)

# ============== plots ==================
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
  facet_wrap(.~variable, ncol = 1, nrow = 2) + 
  geom_abline() + 
  ggtitle("Comparison of true and estimated regression parameters") +
  xlab("True parameter value") +
  ylab("Estimated parameter value")
gg2
saveRDS(gg1, "cs12/data/gg1.rds")
saveRDS(gg2, "cs12/data/gg2.rds")

## sig_b0
sig_b0 <- tibble(
  draw = c(seq(1, length(stan_sig_b0)), seq(1, length(jags_sig_b0))),
  model = c(rep("stan", length(stan_sig_b0)), rep("jags", length(jags_sig_b0))),
  value = c(stan_sig_b0, jags_sig_b0)
)
true_sig_b0 <- 1
gg4 <- ggplot(sig_b0, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = true_sig_b0, ) +
  xlim(0, NA) + 
  theme_minimal() + 
  ggtitle("Posterior density for sigma_b0")  +
  xlab("value") +
  annotate(x=true_sig_b0,y=+Inf,label="True sigma_b0",vjust=2,geom="label")
gg4
saveRDS(gg4, "cs12/data/gg4.rds")

## sig_b1
sig_b1 <- tibble(
  draw = c(seq(1, length(stan_sig_b1)), seq(1, length(jags_sig_b1))),
  model = c(rep("stan", length(stan_sig_b1)), rep("jags", length(jags_sig_b1))),
  value = c(stan_sig_b1, jags_sig_b1)
)
gg5 <- ggplot(sig_b1, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = 1, ) +
  annotate(x=1,y=+Inf,label="True sigma_b1",vjust=2,geom="label") +
  theme_minimal() + 
  xlim(0, NA) +
  ggtitle("Posterior density for sigma_b1", subtitle = "The two components of the normal mixture have the same variance")  +
  xlab("value")

gg5
saveRDS(gg5, "cs12/data/gg5.rds")

## u
gg6 <- ggplot(u_df, aes(mean, true, color = model)) + 
  geom_point(alpha = 0.6) +
  theme_minimal() + 
  geom_abline() +
  ggtitle("Comparison of true and estimated random effects") + 
  xlab("Posterior mean u_g") + 
  ylab("True value of u_g")
gg6
saveRDS(gg6, "cs12/data/gg6.rds")

## mu_b0
mu_b0 <- tibble(
  draw = c(seq(1, length(stan_mu_b0)), seq(1, length(jags_mu_b0))),
  model = c(rep("stan", length(stan_mu_b0)), rep("jags", length(jags_mu_b0))),
  value = c(stan_mu_b0, jags_mu_b0)
)
gg7 <- ggplot(mu_b0, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = 4, ) +
  annotate(x=4,y=+Inf,label="True mu_b0",vjust=2,geom="label") +
  theme_minimal() + 
  # xlim(0, ) +
  ggtitle("Posterior density for mu_b0")  +
  xlab("value")

gg7
saveRDS(gg7, "cs12/data/gg7.rds")

## mu_b11
mu_b11 <- tibble(
  draw = c(seq(1, length(stan_mu_b11)), seq(1, length(jags_mu_b11))),
  model = c(rep("stan", length(stan_mu_b11)), rep("jags", length(jags_mu_b11))),
  value = c(stan_mu_b11, jags_mu_b11)
)
gg8 <- ggplot(mu_b11, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = c(-2, 2), ) +
  annotate(x=c(-4, 4),y=+Inf,label="True mu_b1",vjust=2,geom="label") +
  theme_minimal() + 
  xlim(-10, 10) +
  ggtitle("Posterior density for mu_b11")  +
  xlab("value")

gg8
saveRDS(gg8, "cs12/data/gg8.rds")

## mu_b11
mu_b12 <- tibble(
  draw = c(seq(1, length(stan_mu_b12)), seq(1, length(jags_mu_b12))),
  model = c(rep("stan", length(stan_mu_b12)), rep("jags", length(jags_mu_b12))),
  value = c(stan_mu_b12, jags_mu_b12)
)
gg9 <- ggplot(mu_b12, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = c(-2, 2), ) +
  annotate(x=c(-2, 2),y=+Inf,label="True mu_b1",vjust=2,geom="label") +
  theme_minimal() + 
  xlim(-10, 10) +
  ggtitle("Posterior density for mu_b12")  +
  xlab("value")

gg9
saveRDS(gg9, "cs12/data/gg9.rds")

## S
S <- bind_rows(tibble::as_tibble(jags_S) %>%
  mutate(draw = row_number(), model = "jags") %>%
  pivot_longer(cols = -c(draw, model)),
tibble::as_tibble(stan_S) %>% 
  purrr::set_names(paste0("S[", 1:10, "]")) %>%
  mutate(draw = row_number(), model = "stan") %>%
  pivot_longer(cols = -c(draw, model)))
gg10 <- ggplot(S, aes(value, fill = model)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = c(-2, 2), ) +
  annotate(x=c(-2, 2),y=+Inf,label="True mu_b1",vjust=2,geom="label") +
  theme_minimal() + 
  xlim(-10, 10) +
  ggtitle("Posterior density for ")  +
  xlab("value")

gg9
saveRDS(gg9, "cs12/data/gg9.rds")

