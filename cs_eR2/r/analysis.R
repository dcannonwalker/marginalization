library(posterior)
library(dplyr)
library(ggplot2)
library(cmdstanr)

sim_list <- readRDS("cs13/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs13/data/stan_fit.rds")
eR_fit <- readRDS("cs13/data/eR_fit.rds")

stan_smry <- stan_post$summary()

# ============= various sig_ ================
stan_sig_b0 <- extract_variable_matrix(stan_post, "sig_b0")
plot(density(stan_sig_b0))

stan_sig_b1 <- extract_variable_matrix(stan_post, "sig_b1")
plot(density(stan_sig_b1))

# ============= various mu ================
stan_mu_b0 <- extract_variable_matrix(stan_post, "mu_b0")
plot(density(stan_mu_b0))

stan_mu_b11 <- extract_variable_matrix(stan_post, "mu_b11")
plot(density(stan_mu_b11))

stan_mu_b12 <- extract_variable_matrix(stan_post, "mu_b12")
plot(density(stan_mu_b12))

# ============= S ================
stan_S <- extract_variable_array(stan_post, variable = "S")[, 1, ]
plot(density(stan_S[, 8]))

# ================ p ===================

stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(mean, median, sd)
# ============ FIX =====================
# p_df <- rbind(stan_p, jags_p) %>%
#   mutate(model = rep(c("stan", "jags"), each = G),
#          variable = rep(paste0("p[", 1:G, "]"), 2),
#          true = rep(sim_list$D_g, 2),
#          b1 = rep(sim_list$b1_g * sim_list$D_g, 2))

# =============== b ====================
stan_b <- stan_smry %>%
  filter(grepl("b1\\[|b0\\[", variable)) %>%
  select(variable, mean, sd)
# ============== FIX ===================
# b_df <- rbind(stan_b, jags_b) %>%
#   mutate(model = rep(c("stan", "jags"), each = 2 * G),
#          variable = rep(rep(c("b0", "b1"), each = G), 2),
#          true_class = rep(sim_list$D_g, 4),
#          est_class = c(rep(p_df$mean[1:G], 2), rep(p_df$mean[(G + 1):(2 * G)], 2)), 
#          true = rep(c(sim_list$b0_g, sim_list$b1_g * sim_list$D_g), 2)) %>%
#   mutate(mean = if_else(variable == "b1", mean * est_class, mean))

# ============== u ====================
stan_u <- stan_smry %>%
  filter(grepl("^u_g", variable)) %>%
  select(variable, mean, sd) %>%
  mutate(model = "stan")
u_name <- paste0("u_g[", rep(1:50, each = 5), ",", rep(1:5, 50), "]")
true_u <- tibble(true = sim_list$u, variable = u_name)
# ========= FIX ============
# u_df <- rbind(stan_u, jags_u) %>%
#   left_join(true_u)

# ============== FIX ALL PLOTS ==========
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
saveRDS(gg1, "cs13/data/gg1.rds")
saveRDS(gg2, "cs13/data/gg2.rds")

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
saveRDS(gg4, "cs13/data/gg4.rds")

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
saveRDS(gg5, "cs13/data/gg5.rds")

## u
gg6 <- ggplot(u_df, aes(mean, true, color = model)) + 
  geom_point(alpha = 0.6) +
  theme_minimal() + 
  geom_abline() +
  ggtitle("Comparison of true and estimated random effects") + 
  xlab("Posterior mean u_g") + 
  ylab("True value of u_g")
gg6
saveRDS(gg6, "cs13/data/gg6.rds")

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
saveRDS(gg7, "cs13/data/gg7.rds")

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
saveRDS(gg8, "cs13/data/gg8.rds")

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
saveRDS(gg9, "cs13/data/gg9.rds")

## S
S <- bind_rows(tibble::as_tibble(jags_S) %>%
                 mutate(draw = row_number(), model = "jags") %>%
                 tidyr::pivot_longer(cols = -c(draw, model)),
               tibble::as_tibble(stan_S) %>% 
                 purrr::set_names(paste0("S[", 1:10, "]")) %>%
                 mutate(draw = row_number(), model = "stan") %>%
                 tidyr::pivot_longer(cols = -c(draw, model)))
true_S <- tibble(name = paste0("S[", 1:10, "]"),  true = sim_list$sample_effects_g[, 1])
S <- left_join(S, true_S) %>%
  mutate(sample_number = as.numeric(stringr::str_extract(name, "\\d+")))

ggplot(S, aes(value, fill = model)) +
  geom_density(alpha = 0.6) + 
  geom_vline(aes(xintercept = true)) + 
  facet_wrap(.~sample_number, ncol = 2) +
  theme_minimal() + 
  ggtitle("Posterior density for log normfactors")

