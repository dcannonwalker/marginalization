# setup
library(posterior)
library(dplyr)
library(ggplot2)
library(cmdstanr)

# data
sim_list <- readRDS("cs_eR4/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs_eR4/data/stan_fit.rds")
stan_smry <- readRDS("cs_eR4/data/stan_smry.rds")

# true nulls
# reminder: D = 1 <=> non-null tag
tn <- 1 - sim_list$D_g
tb1 <- sim_list$b1_g * sim_list$D_g
tb0 <- sim_list$b0_g
tu <- sim_list$u_g

# 'p' ============
stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(variable, mean)

# 'b0' ===========
stan_b0 <- stan_smry %>%
  filter(grepl("^b0", variable)) %>%
  select(variable, mean)
stan_u <- stan_smry %>%
  filter(grepl("^u_g", variable)) %>%
  select(variable, mean)

u_id <- stringr::str_extract_all(stan_u$variable, pattern = "\\d+", simplify = TRUE)
u_id <- apply(u_id, 2, as.numeric)

u_df <- stan_u %>%
  mutate(tag = u_id[, 1], pair = u_id[, 2]) %>%
  tidyr::pivot_wider(id_cols = tag, names_from = pair, values_from = mean, names_prefix = "u_pair_")

tu_df <- as_tibble(tu)
colnames(tu_df) <- paste0("tu_pair_", 1:5)
tu_df <- tu_df %>%
  mutate(tag = row_number())

# 'b1' ===========
x <- stan_smry %>%
  filter(grepl("^b1\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(tag = row_number(), model = "stan", p = stan_p$mean, tn = tn, b1 = mean * p, tb1 = tb1, b0 = stan_b0$mean, tb0 = tb0, mean = NULL, variable = NULL) %>%
  left_join(u_df) %>%
  left_join(tu_df)
x %>%
  arrange(-p) %>%
  print(n = 75)
