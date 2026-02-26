# ======== de.helpers =========
# functions from de.helpers -- consider installing the package to the project
# later
calc_perf_metrics <- function(p, tn, model, presorted = FALSE, nfdr = NULL, fdr_mtd = c("bfdr", "fdr")) {
  if (!presorted) {
    ord <- order(p)
    p <- p[ord]
    tn <- tn[ord]
  }
  fpr <- sapply(p, calc_rate, p = p, s = tn)
  tpr <- sapply(p, calc_rate, p = p, s = 1 - tn)
  tfdr <- sapply(p, calc_fdr, p = p, tn = tn)
  
  if (is.null(nfdr)) {
    fdr_mtd <- match.arg(fdr_mtd)
    fdr_fn <- switch (fdr_mtd,
                      bfdr = calc_bfdr,
                      fdr = function(p) p.adjust(p, method = "fdr")
    )
    nfdr <- fdr_fn(p)
  }
  tibble::tibble(model = model, p = p, tn = tn, fpr = fpr, tpr = tpr, tfdr = tfdr, nfdr = nfdr)
}
calc_fdr <- function(p0, p, tn) {
  sum(tn[p <= p0]) / sum(p <= p0)
}
calc_rate <- function(p0, p, s) {
  x <- sum(s[p <= p0])
  x / sum(s)
}
#' Calculate the estimated false discovery rate
#' @param p A vector of posterior probabilities
#' @export
calc_bfdr <- function(p, presorted = FALSE) {
  if (!presorted) p <- p[order(p)]
  sapply(p, function(p0) {
    sum(p[p <= p0]) / sum(p <= p0)
  })
}

# setup
library(posterior)
library(dplyr)
library(ggplot2)
library(cmdstanr)

# data
sim_list <- readRDS("cs_BYN1/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs_BYN1/stan_fit_9.rds")
# stan_smry <- stan_post$summary()
stan_smry <- readRDS("cs_BYN1/stan_smry_9.rds")

nb_post <- readRDS("cs_BYN1/nb_fit_9.rds")
# nb_smry <- nb_post$summary()
nb_smry <- readRDS("cs_BYN1/nb_smry_9.rds")

eR_fit <- readRDS("cs_BYN1/erfit.rds")

# true nulls
# reminder: D = 1 <=> non-null tag
tn <- sim_list$ds$mu1 == sim_list$ds$mu2
sum(tn)

# 'p' ============
stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(mean, median, sd)

ps <- 1 - stan_p$mean

nb_p <- nb_smry %>%
  filter(grepl("^p\\[", variable)) %>%
  select(mean, median, sd)

pnb <- 1 - nb_p$mean

pe <- eR_fit$tt$table$PValue

pme <- calc_perf_metrics(pe, tn, model = "edger", fdr_mtd = "fdr")

pms <- calc_perf_metrics(ps, tn, model = "stan", fdr_mtd = "bfdr")

pmnb <- calc_perf_metrics(pnb, tn, model = "stan_nb", fdr_mtd = "bfdr")

pm <- rbind(pme, pms, pmnb)

roc_plot <- ggplot(pm, aes(fpr, tpr, color = model)) + 
  geom_line() +
  theme_minimal()  + 
  ggtitle("ROC curve")
# facet_wrap(.~model)
fdr_plot <- ggplot(pm, aes(nfdr, tfdr, color = model)) + 
  geom_line() + 
  theme_minimal() +
  coord_fixed() + 
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  geom_abline() +
  ggtitle("FDR control") + 
  xlab("True FDR") + 
  ylab("Nominal FDR")

# 'b' ===========
stan_b <- stan_smry %>%
  filter(grepl("^b1\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan", p = stan_p$mean, estimate = mean * p, mean = NULL, p = NULL) 

nb_b <- nb_smry %>%
  filter(grepl("^b1\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan_nb", p = nb_p$mean, estimate = mean * p, mean = NULL, p = NULL) 

eR_b <- tibble::tibble(estimate = eR_fit$tt$table$logFC, variable = stan_b$variable, model = "edger")

D_g <- 1 - tn
b1_g <- log(sim_list$ds$mu2 / sim_list$ds$mu1)

b_df <- rbind(stan_b, eR_b, nb_b) %>%
  mutate(
    true_class = rep(D_g, 3),
    true = rep(b1_g * D_g, 3)) 

b1_plot <- ggplot(b_df, aes(true, estimate, color = model)) +
  geom_point(alpha = 0.5) +
  # coord_fixed() +
  geom_abline() +
  theme_minimal() + 
  ggtitle("logFC estimates against true values")

# intercept =========
stan_b0 <- stan_smry %>%
  filter(grepl("^b0\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan", estimate = mean, mean = NULL) 

nb_b0 <- nb_smry %>%
  filter(grepl("^b0\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan_nb", estimate = mean, mean = NULL) 

# eR_b0 <- tibble::tibble(estimate = eR_fit$tt$table$, variable = stan_b0$variable, model = "edger")

D_g <- 1 - tn
b0_g <- log(sim_list$ds$mu1)

b0_df <- rbind(stan_b0, nb_b0) %>%
  mutate(
    true = rep(b0_g, 2)) 

b0_plot <- ggplot(b0_df, aes(true, estimate, color = model)) +
  geom_point(alpha = 0.5) +
  # coord_fixed() +
  geom_abline() +
  theme_minimal() + 
  ggtitle("Log-scale intercept estimates against true values")

# S =============
stan_S <- extract_variable_array(stan_post, variable = "S")[, 1, ]
nb_S <- extract_variable_array(nb_post, variable = "S")[, 1, ]
S <- bind_rows(tibble::as_tibble(stan_S) %>% 
                 purrr::set_names(paste0("S[", 1:10, "]")) %>%
                 mutate(draw = row_number(), model = "stan") %>%
                 tidyr::pivot_longer(cols = -c(draw, model)),
               tibble::as_tibble(nb_S) %>% 
                 purrr::set_names(paste0("S[", 1:10, "]")) %>%
                 mutate(draw = row_number(), model = "stan_nb") %>%
                 tidyr::pivot_longer(cols = -c(draw, model)))
true_S <- tibble(name = paste0("S[", 1:10, "]"),  true = log(c(sim_list$ds$S1, sim_list$ds$S2)))
S <- left_join(S, true_S) %>%
  mutate(sample_number = as.numeric(stringr::str_extract(name, "\\d+")))

S_plot <- ggplot(S %>% filter(sample_number < 3), aes(value, fill = model)) +
  geom_density(alpha = 0.6) + 
  geom_vline(aes(xintercept = true), color = "red") + 
  facet_wrap(.~sample_number) +
  theme_minimal() + 
  ggtitle("Posterior density for log normfactors - samples 1 & 2", subtitle = "Posteriors for the other eight look very similar") 
# 'phi' ==========
nb_phi <- nb_smry %>%
  filter(grepl("^phi\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan", estimate = mean, mean = NULL) 

er_phi <- tibble::tibble(estimate = 1 / eR_fit$fit$dispersion, variable = nb_phi$variable, model = "edger")

true_phi <- 1 / sim_list$ds$phi1
phi_df <- rbind(nb_phi, er_phi) %>%
  mutate(
    true_class = rep(D_g, 2),
    true = rep(true_phi, 2)) 

phi_plot <- ggplot(phi_df, aes(true, estimate, color = model)) +
  geom_point(alpha = 0.5) +
  # coord_fixed() +
  geom_abline() +
  theme_minimal() + 
  ggtitle("Dispersion estimates against true values")

saveRDS(b1_plot, file = "cs_BYN1/b1_plot_9.rds")
saveRDS(b0_plot, file = "cs_BYN1/b0_plot_9.rds")
saveRDS(S_plot, file = "cs_BYN1/S_plot_9.rds")
saveRDS(roc_plot, file = "cs_BYN1/roc_plot_9.rds")
saveRDS(fdr_plot, file = "cs_BYN1/fdr_plot_9.rds")
saveRDS(phi_plot, file = "cs_BYN1/phi_plot_9.rds")

