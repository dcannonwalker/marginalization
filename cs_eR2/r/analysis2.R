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
sim_list <- readRDS("cs_eR2/data/sim_list.rds")
G <- sim_list$G
stan_post <- readRDS("cs_eR2/data/stan_fit.rds")
eR_fit <- readRDS("cs_eR2/data/eR_fit.rds")
stan_smry <- stan_post$summary()

# true nulls
# reminder: D = 1 <=> non-null tag
tn <- 1 - sim_list$D_g

# 'p' ============
stan_p <- stan_smry %>%
  filter(grepl("^p", variable)) %>%
  select(mean, median, sd)

ps <- 1 - stan_p$mean

pe <- eR_fit$tt$table$PValue

pme <- calc_perf_metrics(pe, tn, model = "edger", fdr_mtd = "fdr")

pms <- calc_perf_metrics(ps, tn, model = "stan", fdr_mtd = "bfdr")

pm <- rbind(pme, pms)

roc_plot <- ggplot(pm, aes(fpr, tpr, color = model)) + 
  geom_line() +
  theme_minimal() 
  # facet_wrap(.~model)
fdr_plot <- ggplot(pm, aes(tfdr, nfdr, color = model)) + 
  geom_line() + 
  theme_minimal() +
  coord_fixed() + 
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  geom_abline()

# 'b' ===========
stan_b <- stan_smry %>%
  filter(grepl("b1\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(model = "stan", p = stan_p$mean, estimate = mean * p, mean = NULL, p = NULL) 
eR_b <- tibble::tibble(estimate = eR_fit$tt$table$logFC, variable = stan_b$variable, model = "edger")
b_df <- rbind(stan_b, eR_b) %>%
  mutate(
    true_class = rep(sim_list$D_g, 2),
    true = rep(sim_list$b1_g * sim_list$D_g, 2)) 

b1_plot <- ggplot(b_df, aes(true, estimate, color = model)) +
  geom_point(alpha = 0.5) +
  coord_fixed() +
  geom_abline() +
  theme_minimal()

saveRDS(b1_plot, file = "cs_eR2/data/a2_b1_plot.rds")
saveRDS(roc_plot, file = "cs_eR2/data/a2_roc_plot.rds")
saveRDS(fdr_plot, file = "cs_eR2/data/a2_fdr_plot.rds")
