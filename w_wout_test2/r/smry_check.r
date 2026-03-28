# setup ----
library(dplyr)
library(ggplot2)
library(tidyr)
source("w_wout_test2/r/calc_perf_metrics.r")
erfit <- readRDS("w_wout_test2/data/erfit.rds")
tt <- erfit$tt
smryw <- readRDS("w_wout_test2/data/smryw.rds")
smrywout <- readRDS("w_wout_test2/data/smrywout.rds")
simlist <- readRDS("w_wout_test2/data/sim_list.rds")
names(simlist)
beta_true <- -log(simlist$ds$mu1 / simlist$ds$mu2)
hist(beta_true, breaks = 20)

# rhat checks ----
smryw %>%
  arrange(-rhat) %>%
  # filter(grepl("bmix", variable)) %>%
  print(n = 100)

smrywout %>%
  arrange(-rhat) %>%
  print(n = 100)

combo <- bind_rows(list(
  w = smryw,
  wout = smrywout
), .id = "model")

# var name vecs ----
betanames <- paste0("beta[", 1:1000, "]")
bmixnames <- paste0("bmix[", 1:1000, "]")
allnames  <- c(betanames, bmixnames)
alphanames <- paste0("alpha[", 1:1000, "]")
pnames <- paste0("p[", 1:1000, "]")

# bmix vs betawout and betaw vs betawout ----
combo %>%
  filter(variable %in% allnames) %>%
  select(model, variable, mean) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  mutate(id = factor(paste0(vartype, model))) %>%
  select(tag, id, mean) %>%
  pivot_wider(names_from = id, values_from = mean) %>%
  ggplot(aes(betaw, betawout)) + 
  geom_point(alpha = 0.5, color = 'blue') +
  geom_point(aes(bmixw, betawout), alpha = 0.5, color = 'red') + 
  theme_minimal() + 
  geom_abline(intercept = 0, slope = 1, color = 'purple', linetype = 2)

# alpha plot ----
combo %>%
  filter(variable %in% alphanames) %>%
  arrange(variable) %>% 
  select(model, variable, mean) %>%
  pivot_wider(names_from = model, values_from = mean) %>%
  ggplot(aes(w, wout)) + 
  geom_point() + theme_minimal() + 
  geom_abline(intercept = 0, slope = 1)

# edger alpha plot
lrt <- erfit$lrt
smryw %>%
  filter(variable %in% alphanames) %>%
  select(variable, mean) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  mutate(id = factor(paste0(vartype))) %>%
  select(tag, id, mean) %>%
  mutate(edger_est = lrt$coefficients[, 1]) %>%
  ggplot(aes(mean, edger_est)) + 
  geom_point() + 
  coord_fixed() + 
  geom_abline()

# beta vs. p for w model ----
smryw %>%
  filter(variable %in% c(pnames, betanames)) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  select(tag, vartype, mean) %>%
  pivot_wider(names_from = vartype, values_from = mean) %>%
  ggplot(aes(p, beta)) +
  geom_point() + theme_minimal()

# comparisons to true values
combo %>%
  filter(variable %in% bmixnames) %>%
  select(variable, mean) %>%
  mutate(true = beta_true) %>%
  ggplot(aes(true, mean)) + 
  geom_point() + 
  geom_abline(intecept = 0, slope = 1)

smryw %>%
  filter(variable %in% c(pnames, betanames)) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  select(tag, vartype, mean) %>%
  pivot_wider(names_from = vartype, values_from = mean) %>%
  mutate(true = beta_true) %>%
  ggplot(aes(true, beta, size = p)) + 
  geom_point(alpha = 0.1)


qdf <- smryw %>%
  filter(variable %in% c(pnames, betanames)) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  select(tag, vartype, mean) %>%
  pivot_wider(names_from = vartype, values_from = mean) 

ngs_perf_metrics <- calc_perf_metrics(qdf$p, beta_true == 0, model = "stan_nb", fdr_mtd = "bfdr")
er_perf_metrics <- calc_perf_metrics(tt$table$PValue, beta_true == 0, model = "edger", fdr_mtd = "fdr")
rbind(ngs_perf_metrics, er_perf_metrics) %>%
  ggplot(aes(nfdr, tfdr, color = model)) + 
  geom_line() + 
  coord_fixed()

qdf %>%
  mutate(true = beta_true,
         er_beta = tt$table$logFC, 
         er_fdr = tt$table$FDR) %>%
  arrange(p) %>%
  mutate(ngs_fdr = ngs_perf_metrics$nfdr)

