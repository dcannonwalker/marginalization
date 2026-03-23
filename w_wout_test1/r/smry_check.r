# setup ----
library(dplyr)
library(ggplot2)
library(tidyr)
smryw <- readRDS("w_wout_test1/data/smryw.rds")
smrywout <- readRDS("w_wout_test1/data/smrywout.rds")

# rhat checks ----
smryw %>%
  arrange(-rhat)

smrywout %>%
  arrange(-rhat)

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

# beta vs. p for w model ----
smryw %>%
  filter(variable %in% c(pnames, betanames)) %>%
  mutate(vartype = stringr::str_split_i(variable, "\\[", 1), 
         tag = stringr::str_extract(variable, "\\d+")) %>%
  select(tag, vartype, mean) %>%
  pivot_wider(names_from = vartype, values_from = mean) %>%
  ggplot(aes(p, beta)) +
  geom_point() + theme_minimal()
