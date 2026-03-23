library(cmdstanr)
library(dplyr)
library(ggplot2)
mwout <- cmdstan_model("priors/wout_fixed_s.stan")
data(kidney, package = "SimSeq")
str(kidney)
keep <- edgeR::filterByExpr(kidney$counts, min.count = 100)
counts <- kidney$counts[keep, ]
G <- 1000
N <- 144
P <- N / 2
x <- as.numeric(kidney$treatment[1:N]) - 1
z <- model.matrix(~ 0 + factor(kidney$replic[1:N]))
po <- 0
y <- counts[1:G, 1:N]
s <- log(edgeR::normLibSizes(y) * colSums(y))

standata <- list(
  G = G,
  N = N,
  P = P,
  x = x,
  z = z,
  po = po,
  y = y,
  mao = 0,
  s = s
)

fit <- mwout$sample(data = standata, chains = 1)
fit$save_object("priors/wout_fixed_s_fit.rds")
smry <- fit$summary()

unames <- smry %>%
  filter(grepl("u\\[", variable)) %>%
  pull(variable)
unumbers <- stringr::str_extract_all(unames, "\\d+", simplify = TRUE) 
unumbers <- apply(unumbers, 2, as.numeric)

smry %>%
  filter(grepl("u\\[", variable)) %>%
  select(variable, mean) %>%
  mutate(tag = unumbers[, 1], pair = unumbers[, 2]) %>%
  ggplot(aes(mean, fill = factor(tag))) + 
  geom_density() +
  theme_minimal() +
  guides(fill = "none")
  
smry %>%
  filter(grepl("s\\[", variable)) %>%
  select(variable, mean) %>%
  ggplot(aes(mean)) + 
  geom_histogram()
smry %>%
  filter(variable == "ma")
smry %>%
  select(variable, rhat) %>%
  filter(variable %in% unames) %>%
  arrange(-rhat)


