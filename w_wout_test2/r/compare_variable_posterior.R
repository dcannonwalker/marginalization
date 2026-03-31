f1 <- readRDS("w_wout_test1/data/fitw.rds")
f2 <- readRDS("w_wout_test1/data/fitwout.rds")
var <- "beta"
tag <- 1
library(posterior)
library(dplyr)
library(ggplot2)
compare_variable_posterior <- function(f1, f2, var, tag, mnames = c("w", "wout")) {
  d1 <- merge_chains(f1$draws())
  d2 <- merge_chains(f2$draws())
  var1 <- extract_variable_array(d1, var) 
  var2 <- extract_variable_array(d2, var)
  arrdim <- dim(var1) # assume the third dim is tags
  if (arrdim[3] > 1) {
    if (length(arrdim) == 3) {
      var1 <- var1[, , tag]
      var2 <- var2[, , tag]
    } else if (length(arrdim) == 4) {
      var1 <- var1[, , tag, ]
      var2 <- var2[, , tag, ]
    }
  }
  out <- cbind(var1, var2)
  colnames(out) <- paste0(var, "_", mnames)
  out
}
varmat <- compare_variable_posterior(f1 = f1, f2 = f2, var = "pi0", tag = 5)
as_tibble(varmat) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(cols = -draw) %>%
  ggplot(aes(value, color = factor(name))) + 
  geom_freqpoly() + 
  theme_minimal() 

