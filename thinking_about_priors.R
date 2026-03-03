library(MASS)
library(edgeR)
library(dplyr) 
library(tidyr)
library(ggplot2)
data(kidney, package = 'SimSeq')
trt <- factor(kidney$treatment)
pair <- factor(kidney$replic)
design <- model.matrix(~ pair + trt)
design0 <- model.matrix(~ trt)

samples <- data.frame(
  treatment = trt,
  replic = pair
)

kirc <- DGEList(counts = kidney$counts, samples = samples, group = trt)

kirc <- normLibSizes(kirc)

r2u <- filterByExpr(kirc)
kirc <- kirc[r2u, ]

kirc <- estimateDisp(kirc)



kirc_glm_fit <- glmFit(kirc, design = design)
kirc_glm_fit0 <- glmFit(kirc, design = design0)

# comparing intercept estimates w/ and w/out pair effect
cf <- kirc_glm_fit$coefficients
cf0 <- kirc_glm_fit0$coefficients

plot(cf[, 1], cf0[, 1], asp = 1)
abline(a = 0, b = 1, col = "red")

# look at some histograms that approximate the vs
hist(cf[, 1] - cf0[, 1], breaks = 1e3)
hist(cf[, 1] + cf[, 5] - cf0[, 1], breaks = 1e3)

# look at hist to approx b0
hist(cf[, 1], breaks = 1000)

# look at hist to approx b1
hist(cf[, ncol(cf)], breaks = 1000)
hist(cf0[, ncol(cf0)], breaks = 1000)
plot(cf[, ncol(cf)], cf0[, ncol(cf0)], asp = 1)
abline(a = 0, b = 1)

# look at hist to approx 1 / sqrt(phi)
hist(kirc$tagwise.dispersion)

b0_comp <- tibble(b01 = cf[, 1], b00 = cf0[, 1])
b0_comp %>%
  mutate(id = row_number(), v1 = b01 - b00) %>%
  arrange(v1)




####### extra stuff
cf <- kirc_glm_fit$coefficients
ucf <- kirc_glm_fit$unshrunk.coefficients
str(cf)
str(ucf)
plot(cf[, 1], ucf[, 1], asp = 1)
abline(a = 0, b = 1, col = "red")

plot_cfs <- function(i) {
  plot(cf[, i], ucf[, i], asp = 1)
  abline(a = 0, b = 1, col = "red") 
}

plot_cfs(2)
plot_cfs(4)
hist(cf[, 2], breaks = 1000)

plotBCV(kirc)

hist(1 / sqrt(kirc$tagwise.dispersion))