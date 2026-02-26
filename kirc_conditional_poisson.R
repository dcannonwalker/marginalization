library(MASS)
library(edgeR)
library(dplyr) 
library(tidyr)
library(ggplot2)
data(kidney, package = 'SimSeq')
trt <- factor(kidney$treatment)
pair <- factor(kidney$replic)

G <- 1000
N <- 72
counts <- kidney$counts
row2keep <- filterByExpr(counts, min.count = 100)
counts <- counts[row2keep, ]
y <- DGEList(counts = counts)
y <- normLibSizes(y)
offsets <- log(y$samples$norm.factors * y$samples$lib.size)
sS <- offsets[1:N]
st <- trt[1:N]
sp <- factor(pair[1:N])
sc <- counts[1:G, 1:N]

sf_nb <- lapply(1:G, function(g) {
  tryCatch(expr = {
    ff <- glm.nb(sc[g, ] ~ st + sp + offset(sS))
  list(fit = ff, aic = ff$aic)
  }, error = function(e) {
    message("Fit errored out")
    list(fit = NULL, aic = NA)
  })
})
names(sf_nb) <- paste0("gene", 1:G)

sf_po <- lapply(1:G, function(g) {
  ff <- glm(sc[g, ] ~ st + sp + offset(sS), family = poisson)
  list(fit = ff, aic = ff$aic)
})
names(sf_po) <- paste0("gene", 1:G)

augment_aic <- function(aic1, aic2, cnames) {
  probs <- t(apply(cbind(aic1, aic2), 1, function(aic) exp(-(aic - min(aic)) / 2)))
  colnames(probs) <- paste0("prob_", cnames)
  out <- cbind(aic1, aic2)
  colnames(out) <- cnames
  out <- as_tibble(cbind(out, probs)) 
}

aic_nb <- sapply(sf_nb, function(ff) ff$aic)
aic_po <- sapply(sf_po, function(ff) ff$aic)

aic <- augment_aic(aic_nb, aic_po, c("nb", "pois")) %>% 
  mutate(preferred_model = factor(if_else(prob_nb > prob_pois, "nb", "pois"))) %>%
  arrange(-nb)
summary(aic$preferred_model)

aic %>%
  mutate(gene = row_number()) %>%
  pivot_longer(cols = c(nb, pois), names_to = "fit", values_to = "aic") %>%
  ggplot(aes(aic, fill = fit)) + 
  geom_histogram() + facet_wrap(.~fit)

p1 <- aic %>%
  ggplot(aes(nb, pois)) + 
  geom_point() + 
  theme_minimal() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_transform(y = "log") +
  ggtitle("Comparison of AIC for NB and Poisson GLMs", subtitle = "KIRC kidney cancer data set") +
  xlab("AIC for NB") +
  ylab("AIC for Poisson (log-scale axis marks)")
p1
saveRDS(p1, "kirc_aic_plot.rds")
