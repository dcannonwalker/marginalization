library(ggplot2)
sim_list <- readRDS("cs4_1/data/sim_list.rds")
sim_list$b1_g

b1 <- seq(-10, 10, by = 0.01)
dens <- 0.5 * dnorm(b1, -4, 1) + 0.5 * dnorm(b1, 4, 1)
df <- tibble::tibble(b1 = b1, dens = dens)
b1_plot <- ggplot(df, aes(b1, dens)) + 
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  ggtitle("Simulation distribution for b1") + 
  ylab("Density")

saveRDS(b1_plot, "cs4_1/data/sim_b1.rds")
