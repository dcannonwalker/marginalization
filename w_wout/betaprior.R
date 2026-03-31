x <- seq(0, 1, length = 10000)
pi0 <- 0.8 
a <- 4
b <- (1 - pi0) * a / pi0
plot(x, dbeta(x, a, b))
