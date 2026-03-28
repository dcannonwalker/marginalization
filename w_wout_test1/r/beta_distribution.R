prob <- 0.9
a <- 4
b <- (1 - prob) * a / prob

x <- seq(0, 1, length = 10000)
plot(x, dbeta(x, a, b))

mu <- 0.9
bss <- 50
a <- mu * bss
b <- bss - mu * bss
x <- seq(0, 1, length = 10000)
plot(x, dbeta(x, a, b))

dbeta(0.5, a, b)
dbeta(0.9, a, b)
