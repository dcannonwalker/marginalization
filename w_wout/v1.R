# first test comparison

# setup
library(posterior)

mwout <- cmdstanr::cmdstan_model("w_wout/wout.stan")
mw <- cmdstanr::cmdstan_model("w_wout/w.stan")
dpo <- list(
  N = 10, 
  G = 1000,
  P = 0,
  x = rep(c(0, 1), each = 5),
  # z = model.matrix(~0 + factor(rep(1:5, 2))),
  z = numeric(0),
  po = 1,
  y = numeric(0),
  mao = 6 
)
dpow <- dpo
dpowout <- dpo
dpow$prob <- 0.8 

# simulation from w
pow <- mw$sample(dpow, chains = 1, iter_sampling = 500, iter_warmup = 1000)
drpow <- pow$draws()
ysimw <- extract_variable_array(drpow, "ysim")
bsimw <- extract_variable_array(drpow, "beta")
asimw <- extract_variable_array(drpow, "alpha")
usimw <-  extract_variable_array(drpow, "u")
Dsimw <- extract_variable_array(drpow, "D")
pi0simw <-extract_variable_array(drpow, "pi0") 
bsimw[1, 1, ]
median(pi0simw)
ysimw1 <- ysimw[1, 1, , ]
sum(Dsimw[1, 1, ])

# simulation from wout
powout <- mwout$sample(dpowout, chains = 1, iter_sampling = 500, iter_warmup = 1000)
drpowout <- powout$draws()
ysimwout <- extract_variable_array(drpowout, "ysim")
ysimwout1 <- ysimwout[1, 1, , ]

# fits to the w simulation
dlw <- dpow
dlw$po <- 0
dlw$y <- ysimw1
dlwout <- dpowout
dlwout$po <- 0
dlwout$y <- ysimw1

wfw <- mw$sample(dlw, chains = 1)
wfw$save_object("w_wout/data/wfw.rds")
wfwout <- mwout$sample(dlwout, chains = 1)
wfwout$save_object("w_wout/data/wfwout.rds")

# fits to the w simulation
dlw <- dpow
dlw$po <- 0
dlw$y <- ysimwout1
dlwout <- dpowout
dlwout$po <- 0
dlwout$y <- ysimwout1

woutfw <- mw$sample(dlw, chains = 1)
woutfw$save_object("w_wout/data/woutfw.rds")
woutfwout <- mwout$sample(dlwout, chains = 1)
woutfwout$save_object("w_wout/data/woutfwout.rds")

pw <- extract_variable_array(wfw, "p")
pinvw <- extract_variable_array(wfw, "pinv")
betaww <- extract_variable_array(wfw, "beta")
sbww <- extract_variable_array(wfw, "sb")
hist(betaww[, , 2])
alphaww <- extract_variable_array(wfw, "alpha")
str(betaww)
plot(apply(betaww, 3, mean), bsimw[1, 1, ])
plot(apply(alphaww, 3, mean), asimw[1, 1, ])
draw_hist <- function(var_arr, id) {
  hist(var_arr[, 1, id])
}
pinvw
par(mfrow = c(1, 1))
for (i in 1:16) {
  draw_hist(pinvw, i)
}
for (i in 81:96) {
  draw_hist(pw, i)
}
  
pwout <- extract_variable_array(wfw, "p")
pwout
ll <- extract_variable_array(wfw, "ll")
str(pw)

pscr <- apply(ll, 1:3, function(l) exp(l[2]) / sum(exp(l)))
pinv <- apply(ll, 1:3, function(l) 1 / (1 + exp(l[1] - l[2])))
str(pscr)
hist(pscr[, 1, 7])
pscr[, 1, 1]
pinv[, 1, 1]
