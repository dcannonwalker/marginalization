data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> P;
  vector[N] x;
  matrix[N, P] z;
  int<lower=0, upper=1> po; // prior only?
  array[po ? 0:G, N] int<lower=0> y; 
  real mao; // positive offset for ma, the mean of alpha
  vector[N] s; // log effective library size 
}

parameters {
  array[G] real alpha; // intercept
  array[G] real beta; // trt effect
  array[G] real<lower=0> bcv; // 1 / sqrt(size) for nb
  array[G] vector[P] u; // ranefs
  real ma; // mean for alpha
  real mb; // mean for beta
  // vector[N] s; // sample effects
  real<lower=0> sa; // scale for alpha 
  real<lower=0> sb; // scale for beta
  real<lower=0> sp; // scale for bcv (named for phi)
  real<lower=0> su; // scale for u
  // real<lower=0> ss; // scale for s
}

transformed parameters {
  array[G] real<lower=0> phi;
  array[G] vector[N] eta;
  array[G] real u_contr;
  array[G] real ll;
  for (i in 1:G) {
    phi[i] = 1 / bcv[i]^2;
    u_contr[i] = normal_lpdf(u[i] | 0, su);
    eta[i] = alpha[i] + x * beta[i] + z * u[i] + s;
    if (po == 0) {
      ll[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi[i]);
    }
  }
}

model {
  alpha ~ normal(ma, sa);
  beta ~ normal(mb, sb);
  bcv ~ normal(0, sp);
  // s ~ normal(0, ss);
  ma ~ normal(mao, 1);
  mb ~ normal(0, 1); 
  sa ~ normal(0, 1);
  sb ~ normal(0, 1); 
  sp ~ normal(0, 1); 
  su ~ normal(0, 1);
  // ss ~ normal(0, 1);
  if (po == 0) target += sum(ll);
  target += sum(u_contr);
}

generated quantities {
  array[G, N] int<lower=0> ysim;
  for (i in 1:G) {
    ysim[i] = neg_binomial_2_log_rng(eta[i], phi[i]);
  }
}
