data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> P;
  vector[N] x;
  matrix[N, P] z;
  // array[G, N] int<lower=0> y;
  int<lower=0, upper=1> po;
  array[po ? 0:G, N] int<lower=0> y;
  real mao; 
}

parameters {
  array[G] real alpha;
  array[G] real beta;
  array[G] real<lower=0> bcv;
  array[G] vector[P] u;
  real ma;
  real mb;
  real<lower=0> sa;
  real<lower=0> sb;
  real<lower=0> sp;
  real<lower=0> su;
}

transformed parameters {
  array[G] real<lower=0> phi;
  array[G] vector[N] eta;
  array[G] real u_contr;
  array[G] real ll;
  for (i in 1:G) {
    phi[i] = 1 / bcv[i]^2;
    u_contr[i] = normal_lpdf(u[i] | 0, su);
    eta[i] = alpha[i] + x * beta[i] + z * u[i];
    if (po == 0) {
      ll[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi[i]);
    }
  }
}

model {
  alpha ~ normal(ma, sa);
  beta ~ normal(mb, sb);
  bcv ~ normal(0, sp);
  ma ~ normal(mao, 1);
  mb ~ normal(0, 1); 
  sa ~ normal(0, 1);
  sb ~ normal(0, 1); 
  sp ~ normal(0, 1); 
  su ~ normal(0, 1);
  if (po == 0) target += sum(ll);
  target += sum(u_contr);
}

generated quantities {
  array[G, N] int<lower=0> ysim;
  for (i in 1:G) {
    ysim[i] = neg_binomial_2_log_rng(eta[i], phi[i]);
  }
}
