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
  real<lower=0, upper=1> prob; // prior mean proportion null
  vector[N] s; // log effective library size 
}

transformed data {
  real a = 4; 
  real b = (1 - prob) * a / prob;
}

parameters {
  array[G] real alpha;
  array[G] real beta;
  array[G] real<lower=0> bcv;
  array[G] vector[P] u;
  real ma;
  real mb;
  real<lower=0, upper=1> pi0; // overall proportion of nulls
  real<lower=0> sa;
  real<lower=0> sb;
  real<lower=0> sp;
  real<lower=0> su;
}

transformed parameters {
  array[G] real<lower=0> phi;
  array[G, 2] vector[N] eta; //eta[i, 1] log mean for null status
  array[G] real u_contr;
  vector[2] log_pi0;
  log_pi0[1] = log(pi0); // log overall proportion nulls
  log_pi0[2] = log(1 - pi0); // log overall proportion non-nulls
  array[G] vector[2] ll; // ll[i, 1]] log likelihood for null status
  vector[G] lse;
  for (i in 1:G) {
    phi[i] = 1 / bcv[i]^2;
    u_contr[i] = normal_lpdf(u[i] | 0, su);
    for (d in 1:2) { // d = 2 => null
      eta[i, d] = alpha[i] + (2 - d) * x * beta[i] + z * u[i] + s;
      if (po == 0) {
        ll[i, d] = log_pi0[d] + neg_binomial_2_log_lpmf(y[i] | eta[i, d], phi[i]);
      }
    } 
    lse[i] = log_sum_exp(ll[i]);
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
  pi0 ~ beta(a, b);
  if (po == 0) target += sum(lse);
  target += sum(u_contr);
}

generated quantities {
  array[G, N] int<lower=0> ysim;
  array[G] real p; 
  array[G] real pinv;
  array[G] int<lower=0, upper=1> D;
  array[G] real bmix;
  for (i in 1:G) {
    if (po == 0) {
      p[i] = exp(ll[i, 2]) / sum(exp(ll[i]));
      pinv[i] = 1 / (1 + exp(ll[i, 1] - ll[i, 2]));
    }
    else p[i] = pi0;
    bmix[i] = (1 - p[i]) * beta[i];
    D[i] = bernoulli_rng(pi0); // is gene i null? 
    ysim[i] = neg_binomial_2_log_rng(eta[i, 1 + D[i]], phi[i]); 
  }
}
