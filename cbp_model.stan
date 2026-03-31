#include custom_beta_prior.stan
data {
  int<lower=0> N;
  int<lower=0> G;
  vector[N] x; // trt group labels (binary vector)
  real<lower=0> c; // log fold change cutoff
  real<lower=0> sb; // scale parameter for b
  real<lower=0> sa; // scale parameter for a
  real ma; // location parameter for a
  real<lower=0, upper=1> p0; // prior proportion null
  vector[N] s; // log-scale sample normalization factors
  array[G, N] int<lower=0> y;
}
parameters {
  array[G] real b;
  array[G] real a;
  real<lower=0> bcv;
} 
transformed parameters {
  array[G] vector[2] ll;
  vector[G] b_contr;
  for (g in 1:G) {
    b_contr[g] = cbp_lpdf(b[g] | sb, c);
    for (d in 0:1) {
      vector[N] mutemp = a[g] + x * d * b[g] + s;
      ll[g, d + 1] = p0 * (1 - d) + (1 - p0) + neg_binomial_2_lpmf(y[g] | mutemp, 1 / bcv^2);
    }
  }
}
model {
  target += sum(b_contr);
  bcv ~ normal(0, 1);
  a ~ normal(ma, sa);
}
generated quantities {
  array[G] real<lower=0, upper=1> p; // posterior probabilities of null expression
  for (g in 1:G) {
    p[g] = exp(ll[g, 1]) / sum(exp(ll[g]));
  }
}
