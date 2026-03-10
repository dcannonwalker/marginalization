data {
  int<lower=0> N;
  int<lower=0> G;
  array[G, N] int<lower=0> y;
  vector[N] x;
  real<lower=0> sig_S;
  matrix[N, N] sample_design;
}

parameters {
  array[G] alpha; // log-scale intercept
  array[G] beta; // log-scale treatment effect
  array[N] S; // log-scale sample effect
  array[G] real<lower=0> bcv; // 1 / sqrt(phi), where phi is the inverse overdispersion,
  // i.e. var(y) = mu + mu^2 / phi
  // named for the Biological Coef. of Variation
}

transformed parameters {
  array[G, N] eta; 
  array[G] ll;
  array[G] phi;
  for (i in 1:G) {
    phi[i] = (1 / bcv[i])^2
    eta[i] = alpha[i] + x * beta[i] + S;
    ll[i] = neg_binomial_2_log_lpmf(y[i] | eta[i], phi[i]);
  }
}

model {
  alpha ~ normal(0, 1);
  beta ~ normal(0, 1); 
  S ~ normal(0, sig_S);
  bcv ~ normal(0, 1);
}
