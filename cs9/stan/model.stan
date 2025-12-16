//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_g;
  int<lower=0> N_u_g;
  int<lower=0> G;
  array[G, N_g] int<lower=0> y_g;
  vector[N_g] x_g;
  matrix[N_g, N_u_g] z_g;
  // real<lower=0> sig;
  real<lower=0, upper=1> pi0;
  real<lower=0> sig_u;
}

transformed data {
  vector[2] log_p_i;
  log_p_i[1] = log(pi0);
  log_p_i[2] = log(1 - pi0);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[G] b0;
  vector[G] b1; 
  array[G] vector[N_u_g] u_g;
  real<lower=0> sig;
  real<lower=0> sig_b0;
  real<lower=0> sig_b1;
}

transformed parameters {
  array[G] vector[2] lp;
  array[G, 2] vector[N_g] mu;
  vector[G] lse;
  vector[G] u_contr;
  for (i in 1:G) {
   u_contr[i] = normal_lpdf(u_g[i, ] | 0, sig_u);
   for (d in 1:2) {
     mu[i, d] = b0[i] + (d - 1) * x_g * b1[i] + z_g * u_g[i];
     lp[i, d] = log_p_i[d] + poisson_log_lpmf(y_g[i] | mu[i, d]);
   } 
   lse[i] = log_sum_exp(lp[i]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // sig, sig_b0, sig_b1 all have a uniform improper prior on the positive reals
  b0 ~ normal(0, sig_b0);
  b1 ~ normal(0, sig_b1);
  target += sum(u_contr);
  target += sum(lse);
}

generated quantities {
  array[G] real p; 
  for (i in 1:G) {
    p[i] = exp(lp[i, 2]) / sum(exp(lp[i]));
  }
}

