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
  int<lower=0, upper=1> two_comp_mu_b1;
  array[G, N_g] int<lower=0> y_g;
  vector[N_g] x_g;
  matrix[N_g, N_u_g] z_g;
  // real<lower=0> sig;
  real<lower=0, upper=1> pi0;
  real<lower=0> sig_u;
  real<lower=0> sig_S;
  matrix[N_g, N_g] sample_design_g;
  vector<lower=0>[2] phi_bounds;
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
  vector[N_g] S;
  array[G] vector[N_u_g] u_g;
  real<lower=0> sig;
  vector[two_comp_mu_b1 + 1] mu_b1;
  real mu_b0;
  real<lower=0> sig_b0;
  real<lower=0> sig_b1;
  vector<lower=phi_bounds[1], upper=phi_bounds[2]>[G] phi;
}

transformed parameters {
  array[G] vector[2] lp;
  array[G, 2] vector[N_g] mu;
  vector[G] lse;
  vector[G] u_contr;
  vector[G] b1_contr;
  for (i in 1:G) {
   u_contr[i] = normal_lpdf(u_g[i, ] | 0, sig_u);
   if (two_comp_mu_b1 == 1) {
     b1_contr[i] = log_sum_exp(normal_lpdf(b1[i] | mu_b1[1], sig_b1),
                             normal_lpdf(b1[i] | mu_b1[2], sig_b1));
   } 
   else {
     b1_contr[i] = normal_lpdf(b1[i] | mu_b1[1], sig_b1);
   }
   for (d in 1:2) {
     mu[i, d] = b0[i] + (d - 1) * x_g * b1[i] + z_g * u_g[i] + sample_design_g * S;
     lp[i, d] = log_p_i[d] + neg_binomial_2_log_lpmf(y_g[i] | mu[i, d], phi[i]);
   } 
   lse[i] = log_sum_exp(lp[i]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // sig, sig_b0, sig_b1 all have a uniform improper prior on the positive reals
  b0 ~ normal(mu_b0, sig_b0);
  S ~ normal(0, sig_S);
  mu_b1 ~ normal(0, sqrt(10));
  mu_b0 ~ normal(0, sqrt(10));
  target += sum(b1_contr);
  target += sum(u_contr);
  target += sum(lse);
}

generated quantities {
  array[G] real p; 
  for (i in 1:G) {
    p[i] = exp(lp[i, 2]) / sum(exp(lp[i]));
  }
}

