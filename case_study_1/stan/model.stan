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
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real<lower=0> sig;
  real<lower=0, upper=1> pi0;
}

transformed data {
  vector[2] log_p_i;
  log_p_i[1] = log(pi0);
  log_p_i[2] = log(1 - pi0);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real b0;
  real b1; 
}

transformed parameters {
  array[N] vector[2] lp;
  array[N] vector[2] mu;
  vector[N] lse;
  for (i in 1:N) {
   for (d in 1:2) {
     mu[i, d] = b0 + (d - 1) * x[i] * b1;
     lp[i, d] = log_p_i[d] + normal_lpdf(y[i] | mu[i, d], sig);
   } 
   lse[i] = log_sum_exp(lp[i]);
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  b0 ~ normal(0, 1);
  b1 ~ normal(0, 1);
  target += sum(lse);
}

generated quantities {
  array[N] real p; 
  for (i in 1:N) {
    p[i] = exp(lp[i, 2]) / sum(exp(lp[i]));
  }
}

