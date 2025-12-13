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
  int<lower=0> G;
  array[G, N_g] real y_g;
  vector[N_g] x_g;
  real<lower=0> sig;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[G] b0;
  vector[G] b1; 
}

transformed parameters {
  array[G] real ll;
  array[G] vector[N_g] mu;
  for (i in 1:G) {
    mu[i] = b0[i] + x_g * b1[i];
    ll[i] = normal_lpdf(y_g[i] | mu[i], sig);
  } 
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  b0 ~ normal(0, sqrt(10));
  b1 ~ normal(0, sqrt(10));
  target += sum(ll);
}
