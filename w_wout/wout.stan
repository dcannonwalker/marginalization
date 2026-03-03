data {
  int<lower=0> N;
  int<lower=0> G;
  int<lower=0> P;
  vector[N] x;
  matrix[N, P] z;
  array[G, N] int<lower=0> y;
}

parameters {
  vector[G] alpha;
  vector[G] beta;
  vector<lower=0>[G] phi;
  real ma;
  real mb;
  real<lower=0> sa;
  real<lower=0> sb;
  real<lower=0> sp;
}

model {
  alpha ~ normal(ma, sa);
  beta ~ normal(mb, sb);
  phi ~ normal(0, sp);
}
