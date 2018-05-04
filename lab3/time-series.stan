data {
  int<lower=0> N;
  vector[N] x;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real phi;
}

model {
  x[2:N] ~ normal(mu + phi * x[1:(N - 1)], sigma);
}