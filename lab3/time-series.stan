data {
  int<lower=0> N;
  vector[N] x;
  real mu;
  real phi;
  real<lower=0> sigma;
}

#parameters {
#}

model {
  x[2:N] ~ normal(mu + phi * x[1:(N - 1)], sigma);
}