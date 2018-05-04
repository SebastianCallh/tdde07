data {
  int<lower=0> N;
  int c[N];
}

parameters {
  real mu;
  real<lower=0> sigma;
  real phi;
  vector[N] x;
}

model {
  // sigma  ~ normal(0, 0.02);
  phi    ~ normal(0, 0.6);
  x[2:N] ~ normal(mu + phi * x[1:(N - 1)], sigma);
  for (n in 1:N)
    c[n] ~ poisson(exp(x[n]));
}