data {
  int<lower=0> N;
  int c[N];
}

parameters {
  real mu;
  real<lower=0> sigma2;
  real phi;
  vector[N] x;
}
transformed parameters {
  real sigma;
  sigma = sqrt(sigma2);
}
model {
  // sigma  ~ normal(0, 0.02);
  // phi    ~ normal(0, 0.6);
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
  for (n in 1:N)
    c[n] ~ poisson(exp(x[n]));
}