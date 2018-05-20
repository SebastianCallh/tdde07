data {
  int<lower=1> N;
  real x[N];
}
parameters {
  real mu;
  real<lower=0> sigma2;
  real<lower=-1,upper=1> phi;
}
transformed parameters {
  real sigma;
  sigma = sqrt(sigma2);
}
model {
 for (n in 2:N)
   x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
}
