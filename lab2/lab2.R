dinvgamma <- function(x, a, b) {
  (b^a)/gamma(a) * x^(-a-1) * exp(-b/x)
}

dinvchisq2 <- function(x, n, tau2) {
  a <- n / 2
  b <- n*tau2 / 2
  dinvgamma(x, a, b)
}

sigma.prior <- function (x) {
  v0     <- 20
  sigma0 <- 2
  sigma.delta <- 0.05
  dinvchisq2(x, v0, sigma0)
}

theta.prior <- function(mu, sigma) {
  theta.delta <- 0.05
  dnorm(x, mu, sigma)
}

sigma.grid  <- seq(0, 20, sigma.delta)
plot(sigma.grid, sapply(sigma.grid, sigma.prior),
     main = "Sigma^2 prior",
     xlab = expression(sigma^2),
     ylab = "Density",
     type = 'l')  

theta.grid  <- seq(-5, 5, theta.delta)
plot(theta.grid, sapply(theta.grid, theta.prior),
     main = "Theta prior",
     xlab = expression(theta),
     ylab = "Density",
     type = 'l')

nDraws  <- 5000 
v0      <- 20
sigma0  <- 2
sigma2s <- v0*sigma0/rchisq(nDraws, v0)
plot(density(sigma2))

mu0    <- 8 # celcius from TDDE01 
kappa0 <- 2 # celcius from TDDE01 
thetas <- sapply(sigma2s, function(sigma2) { rnorm(1, mu0, sqrt(sigma2/kappa0)) })
plot(density(thetas))
mean(thetas)
