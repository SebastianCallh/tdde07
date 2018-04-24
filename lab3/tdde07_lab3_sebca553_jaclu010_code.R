data <- read.table("rainfall.dat", header=FALSE)
mu0    <- 20 # 0.5 cm
tau0   <- 5  # Variance of expected rainfall
nu0    <- 3  # Variance of expected value of variance
sigma0 <- 5  # Expected value of variance

nIter  <- 5000
theta0 <- c(20, 2)
x      <- data[,1]
n      <- length(x)

rinvchisq <- function(vx, sigmax, draws) {
  vx*sigmax/rchisq(draws, vx)
}

Mu.Conditional.Posterior.Draw <- function(sigma2) {
  taun <- sqrt(1/(n/sigma2 + 1/tau0^2))
  w    <- (n/sigma2) / (n/sigma2 + 1/tau0^2)
  mun  <- w*mean(x) + (1 - w) * mu0
  rnorm(1, mun, taun)
}

Sigma.Conditional.Posterior.Draw <- function(mu) {
  nun    <- n + nu0
  sigman <- (nu0*sigma0^2 + sum((x-mu)^2)) / nun
  rinvchisq(nun, sigman, 1)
}

Gibbs <- function(theta_t) {
  mu     <- Mu.Conditional.Posterior.Draw(theta_t[2])
  sigma2 <- Sigma.Conditional.Posterior.Draw(mu)
  c(mu, sigma2)
}

thetas     <- matrix(rep(0, nIter*2), nrow = nIter)
thetas[1,] <- theta0

for(i in 2:nIter) {
  thetas[i,] <- Gibbs(thetas[i-1,])
}

plot(thetas[,2], type='l')

#Reduce(Gibbs, rep(2, 10), theta0, accumulate = TRUE)


