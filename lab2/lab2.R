library(mvtnorm)
data <- read.csv("TempLinkoping.txt", sep = "\t")
X <- data$temp

# 1. Linear and polynomial regression

#a) Determining the prior distribution
delta <- 0.01
times <- seq(delta, 1, delta)
nDraws  <- 5000 
mu0    <- c(-10, 150, -150)   # From tuning the plotted polynomial
omega0 <- 1*diag(3)           # same variance for betas, we dont know much
v0      <- 20
sigma0  <- 2
sigma2s <- v0*sigma0/rchisq(nDraws, v0)
n       <- length(X)

# Plot of prior variance
plot(density(sigma2s),
     main = "Sigma^2 prior",
     xlab = expression(sigma^2),
     ylab = "Density",
     type = 'l')  

# Thetas are matrix of theta - all parameters with priors
thetas <- cbind(t(sapply(sigma2s, function(sigma2) { 
  rmvnorm(1, mu0, sigma2*solve(omega0))
})), sigma2s)

# y is the prior distribution of the model
y <- apply(thetas[seq(1, 10),], 1, function(theta) {
  sapply(times, function(time) {
    theta[1] + theta[2]*time + theta[3]*time^2 + rnorm(1, 0, theta[4])
  })
})


#b) Check if prior distribution is sensible
plot(y[,1], type = 'l', ylim = c(-15, 40))
for (i in 2:10) {
  lines(y[,i])
}

# c) Simulating from posterior distribution
A   <- as.numeric(t(X) %*% X)
beta_hat <- apply(thetas[, 1:3], 2, mean) # minimize sqaure loss
mun      <- solve(A + omega0) %*% (A*beta_hat + omega0*mu0)
omegan   <- A + omega0
vn       <- v0 + n
