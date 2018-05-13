
# 1 - Poisson regression - the MCMC way

# a)
ebay.data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
ebay.glm <- glm(nBids ~ .-Const, data = ebay.data, family = poisson(link=log))
summary(ebay.glm)

# Estimate tells us: MinBidShare most important -1.89, sealed 0.44, veriID 0.39, majBlem -0.22 ,logBook -0.12, bleh

# b)
library(optimr)
library(mvtnorm)

n <- dim(ebay.data)[1]
p <- dim(ebay.data)[2]
X <- scale(as.matrix(ebay.data[,3:10]))
Y <- ebay.data[,1]
sigma  <- 100*t(X)%*%X
log.poisson.posterior.of.beta <- function(beta) { 
  lambda <- exp(X %*% beta)
  beta.log.prior         <- dmvnorm(beta, sigma = sigma, log = TRUE)
  poisson.log.likelihood <- sum(dpois(Y, lambda, log = TRUE)) # log(lambda^Y*exp(-lambda))) #
  poisson.log.likelihood + beta.log.prior
}

beta.initial = rep(0, dim(X)[2])

## This gives values around 10^280
optim.beta <- optim(beta.initial, 
                    log.poisson.posterior.of.beta,
                    gr      = NULL,
                    hessian = TRUE,
                    method  = c("BFGS"),
                    control = list(fnscale=-1))

# beta.prior <- dnorm(grid, 0, 100* solve(t(C) %*% C)
# Ta bort / ha med intercept?

# c)
# Assumes diagonal covariance matrix, is that fine?
metropolis <- function(c, logPostFunc, theta, ...) {
  n       <- length(theta)
  sigma   <- c*diag(n)
  theta_p <- as.vector(rmvnorm(1, theta, sigma, ...))
  log.prob.theta   <- logPostFunc(theta, ...)
  log.prob.theta_p <- logPostFunc(theta_p, ...)
  
  alpha  <- min(1, exp(log.prob.theta_p - log.prob.theta))
  accept <- runif(1, 0, 1) > alpha
  if (accept) {
    list(theta = theta_p, accept = accept)
  } else {
    list(theta = theta, accept = accept)
  }
}

nSamples <- 3000
d <- dim(X)[2]
c <- 0.0000001
betas     <- matrix(rep(rep(0, d), nSamples), ncol = d)
betas[1,] <- rmvnorm(1, optim.beta$par, 0.1*diag(d))
accepts <- rep(0, nSamples-1)
log.approx.normal.posterior <- function (x) { 
  dmvnorm(x, optim.beta$par, -solve(optim.beta$hessian), log = TRUE) 
}

for (i in 2:nSamples) {
  result <- metropolis(c, log.approx.normal.posterior, betas[i-1,])
  betas[i,]  <- result$theta
  accepts[i-1] <- result$accept
}
burnin <- nSamples / 2
accepts.reduced <- accepts[burnin:nSamples-1]
accept.ratio <- sum(accepts.reduced) / length(accepts.reduced)
approx.post.mean <- apply(betas[burnin:nSamples,], 2, mean)

accept.ratio
approx.post.mean

plot(optim.beta$par, col ="blue", type ="l", ylim = c(-2.5,1.5))
points(ebay.glm$coefficients[2:9], col ="red",type ="l")
points(approx.post.mean, type="l", col="green")

# These shuld be noise, right? Too bad.
plot(betas[,1], type = "l", col="red")

plot(betas[,2], type = "l", col="red")
plot(betas[,3], type = "l", col="red")
plot(betas[,4], type = "l", col="blue")
plot(betas[,5], type = "l", col="blue")
plot(betas[,6], type = "l", col="blue")
plot(betas[,7], type = "l", col="blue")
plot(betas[,8], type = "l", col="blue")

