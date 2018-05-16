
# 1 - Poisson regression - the MCMC way

# a)
ebay.data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
ebay.glm <- glm(nBids ~ .-Const, data = ebay.data, family = poisson(link=log))
summary(ebay.glm)

# Estimate tells us: MinBidShare most important -1.89, sealed 0.44, veriID 0.39, majBlem -0.22 ,logBook -0.12, bleh

# b)
library(optimr)
library(mvtnorm)

log.posterior.of.beta <- function(beta, mu0, sigma0, X, Y) { 
  beta.log.prior         <- dmvnorm(beta, mu0, sigma0, log = TRUE)
  poisson.log.likelihood <- sum(dpois(Y, exp(X %*% beta), log = TRUE))
  poisson.log.likelihood + beta.log.prior
}

Y      <- ebay.data[,1]
X      <- as.matrix(ebay.data[,2:10])
n      <- dim(X)[1]
p      <- dim(X)[2]
mu0       <- rep(0, p)
sigma0    <- 100*solve(t(X)%*%X)
beta.init <- rep(0, p)

optim.beta <- optim(beta.init, 
                    log.posterior.of.beta,
                    gr      = NULL,
                    mu      = mu0,
                    sigma   = sigma0,
                    X       = X,
                    Y       = Y,
                    hessian = TRUE,
                    method  = c("BFGS"),
                    control = list(fnscale=-1))

log.approx.normal.posterior <- function (x) { 
  dmvnorm(x, optim.beta$par, -solve(optim.beta$hessian), log = TRUE) 
}

# c)
metropolis.sample <- function(sigma, logPostFunc, theta_c, ...) {
  theta_p <- as.vector(rmvnorm(1, theta_c, sigma))
  log.prob.theta_c <- logPostFunc(theta_c, ...)
  log.prob.theta_p <- logPostFunc(theta_p, ...)
  alpha  <- min(1, exp(log.prob.theta_p - log.prob.theta_c))
  accept <- runif(1, 0, 1) < alpha
  
  if (accept) {
    list(theta = theta_p, accept = accept)
  } else {
    list(theta = theta_c, accept = accept)
  }
}

# Initial beta is zero vector
metropolis <- function (nSamples, burnin) {
  c <- 0.5
  sigma.proposal <- -c * solve(optim.beta$hessian)
  betas     <- matrix(rep(rep(0, p), nSamples), ncol = p)
  accepts <- rep(0, nSamples-1)
  for (i in 2:nSamples) {
    result <- metropolis.sample(sigma.proposal, 
                         log.posterior.of.beta, 
                         betas[i-1,],
                         mu0,
                         sigma0,
                         X,
                         Y)
    betas[i,]  <- result$theta
    accepts[i-1] <- result$accept
  }
  
  valid.accepts <- accepts[burnin:nSamples]
  accept.ratio <- sum(valid.accepts) / length(valid.accepts)
  valid.betas <- betas[burnin:nSamples,]
  list(betas = valid.betas, accept.ratio = accept.ratio)
}
res <- metropolis(5000, 1000)

plot(optim.beta$par, col ="blue", type ="l", ylim = c(-2.5,1.5))
points(ebay.glm$coefficients[2:9], col ="red",type ="l")
points(approx.post.mean, type="l", col="green")

plot(valid.betas[,1], type = "l", col="red")
plot(valid.betas[,2], type = "l", col="red")
plot(valid.betas[,3], type = "l", col="red")
plot(valid.betas[,4], type = "l", col="blue")
plot(valid.betas[,5], type = "l", col="blue")
plot(valid.betas[,6], type = "l", col="blue")
plot(valid.betas[,7], type = "l", col="blue")
plot(valid.betas[,8], type = "l", col="blue")

# d)

# New data point with const
x <- c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)
beta <- metropolis(5000, 1000)$betas
prob <- dpois(0, exp(x %*% t(beta)))
plot(density(prob))
mean(prob)     
