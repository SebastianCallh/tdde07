
# 1 - Poisson regression - the MCMC way

# a)
ebay.data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
ebay.glm <- glm(nBids ~ .-Const, data = ebay.data, family = poisson)
summary(ebay.glm)

# Estimate tells us: MinBidShare most important -1.89, sealed 0.44, veriID 0.39, majBlem -0.22 ,logBook -0.12, bleh

# b)
library(optimr)
library(mvtnorm)

n <- dim(ebay.data)[1]
p <- dim(ebay.data)[2]
X <- scale(as.matrix(ebay.data[,3:10]))
Y <- ebay.data[,1]
C <- cov(X)
lambda <- 5
sigma  <- 100*solve(t(C)%*%C)

log.poisson.posterior.of.beta <- function(beta) { 
  print(length(beta))
  print(dim(sigma))
  lambda <- X %*% beta
  beta.log.prior         <- dmvnorm(beta, sigma = sigma, log = TRUE)
  poisson.log.likelihood <- sum(lambda^Y*exp(-lambda)) #dpois(Y, z, log = TRUE))
  poisson.log.likelihood + beta.log.prior
}

beta.initial = rep(0, dim(X)[2])
log.poisson.posterior.of.beta(c(1, 2, 1, 2, 1, 1, 1, 1))

## This gives values around 10^280
optim.beta <- optim(beta.initial, 
                    log.poisson.posterior.of.beta,
                    gr      = NULL,
                    hessian = TRUE,
                    method  = c("BFGS"),
                    control = list(fnscale=-1))

# beta.prior <- dnorm(grid, 0, 100* solve(t(C) %*% C)

