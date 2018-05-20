
# 1 - Poisson regression - the MCMC way

# a)
ebay.data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)
ebay.glm <- glm(nBids ~ .-Const, data = ebay.data, family = poisson(link=log))
summary(ebay.glm)

# Estimate tells us: MinBidShare most important -1.89, sealed 0.44, veriID 0.39, majBlem -0.22 ,logBook -0.12, bleh

# b)
library(optimr)
library(mvtnorm)
library(ggplot2)
library(gridExtra)

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
metropolis.sample <- function(c, sigma, logPostFunc, theta_c, ...) {
  theta_p <- as.vector(rmvnorm(1, theta_c, c*sigma))
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
  sigma.proposal <- -solve(optim.beta$hessian)
  betas   <- matrix(rep(rep(0, p), nSamples), ncol = p)
  accepts <- rep(0, nSamples-1)
  for (i in 2:nSamples) {
    result <- metropolis.sample(c, sigma.proposal, 
                         log.posterior.of.beta, 
                         betas[i-1,],
                         mu0,
                         sigma0,
                         X,
                         Y)
    betas[i,]    <- result$theta
    accepts[i-1] <- result$accept
  }
  
  valid.accepts <- accepts[burnin:nSamples]
  accept.ratio <- sum(valid.accepts) / length(valid.accepts)
  list(betas = betas, accept.ratio = accept.ratio)
}

nSamples <- 5000
burnin   <- 1000
res <- metropolis(nSamples, burnin)
betas <- res$betas
beta.post.mean <- apply(betas, 2, mean)

plot.graphical.comparison <- function () {
  plot(optim.beta$par, col ="gray", 
       ylim = c(-2.5,1.5), 
       xlab = expression(beta), 
       ylab = "Value",
       main = "Graphical comparison of the differently estimated betas")
  points(ebay.glm$coefficients, col ="green")
  points(beta.post.mean, col="orange")
  legend("topright",
         c("Point estimate", "Normal approximation", "Mean of MCMC samples"),
                   lty=1,
                   col=c("gray", "green", "orange"),
                   bty='n',
                   cex=.75)
}

plot.posterior.beta.samples <- function () {
  x.axis <- seq(1, dim(betas)[1], 1)
  beta.df <- data.frame(betas, x.axis)
  pb1 <- ggplot(beta.df, mapping = aes(x.axis, X1)) + geom_line() + labs(x = expression(beta[1]), y = "Value")
  pb2 <- ggplot(beta.df, mapping = aes(x.axis, X2)) + geom_line() + labs(x = expression(beta[2]), y = "Value")
  pb3 <- ggplot(beta.df, mapping = aes(x.axis, X3)) + geom_line() + labs(x = expression(beta[3]), y = "Value")
  pb4 <- ggplot(beta.df, mapping = aes(x.axis, X4)) + geom_line() + labs(x = expression(beta[4]), y = "Value")
  pb5 <- ggplot(beta.df, mapping = aes(x.axis, X5)) + geom_line() + labs(x = expression(beta[5]), y = "Value")
  pb6 <- ggplot(beta.df, mapping = aes(x.axis, X6)) + geom_line() + labs(x = expression(beta[6]), y = "Value")
  pb7 <- ggplot(beta.df, mapping = aes(x.axis, X7)) + geom_line() + labs(x = expression(beta[7]), y = "Value")
  pb8 <- ggplot(beta.df, mapping = aes(x.axis, X8)) + geom_line() + labs(x = expression(beta[8]), y = "Value")
  grid.arrange(pb1, pb2, pb3, pb4, pb5, pb6, pb7, pb8, ncol = 2, top ="Beta posterior samples")
}

plot.posterior.phis <- function () {
  phis <- exp(betas)
  phi.df <- data.frame(phis)
  pp1 <- ggplot(phi.df, aes(x = X1)) + geom_density() + labs(x = expression(phi[1]), y = "Density")
  pp2 <- ggplot(phi.df, aes(x = X2)) + geom_density() + labs(x = expression(phi[2]), y = "Density")
  pp3 <- ggplot(phi.df, aes(x = X3)) + geom_density() + labs(x = expression(phi[3]), y = "Density")
  pp4 <- ggplot(phi.df, aes(x = X4)) + geom_density() + labs(x = expression(phi[4]), y = "Density")
  pp5 <- ggplot(phi.df, aes(x = X5)) + geom_density() + labs(x = expression(phi[5]), y = "Density")
  pp6 <- ggplot(phi.df, aes(x = X6)) + geom_density() + labs(x = expression(phi[6]), y = "Density")
  pp7 <- ggplot(phi.df, aes(x = X7)) + geom_density() + labs(x = expression(phi[7]), y = "Density")
  pp8 <- ggplot(phi.df, aes(x = X8)) + geom_density() + labs(x = expression(phi[8]), y = "Density")
  grid.arrange(pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, ncol = 2, top = "Phi posterior densities")
}

# d)

# New data point with const
x <- c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)
beta <- metropolis(5000, 1000)$betas
posterior.prob <- dpois(0, exp(x %*% t(beta)))
plot.posterior.predictive.density <- function () {
  df <- data.frame(t(posterior.prob))
  ggplot(df, aes(x = t.posterior.prob.)) +
    geom_density() +
    xlab("Probability") +
    ylab("Density") +
    ggtitle("Plot of posterior probability of zero bidders")
}
mean(posterior.prob)
