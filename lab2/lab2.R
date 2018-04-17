library(mvtnorm)
data <- read.csv("TempLinkoping.txt", sep = "\t")
X <- data$time
Y <- data$temp

# 1. Linear and polynomial regression

rinvchisq <- function(vx, sigmax, draws) {
  vx*sigmax/rchisq(draws, vx)
}

#a) Determining the prior distribution
delta <- 0.01
times <- seq(delta, 1, delta)
nDraws  <- 5000 
posteriorDraws <- 5000
mu0    <- c(-10, 150, -150)   # From tuning the plotted polynomial
omega0 <- 1*diag(3)           # same variance for betas, we dont know much
v0      <- 20
sigma0  <- 2
sigma2s <- rinvchisq(v0, sigma0, nDraws)#v0*sigma0/rchisq(nDraws, v0)
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
prior.temp <- apply(thetas[seq(1, 10),], 1, function(theta) {
  sapply(times, function(time) {
    theta[1] + theta[2]*time + theta[3]*time^2 + rnorm(1, 0, theta[4])
  })
})


#b) Check if prior distribution is sensible
plot(prior.temp[,1], type = 'l', ylim = c(-15, 40))
for (i in 2:10) {
  lines(prior.temp[,i])
}

# c) Simulating from posterior distribution
A   <- as.numeric(t(X) %*% X)
B   <- as.numeric(t(Y) %*% Y)
beta_hat <- apply(thetas[, 1:3], 2, mean) # minimize sqaure loss
mun      <- solve(A + omega0) %*% (A*beta_hat + omega0%*%mu0)
omegan   <- A + omega0
vn       <- v0 + n
sigma2n  <- as.numeric((v0%*%sigma0 + (B + t(mu0)%*%omega0%*%mu0 - t(mun)%*%omegan%*%mun))/vn)
# Super big negative values gives negative variance in posterior beta distrubution ^- 


post.sigma2s <- rinvchisq(vn, sigma2n, posteriorDraws)
post.betas   <- sapply(post.sigma2s, function(post.sigma2) {
  rmvnorm(1, mun, post.sigma2 * solve(omegan)) 
})




# Posterior approximation for classification with logistic regression

# a)
data <- read.table("WomenWork.dat", header=TRUE)
y    <- as.vector(data[,1])
X    <- as.matrix(data[,-1])

glmModel <- glm(Work ~0 + ., data=data, family = binomial)

# b)
# want posterior
# to get that first we need prior and likelihood

beta.log.posterior <- function(betas, X, y, mu, sigma) {
  d    <- length(betas)
  pred <- X %*% betas
  
  log.likelihood <- sum(y*pnorm(pred, log.p = TRUE) + (1-y)*pnorm(pred, log.p = TRUE, lower.tail = FALSE))
  log.prior      <- dmvnorm(betas, matrix(mu, d, 1), sigma*diag(d), log = TRUE)
  
  log.likelihood + log.prior
}

tau2 <- 10^2
mu   <- 0
d <- dim(X)[2]
beta.initial <- rep(0, d)
result       <- optim(beta.initial, beta.log.posterior, gr = NULL, 
                      X, y, mu, tau2, 
                      hessian = TRUE,
                      method  = c("BFGS"),
                      control = list(fnscale=-1))

beta_hat <- result$par
J        <- result$hessian

# numerical values should be in the report
beta.post <- rmvnorm(10000, beta_hat, -solve(J))

# compute 95% credible interval for NSmallChild?
x <- c(1, 10, 8, 10, 1, 40, 1, 1)
y_hat <- beta.post %*% x
p <- 1 / (1 + exp(-y_hat))
hist(p, breaks = 100)
